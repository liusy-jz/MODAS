import numpy as np
import pandas as pd
from sklearn.decomposition import MiniBatchSparsePCA
from sklearn.metrics import pairwise_distances
from sklearn.linear_model import LinearRegression
from scipy.stats import t
import statsmodels.api as sm
from statsmodels.formula.api import ols


def standardize_phe(phe):
    standardized_phe = (phe - np.mean(phe, axis=0)) / np.std(phe, axis=0)
    return np.nan_to_num(standardized_phe)


def phe_transform(phe):
    transformed_phe = pd.DataFrame(pairwise_distances(phe + 1, metric='braycurtis'), index=phe.index, columns=phe.index)
    transformed_phe[np.isnan(transformed_phe)] = 0
    return transformed_phe


def scpca(stress, control, phe_name, alpha, n_comp):
    index = stress.index
    stress = stress - np.mean(stress, axis=0)
    stress = standardize_phe(stress)
    control = control - np.mean(control, axis=0)
    control = standardize_phe(control)
    stress_cov = stress.T.dot(stress) / (stress.shape[0] - 1)
    control_cov = control.T.dot(control) / (control.shape[0] - 1)
    spca = MiniBatchSparsePCA(n_components=n_comp, random_state=3)
    scpca_res = spca.fit_transform(stress_cov - alpha * control_cov)
    scpca_res = pd.DataFrame(scpca_res, index=index, columns=[phe_name+'_CPCA_PC'+str(i) for i in range(1, n_comp + 1)])
    return scpca_res


def beta_test(scpca_pc, stress_phe, control_phe):
    lr = LinearRegression()
    lr.fit(scpca_pc.to_frame(), stress_phe)
    stress_beta = lr.coef_[0]
    rss_stress = ((stress_phe - lr.predict(scpca_pc.to_frame()))**2).sum()
    lr.fit(scpca_pc.to_frame(), control_phe)
    control_beta = lr.coef_[0]
    rss_control = ((control_phe - lr.predict(scpca_pc.to_frame()))**2).sum()
    s_scpca_phe = (rss_stress + rss_control) / ((stress_phe.shape[0] - 2) + (control_phe.shape[0] - 2))
    s_con = np.sqrt(s_scpca_phe / ((scpca_pc - scpca_pc.mean())**2).sum() + s_scpca_phe / ((scpca_pc - scpca_pc.mean())**2).sum())
    t_stat = (stress_beta - control_beta) / s_con
    pvalue = (1 - t.cdf(abs(t_stat), (stress_phe.shape[0] - 2) + (control_phe.shape[0] - 2))) * 2
    return [stress_beta, control_beta, t_stat, pvalue]


def scpca_omics(stress_omics, control_omics, alpha, n_comp, beta_test_pvalue):
    idx = list(set(stress_omics.columns.tolist()) & set(control_omics.columns.tolist()))
    index = stress_omics.index.intersection(control_omics.index)
    stress_omics = stress_omics.reindex(index)
    control_omics = control_omics.reindex(index)
    omics_scpca_phe = list()
    omics_scpca_test = list()
    for mol_feature in idx:
        stress = phe_transform(stress_omics[[mol_feature]])
        control = phe_transform(control_omics[[mol_feature]])
        scpca_res = scpca(stress, control, mol_feature, alpha, n_comp)
        omics_scpca_phe.append(scpca_res)
        for scpca_phe_id in scpca_res:
            omics_scpca_test.append([scpca_phe_id, mol_feature] + beta_test(scpca_res[scpca_phe_id], stress_omics[mol_feature], control_omics[mol_feature]))
    omics_scpca_phe = pd.concat(omics_scpca_phe, axis=1)
    omics_scpca_test = pd.DataFrame(omics_scpca_test)
    omics_scpca_test.columns = ['scpca_phe', 'molecular_feature', 'stress_beta', 'control_beta', 'beta_test_stat', 'beta_test_pvalue']
    if beta_test_pvalue is None:
        beta_test_pvalue = 1.0 / omics_scpca_test.shape[0]
    omics_scpca_phe = omics_scpca_phe.loc[:, omics_scpca_test.loc[omics_scpca_test.beta_test_pvalue <= beta_test_pvalue, 'scpca_phe']]
    return omics_scpca_phe, omics_scpca_test


def anova(stress_omics, control_omics, scpca_qtl_res, geno):
    index = stress_omics.index.intersection(geno.index)
    stress_omics = stress_omics.reindex(index)
    control_omics = control_omics.reindex(index)
    geno = geno.reindex(index)
    anova_res = list()
    for phe_name in scpca_qtl_res.phe_name.unique():
        mol_feature = '_'.join(phe_name.split('_')[:-2]).replace('m.z', 'm/z')
        res = pd.concat([stress_omics[mol_feature].to_frame().reset_index(), control_omics[mol_feature].to_frame().reset_index()])
        res.columns = ['RIL', 'exp']
        res['con'] = ['stress'] * stress_omics.shape[0] + ['control'] * control_omics.shape[0]
        for snp in scpca_qtl_res.loc[scpca_qtl_res.phe_name == phe_name, 'SNP']:
            res['genotype'] = 0
            res.loc[res['RIL'].isin(geno.index[geno[snp] == 2]), 'genotype'] = 1
            model = ols('exp ~ C(con) + C(genotype) +\
                        C(con):C(genotype)', data=res).fit()
            pvalue = sm.stats.anova_lm(model, typ=2).iloc[2, 3]
            anova_res.append([phe_name, snp, pvalue])
    anova_res = pd.DataFrame(anova_res)
    anova_res.columns = ['phe_name', 'SNP', 'pvalue']
    scpca_qtl_res = pd.merge(scpca_qtl_res, anova_res, on=['phe_name', 'SNP'])
    return scpca_qtl_res
