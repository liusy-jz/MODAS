#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: setup.py
Author: CrazyHsu @ crazyhsu9527@gmail.com 
Created on: 2020-09-25 15:46:07
Last modified: 2022-01-08 02:00:57
'''

from setuptools import setup

setup(name='MODAS',
    python_requires='>3.7',
    version='1.0',
    description='MODAS: Multi-omics data association study',
    url='https://github.com/liusy-jz/MODAS.git',
    author='Songyu Liu',
    author_email='',
    license='GPLv3',
    packages=['modas'],
    scripts=['MODAS.py'],
    install_requires=[
        'cython',
        'yattag==1.14.0',
        'matplotlib',
        'numpy==1.21.0',
        'pandas==1.3.5',
        'scipy==1.7.3',
        'scikit-learn==1.0.2',
        'pandas_plink'
    ]
)
