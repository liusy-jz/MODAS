#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: setup.py
Author: CrazyHsu @ crazyhsu9527@gmail.com
Created on: 2020-09-25 15:46:07
Last modified: 2020-09-25 16:06:19
'''

from setuptools import setup

setup(name='MODAS',
    python_requires='>3.7',
    version='1.0',
    description='MODAS: Multi-omics data association study',
    url='https://github.com/CrazyHsu/MODAS.git',
    author='Songyu Liu',
    author_email='',
    license='GPLv3',
    packages=['modas'],
    scripts=['MODAS.py'],
    install_requires=[
        'cython',
        'yattag',
        'pandas_plink',
        'numpy',
        'matplotlib',
        'pandas',
        'scipy',
        'scikit-learn',
    ]
)
