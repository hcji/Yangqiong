# -*- coding: utf-8 -*-
"""
Created on Tue Nov  6 11:25:26 2018

@author: hcji
"""

import rpy2.robjects as robjects
import rpy2.robjects.numpy2ri as numpy2ri
import numpy as np

def init_r():
    """
        1. Install R >= 3.4.1
        2. Install Java SE Runtime Environment 8
        3. Install rcdk
    """
    numpy2ri.activate()
    robjects.r('''source('cdk/cdk.R')''')
    cdk_fingerprint = robjects.globalenv['cdk_fingerprint']
    return cdk_fingerprint


def get_cdk_fingerprints(smi, tp=['pubchem', 'kr']):
    '''
    supported:
    'pubchem', 'maccs', 'standard', 'extended', 'graph', 'hybridization', 'estate', 'kr', 'shortestpath', 'signature', 'circular'
    '''
    cdk_fingerprint = init_r()
    res = np.array(cdk_fingerprint(smi, tp))
    res = res.astype(int)
    return res
    

if __name__ == '__main__':
    fp = get_cdk_fingerprints('CCCN', tp=['pubchem', 'kr'])
    