#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  1 15:37:30 2022

@author: ken
"""

from scipy.stats import skew, kurtosis
import pandas as pd
import numpy as np

def bimodality_coefficient(df):
    s_arr = [0] * df.shape[1]
    k_arr = [0] * df.shape[1]
    n_arr = [0] * df.shape[1]
    tmp_arr = [0] * df.shape[1]
    bc_arr = [0] * df.shape[1]
    for i in range(len(df.columns)):
        tmp_df = pd.DataFrame(df.iloc[:,i].dropna(how='any', axis=0))
        
        n = tmp_df.shape[0]
        s = skew(tmp_df)[0]
        k = kurtosis(tmp_df, bias=False)[0]
        tmp = (n-1)*(n-1)/((n-2)*(n-3))
        bc = (s*s+1) / (k + 3 * tmp)
        
        s_arr[i] = s
        k_arr[i] = k
        n_arr[i] = n
        tmp_arr[i] = tmp
        bc_arr[i] = bc
    return s_arr, k_arr, n_arr, tmp_arr, bc_arr

def bimodality_coefficient_difference(grp1, grp2):
    result = pd.DataFrame(None, index=grp1.columns)
    
    result['case_s'], result['case_k'], result['case_n'], result['case_n_terms'], result['case_bc'] = bimodality_coefficient(grp1)
    result['ctrl_s'], result['ctrl_k'], result['ctrl_n'], result['ctrl_n_terms'], result['ctrl_bc'] = bimodality_coefficient(grp2)
   

    assert len(result['case_bc']) == len(result['ctrl_bc'])
    result['bcd'] = [0] * len(result['case_bc'])
    for i in result.index:
        result.loc[i,'bcd'] = abs(result.loc[i, 'case_bc'] - result.loc[i,'ctrl_bc'])
    
    return result

def clipIQR(df, factor):
    for j in range(len(df.columns)):
        col_data = df.iloc[:,j].dropna(how='any', axis=0)
        q3, q1 = np.percentile(col_data, [75, 25])
        iqr = q3 - q1
        for i in range(df.shape[0]):
            if df.iloc[i,j] > q3 + factor * iqr:
                df.iloc[i,j] = q3
            elif df.iloc[i,j] < q1 - factor * iqr:
               df.iloc[i,j] = q1        
    return df


base_dir = '/media/ken/ExtraSpace/School/Research/DataSets/Myers_expr/'
case_file = base_dir + 'AD_noBrain2_4_cases.tsv'
ctrl_file = base_dir + 'AD_noBrain2_4_controls.tsv'

case_df = pd.read_csv(case_file, sep='\t', index_col=0)
ctrl_df = pd.read_csv(ctrl_file, sep='\t', index_col=0)

case_df = case_df.transpose()
ctrl_df = ctrl_df.transpose()

all_df = pd.concat([case_df,ctrl_df])

all_df = clipIQR(all_df, 3)
    
case_df = all_df.loc[case_df.index,:]
ctrl_df = all_df.loc[ctrl_df.index,:]

ad_bcd = bimodality_coefficient_difference(case_df, ctrl_df)
