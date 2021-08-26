#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Implements functions related to calculating and handeling cortical maps 

@author: jdiedrichsen
"""

import os # to handle path information
import nibabel as nb
import numpy as np
import h5py
import pandas as pd
import surfAnalysisPy as surf

return_subjs = np.array([2,3,4,6,8,9,10,12,14,15,17,18,19,20,21,22,24,25,26,27,28,29,30,31])


baseDir = '/Volumes/diedrichsen_data$/data/super_cerebellum'
surfDir = baseDir + '/sc1/surfaceWB'

hem = ['L','R']
hem_name = ['CortexLeft','CortexRight']

def make_wcon(sn = return_subjs, glm = 7, exp = 'sc1',sess = [1,2]):
    glm_str = f"glm{glm}"
    for s in sn: 
        subj = f"s{s:02}"
        print(subj)
        subjDir = os.path.join(surfDir,glm_str,subj)
        infoFile = os.path.join(baseDir,exp,f"GLM_firstlevel_{glm}",subj,'SPM_info.mat')
        T = read_SPMinfo(infoFile)
        num_cond = np.max(T.cond)
        for i,h in enumerate(hem):
            # Load beta file 
            inp_name = f"{subj}.{h}.{exp}.beta.func.gii"
            S = nb.load(os.path.join(subjDir,inp_name))
            X = np.vstack(S.agg_data())
            X = X[:T.cond.shape[0],:]
            # Get residual mean square file 
            inp_name = f"{subj}.{h}.{exp}.ResMS.func.gii"
            S = nb.load(os.path.join(subjDir,inp_name))
            R = S.agg_data()

            data = np.zeros((X.shape[1],num_cond))
            coln = []
            for c in range(num_cond): 
                con = np.logical_and(T.cond==c+1,np.isin(T.sess,sess))
                data[:,c] = (con @ X)/np.sqrt(R)
                a, =np.where(T.cond==c+1)
                coln.append(T.CN[a[0]])
            G = surf.map.make_func_gifti(data,hem_name[i],coln)
            if (type(sess) is int):
                out_name = f"{subj}.{h}.{exp}.sess{sess}.wcon.func.gii"
            else:
                out_name = f"{subj}.{h}.{exp}.wcon.func.gii"
            nb.save(G,os.path.join(subjDir,out_name))

def read_SPMinfo(filename): 
    INFO = h5py.File(filename)
    dict = {}
    for name,data in INFO.items():
        # print("Name " + name) # Name
        if type(data) is h5py.Dataset:
            # If DataSet pull the associated Data
            # If not a dataset, you may need to access the element sub-items
            if data.dtype == 'O': # Object assume string... 
                dict[name] = _convertobj(file_obj=INFO, key=name)
                pass
            elif data.dtype == '<f8':
                dict[name] = np.array(data).reshape(-1).astype(int)
            else:
                print(f'unknown data type {data.dtype}')
    return pd.DataFrame(dict)

def _convertobj(file_obj, key):
    """converts object reference for `key` in `file_obj`"""
    dataset = file_obj[key]
    tostring = lambda obj: "".join(chr(i[0]) for i in obj[:])
    return [tostring(file_obj[val]) for val in dataset[()].flatten()]

if __name__ == "__main__":
    make_wcon(exp='sc1',sess=[1,2])
    make_wcon(exp='sc1',sess=1)
    make_wcon(exp='sc1',sess=2)
    make_wcon(exp='sc2',sess=[1,2])
    make_wcon(exp='sc2',sess=1)
    make_wcon(exp='sc2',sess=2)