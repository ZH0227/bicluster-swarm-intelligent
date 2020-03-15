import numpy as np 
import pandas as pd 
import os
from glob import glob

datasets = ['BCLL', 'YC', 'PBC', 'RAT']
algs = ['FA', 'CS', 'CSFA', 'PSO', 'QPSO']
metrices = ['GV', 'CV', 'MSR', 'Var']
def helper(dataName):
    cs_df = pd.read_csv(os.path.join('./',dataName, "CS.csv"))
    csfa_df = pd.read_csv(os.path.join('./',dataName, "CSFA.csv"))
    fa_df = pd.read_csv(os.path.join('./',dataName, "FA.csv"))
    pso_df = pd.read_csv(os.path.join('./',dataName, "PSO.csv"))
    qpso_df = pd.read_csv(os.path.join('./',dataName, "QPSO.csv"))

    msr_df = pd.DataFrame()
    cv_df = pd.DataFrame()
    gv_df = pd.DataFrame()
    var_df = pd.DataFrame()

    for alg in algs:
        cv_df[alg+'B'] = eval(alg.lower()+'_df')['CV'].to_numpy()
        gv_df[alg+'B'] = eval(alg.lower()+'_df')['GV'].to_numpy()
        msr_df[alg+'B'] = eval(alg.lower()+'_df')['MSR'].to_numpy()
        var_df[alg+'B'] = eval(alg.lower()+'_df')['Var'].to_numpy()

    for met in metrices:
        eval(met.lower()+"_df").to_csv(os.path.join('./',dataName, met.lower()+'.csv'))

if __name__ == "__main__":
    for dataName in datasets:
        helper(dataName)