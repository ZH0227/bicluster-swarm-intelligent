import pandas as pd 
import numpy as np 
import os
from glob import glob

def _calc(df):
    """
    WEScore = ∑_{i=1}^t{x_i s_i}/k
    meanPValue(B) = ∑_{i=1} ^{t} s_i ⁄ t
    rateGeneTerm(B) = k ⁄ t
    """
    t = df.shape[0]
    k = int(df.iloc[0]['GeneRatio'].split('/')[-1])
    S = -np.log10(df['p.adjust']).to_numpy()
    X = df['Count'].to_numpy()

    return np.sum(S * X) / k, np.sum(S) / t, k / t

if __name__ == "__main__":
    p_threshold = 0.0005
    datasets = ['BCLL', 'PBC', 'YC', 'RAT']
    algs = ['FA', 'CS', 'CSFA', 'PSO', 'QPSO']
    for dataset in datasets:
        weDict = {}
        meanPDict = {}
        rateGTDict = {}
        dataPath = os.path.join('./', dataset)
        print(dataPath)

        for alg in algs:
            meanPList = []
            weList = []
            rateGTList = []

            algPath = os.path.join(dataPath, alg+'_go')
            csvs = glob(os.path.join(algPath, "*.csv"))
            for csv in csvs:
                df = pd.read_csv(csv, index_col=0)
                # df = df[df['p.adjust']<p_threshold]
                we, meanP, rateGT = _calc(df)
                weList.append(we)
                meanPList.append(meanP)
                rateGTList.append(rateGT)
            
            weDict[alg+'B'] = weList
            meanPDict[alg+'B'] = meanPList
            rateGTDict[alg+'B'] = rateGTList

        we_df = pd.DataFrame(weDict)
        we_csv = os.path.join(dataPath, 'we.csv')
        we_df.to_csv(we_csv)

        meanP_df = pd.DataFrame(meanPDict)
        meanP_csv = os.path.join(dataPath, 'meanP.csv')
        meanP_df.to_csv(meanP_csv)

        rateGT_df = pd.DataFrame(rateGTDict)
        rateGT_csv = os.path.join(dataPath, 'rateGT.csv')
        rateGT_df.to_csv(rateGT_csv)
    