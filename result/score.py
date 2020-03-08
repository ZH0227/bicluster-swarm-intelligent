import pandas as pd 
import numpy as np 
import os
from glob import glob
import pickle

def calc_hv(bic_np):
    bIj = np.mean(bic_np, 0, keepdims=True)
    biJ = np.mean(bic_np, 1, keepdims=True)
    bIJ = np.mean(bic_np)

    fenzi = np.sum((bic_np-bIj-biJ+bIJ)**2)
    fenmu = np.sum((bic_np-biJ)**2)
    return fenzi / fenmu

def calc_ri(bic_np, alldata_np):
    sigma_Ij = np.std(bic_np, 0, keepdims=True)
    sigma_j = np.std(alldata_np, 0, keepdims=True)

    return np.sum(1 - sigma_Ij / sigma_j) / bic_np.shape[1]

def calc_smsr(bic_np):
    bIj = np.mean(bic_np, 0, keepdims=True)
    biJ = np.mean(bic_np, 1, keepdims=True)
    bIJ = np.mean(bic_np)
    vol = bic_np.shape[0] * bic_np.shape[1]

    return np.sum(((biJ * bIj - bic_np * bIJ) / biJ * bIj)**2) /vol

def read_txt(txtFile):
    if not os.path.exists(txtFile):
        return
    with open(txtfile, 'r') as f:
        lines = [line.strip() for line in f]

    genes_str = lines[4::6]
    conds_str = lines[5::6]
    genes = [list(map(int, line.split())) for line in genes_str]
    conds = [list(map(int, line.split())) for line in conds_str]
    return genes, conds

def saveScore(filePath, scores):
    with open(filePath, 'wb') as f:
        pickle.dump(scores, f)

    #double check
    if os.path.getsize(filePath) <= 0:  
        with open(filePath, 'wb') as f:
            pickle.dump(scores, f)


def processHV(txtFile, alldata_np):
    if not os.path.exists(txtFile):
        return
    pklFile = txtFile.replace('.txt', '_hvScores.pkl')
    if os.path.exists(pklFile) and os.path.getsize(pklFile) >0:
        with open(pklFile, 'rb') as f:
            return pickle.load(f)
    genes, conds = read_txt(txtFile)

    scores = []
    for i in range(len(genes)):
        bic = alldata_np[genes[i],:][:,conds[i]]
        scores.append(calc_hv(bic))
    # saveScore(pklFile, scores)
    return scores

def processRI(txtFile, alldata_np):
    if not os.path.exists(txtFile):
        return
    pklFile = txtFile.replace('.txt', '_riScores.pkl')
    if os.path.exists(pklFile) and os.path.getsize(pklFile) >0:
        with open(pklFile, 'rb') as f:
            return pickle.load(f)
    genes, conds = read_txt(txtFile)

    scores = []
    for i in range(len(genes)):
        bic = alldata_np[genes[i],:][:,conds[i]]
        scores.append(calc_ri(bic, alldata_np[:,conds[i]]))
    #saveScore(pklFile, scores)
    return scores

def processSMSR(txtFile, alldata_np):
    if not os.path.exists(txtFile):
        return
    pklFile = txtFile.replace('.txt', '_smsrScores.pkl')
    if os.path.exists(pklFile) and os.path.getsize(pklFile) >0:
        with open(pklFile, 'rb') as f:
            return pickle.load(f)
    genes, conds = read_txt(txtFile)

    scores = []
    for i in range(len(genes)):
        bic = alldata_np[genes[i],:][:,conds[i]]
        scores.append(calc_smsr(bic))
    # saveScore(pklFile, scores)
    return scores

if __name__ == "__main__":
    datasets = ['BCLL', 'PBC', 'YC', 'RAT']
    algs = ['FA', 'CS', 'CSFA', 'PSO', 'QPSO']
    for dataset in datasets:
        datafile = glob('../../data/*'+dataset+'_norm.csv')[0]
        alldata_np = pd.read_csv(datafile, index_col=0).values
        print(datafile, alldata_np.shape)
        hvDict = {}
        riDict = {}
        smsrDict = {}
        for alg in algs:
            txtfile = os.path.join(dataset, alg+'.txt')
            print(txtfile)
            hvDict[alg] = processHV(txtfile, alldata_np)
            riDict[alg] = processRI(txtfile, alldata_np)
            smsrDict[alg] = processSMSR(txtfile, alldata_np)
        hvfile = os.path.join('./', dataset, "hv.csv")
        hv_df = pd.DataFrame(hvDict)
        hv_df.to_csv(hvfile)

        rifile = os.path.join('./', dataset, "ri.csv")
        ri_df = pd.DataFrame(riDict)
        ri_df.to_csv(rifile)
        
        smsrfile = os.path.join('./', dataset, "smsr.csv")
        smsr_df = pd.DataFrame(smsrDict)
        smsr_df.to_csv(smsrfile)

    