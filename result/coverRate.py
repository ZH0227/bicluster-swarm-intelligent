from glob import glob
import os
import pandas as pd 

def parse(txtFile:str):
    """
    取基因和条件的集合
    intput: txtFile:str
    return: genes, conds
    """
    if not os.path.exists(txtFile):
        return
    with open(txtFile, 'r') as f:
        lines = [line.strip() for line in f]

    genes_str = lines[4::6]
    conds_str = lines[5::6]
    genes_int = [list(map(int, line.split())) for line in genes_str]
    conds_int = [list(map(int, line.split())) for line in conds_str]

    genes_set = set()
    conds_set = set()
    for i in range(len(genes_int)):
        genes_set.update(set(genes_int[i]))
        conds_set.update(set(conds_int[i]))
    
    return list(genes_set), list(conds_set)


if __name__ == '__main__':
    datasets = ['BCLL', 'YC', 'PBC', 'RAT']
    algs = ['FA', 'CS', 'CSFA', 'PSO', 'QPSO']
    dataVol = {}
    dataVol['YC'] = 5847*50
    dataVol['RAT'] = 7751*122
    dataVol['BCLL'] = 12185*21
    dataVol['PBC'] = 21225*286
    root = os.getcwd()

    result = {}
    for dataName in datasets:
        result[dataName] = {}
        for alg in algs:
            file = os.path.join(root, dataName, alg+".txt")
            genes, conds = parse(file)
            result[dataName][alg] = len(genes)*len(conds)/dataVol[dataName]
    print(result)
    df = pd.DataFrame(result)
    df.to_csv('coverRate.csv')
