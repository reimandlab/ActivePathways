import numpy as np
import scipy.stats as stats
from sys import exit

def hypergeomtric(counts):
    if np.any(counts < 0):
        print("counts contains negative values. Something went very wrong.")
        exit()
    
    # NOTE check this
    m = counts[0,0] + counts[1,0]
    n = counts[0,1] + counts[1,1]
    k = counts[0,0] + counts[0,1]
    x = counts[0,0]

    return 1 - stats.hypergeom(m+n, k, m).cdf(x)

def orderedHypergeometric(gene_arr, background, annotations):
    # Only test subsets of genelist that end with a gene in annotations since these are the only tests for which the p-value can decrease
    genes_indices = np.where(gene_arr, annotations)[0]
    if len(genes_indices) == 0:
        return 1, 1

    gl = gene_arr[1:]
    cl = background[~(np.isin(background, gene_arr))]
    genelist0 = len(gl) - 1
    complement1 = np.sum(np.isin(cl, annotations))
    complement0 = len(cl) - complement1

    counts = np.array([[1, complement1], [genelist0, complement0]])
    scores = [hypergeomtric(counts)]

    if len(genes_indices) == 1:
        return scores[0], genes_indices


    for i in range(1,len(genes_indices)):
        diff = genes_indices[i] - genes_indices[i - 1]

        counts[0,0] = i
        counts[1,0] = counts[1,0] + diff - 1
        counts[0,1] = counts[0,1] - 1
        counts[1,1] = counts[1,1] - diff + 1

        scores.append(hypergeomtric(counts))

    
    min_score = np.min(scores)

    ind = genes_indices[np.max(np.where(scores == min_score))]
    p_val = min_score

    return p_val, ind


