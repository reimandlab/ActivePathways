import pandas as pd
import numpy as np
from sys import exit
import scipy.stats as stats
import statsmodels.stats.multitest as multi
import matplotlib as mpl
import matplotlib.pyplot 
from gmt import *
from statistical_tests import *
from cytoscape import *
from merge_p import *

import time

def ActivePathways(scores, gmt, background = None,
                   geneset_filter = [5,1000], cutoff = 0.1, significant = 0.05,
                   merge_method = ["Brown","Fisher","Stouffer","Strube"],
                   correction_method = ["holm", "hommel","bonferroni", "fdr_bh", "fdr_by", "none"],
                   cytoscape_file_tag = None, color_palette = None, custom_colors = None, 
                   color_integrated_only = "#FFFFF0", scores_direction = None, expected_direction = None):
    gmt_old = gmt
    
    ##### Validation #####
    # scores
    if type(scores) != pd.DataFrame and type(scores) != pd.Series:
        print("scores must be a pandas data frame or series")
        return
    if np.any(pd.isnull(scores)):
        print("scores may not contain missing values")
        return
    try:
        scores = scores.astype("float")
    except:
        print("scores must be numeric")
        return
    if np.any(scores < 0) or np.any(scores > 1):
        print("All values in scores must be in [0,1]")
        return
    if np.any(scores.index.duplicated()):
        print("scores contains duplicated genes - row names must be unique")
        return
    
    # scores_direction and expected_direction
    if np.logical_xor(scores_direction == None, expected_direction == None):
        print("Both scores_direction and expected_direction must be provided")
        return
    if scores_direction != None and type(scores) == pd.DataFrame:
        if type(expected_direction) == list:
            expected_direction = np.array(expected_direction)
        if type(expected_direction) != np.ndarray and type(expected_direction) != pd.Series:
            print("expected_direction must be a numeric type array")
            return
        if not np.all(np.isin(expected_direction,[1,-1,0])):
            print("expected_direction must contain the values: 1, -1 or 0")
            return
        if type(scores_direction) != pd.DataFrame:
            print("scores_direction must be a pandas data frame with column labels matching the column names of the scores data frame")
            return
        if np.any(scores_direction.isnull()):
            print("scores_direction must not contain any missing values. Fill NA's with 1 before processing")
            return
        try:
            scores_direction = scores_direction.astype("float")
        except:
            print("scores_direction must be numeric")
            return
        if np.any(not np.isin(scores_direction.index, scores.index)) or np.any(not np.isin(scores_direction.columns, scores.columns)):
            print("scores_direction index and columns must match scores index and columns")
            return
        if scores_direction.shape[1] != len(expected_direction):
            print("expected_direction should have the same number of entries as columns in scores_direction")
            return    
        if np.any(np.isin(expected_direction,0)) and not np.all(scores_direction.loc[:,np.isin(expected_direction,0)] == 0):
            print("scores_direction entries must be set to 0's for columns that do not contain directional information")
            return

    # cutoff and significant
    if not isinstance(cutoff,(int,float)):
        print("cutoff must be numeric")
        return
    if not isinstance(significant,(int,float)):
        print("significant must be numeric")
        return
    if cutoff < 0 or cutoff > 1:
        print("cutoff must be a value in [0,1]")
        return
    if significant < 0 or significant > 1:
        print("significant must be a value in [0,1]")
        return
    
    # gmt
    background_provided = True
    if background == None:
        background = make_background(read_gmt(gmt))
        background_provided = False
    if not is_gmt(gmt):
        gmt = read_gmt(gmt)
    if type(background) == list:
        background = np.array(background)
    if type(background) != np.ndarray and type(background) != pd.Series:
        print("background must be an array")
        return
    if not np.all(list(map(lambda x: isinstance(x,str),background))):
        print("background must contain string values")
        return
    
    # geneset_filter
    if geneset_filter != None:
        if type(geneset_filter) == list:
            geneset_filter = np.array(geneset_filter)
        if type(geneset_filter) != np.ndarray and type(geneset_filter) != pd.Series:
            print("geneset_filter must be an array")
            return
        if not np.all(list(map(lambda x: isinstance(x,(int, np.integer)), geneset_filter))):
            print("geneset_filter must be numeric")
            return
        if len(geneset_filter) != 2:
            print("geneset_filter must be length 2")
            return
        if np.any([x < 0 for x in geneset_filter]):
            print("geneset_filter limits must be positive")
            return
              
    # custom_colors
    if custom_colors != None and type(scores) == pd.DataFrame:
        if type(custom_colors) == list:
            custom_colors = np.array(custom_colors)
        if type(custom_colors) != np.ndarray and type(custom_colors) != pd.Series:
            print("custom_colors must be an array")
            return
        if not np.all(list(map(lambda x: isinstance(x,str),custom_colors))):
            print("custom_colors must contain string values")
            return
        if len(scores.columns) != len(custom_colors):
            print("incorrect number of colors is provided")
            return
        if color_palette != None:
            print("Both custom_colors and color_palette are provided. Specify only one of these parameters for node coloring.")
            return
        if type(custom_colors) == pd.Series:
            if not np.all(np.isin(custom_colors.index,scores.columns)):
                print("The custom_colors index names should match the scores column names")
                return
    
    # color_palette
    if color_palette != None:
        if color_palette not in mpl.pyplot.colormaps():
            print("The color palette must be one of the colormaps in the matplotlib package")
            return
              
    # color_integrated_only
    if not isinstance(color_integrated_only,str):
        print("color_integrated_only must be provided as a string")
        return
    if len([color_integrated_only]) != 1:
        print("only a single color must be specified")
        return
    
    # contribution
    contribution = True
    if type(scores) == pd.Series or (type(scores) == pd.DataFrame and len(scores.columns) == 1):
        contribution = False
        print("scores contains only one column. Column contributions will not be calculated")
              
    
    ##### Filtering and sorting #####
    # Remove any genes not found in the background
    orig_length = len(scores.index)
    scores = scores[scores.index.isin(background)]
    if scores_direction != None:
        scores_direction = scores_direction[scores_direction.index.isin(background)]
    if len(scores.index) == 0:
        print("scores does not contain any genes in the background")
        return
    if len(scores.index) < orig_length:
        print(f"{orig_length - len(scores.index)} rows were removed from scores because they are not found in the background")
    


    # how many gene sets do we start with
    origin_num_gene_sets = len(gmt.genes_dict.keys())


    # remove gene sets beyond the range of max and min number of genes
    key_lengths = np.array([len(x) for x in gmt.genes_dict.values()])
    mask = np.logical_and(key_lengths >= geneset_filter[0], key_lengths <= geneset_filter[1])

    # a dictionary of keys that we keep
    good_keys = dict(zip(np.array(list(gmt.genes_dict.keys()))[mask], np.zeros(np.sum(mask))))

    gmt.genes_dict = {key:gmt.genes_dict[key] for key in good_keys}
    gmt.name_dict = {key:gmt.name_dict[key] for key in good_keys}

    # if the background is custom, make sure we exclude genes not in background
    if background_provided:
        for key, gene_list in list(gmt.genes_dict.items()):
            common_genes = np.intersect1d(gene_list, background)

            # remove if there are no genes, or if the number of genes is outside parameters
            if len(common_genes) == 0:
                gmt.genes_dict.pop(key)
                gmt.name_dict.pop(key)

            else:
                gmt.genes_dict[key] = common_genes


    # count number of pathways that passed
    if len(gmt.genes_dict.keys()) == 0:
        print("No pathways in gmt made the geneset_filter")
        return
    if origin_num_gene_sets != len(gmt.genes_dict.keys()):
        print(f"{origin_num_gene_sets - len(gmt.genes_dict.keys())} terms were removed from gmt because they did not make the geneset_filter")
    
    # merge p-values to get a single score for each gene and remove any genes that don't make the cutoff
    merged_scores = merge_p_values(scores, merge_method, scores_direction, expected_direction)
    print(merged_scores)
    merged_scores = merged_scores[merged_scores[0] <= cutoff]

    if merged_scores.index.size == 0:
        print("No genes made the cutoff")
        return

    # sort genes by p-value (convert series to np array)
    merged_scores = merged_scores.sort_values(0).index.to_numpy(dtype="<U64")

              
    ##### enrichmentAnalysis and column contribution #####
    # create a dataframe 
    res = enrichmentAnalysis(merged_scores, gmt, background)
    res.loc[:,'adjusted_p_val'] = adjust_p_value(res.adjusted_p_val, method = correction_method)
    
    significant_indices = res[res.adjusted_p_val <= significant].index
    if len(significant_indices) == 0:
        print("No significant terms were found")
        return None
    
    if contribution:
        sig_cols = columnSignificance(scores, gmt, background, cutoff, significant, correction_method, res.adjusted_p_val)
        res = pd.concat([res, sig_cols.iloc[:,-1:]], axis=1)
    else:
        sig_cols = None
    
    # if significant results were found and cytoscape file tag exists
    # proceed with writing files in the working directory
    if len(significant_indices) > 0 and cytoscape_file_tag != None:
        significant_ids = [list(gmt.genes_dict)[i] for i in significant_indices]
        genes_dict_significant = {key : gmt.genes_dict[key] for key in significant_ids}
        name_dict_significant = {key : gmt.name_dict[key] for key in significant_ids}
        gmt_significant = GMT(name_dict_significant, genes_dict_significant)
              
        prepareCytoscape(res.loc[significant_indices, ["term_id","term_name","adjusted_p_val"]],
                         gmt_significant,
                         cytoscape_file_tag,
                         sig_cols.iloc[significant_indices, :], color_palette, custom_colors, color_integrated_only)
              
    return res.iloc[significant_indices, :]

              
# Perform pathway enrichment analysis on an ordered list of genes
def enrichmentAnalysis(genelist, gmt, background):
    df = pd.DataFrame(columns = ['term_id', 'term_name', 'adjusted_p_val', 'term_size', 'overlap']).astype('object')
    # df = pd.DataFrame().astype("object")

    for i,(key, gmt_genes) in enumerate(gmt.genes_dict.items()):
        tmp = orderedHypergeometric(genelist, background, gmt_genes)
        overlap = genelist[0:tmp[1]]
        overlap = np.array(np.intersect1d(overlap,gmt_genes),dtype="<U64")
        df.loc[i,'term_id'] = key
        df.loc[i,'term_name'] = gmt.name_dict[key]
        df.loc[i,'adjusted_p_val'] = tmp[0]
        df.loc[i,'term_size'] = len(gmt_genes)
        df.loc[i,'overlap'] = overlap
    return df
    
# Correct p-values by various multitest correction methods
def adjust_p_value(pvals, method):
    if method == "none":
        return np.array(pvals)
    else:
        return multi.multipletests(pvals,method = method)[1]
              
# Determine which terms are found to be significant using each column individually
def columnSignificance(scores, gmt, background, cutoff, significant, correction_method, pvals):
    df = pd.DataFrame({'term_id':gmt.genes_dict.keys(), 'evidence' : None})
    for col in scores.columns:
        col_scores = scores.loc[:,col]
        col_scores = col_scores[col_scores <= cutoff]
        col_scores = col_scores.sort_values().index.to_numpy(dtype="<U64")
        
        res = enrichmentAnalysis(col_scores, gmt, background)
        res.loc[:,'adjusted_p_val'] = adjust_p_value(res.adjusted_p_val, method = correction_method)
        res.loc[:,'overlap'][res.adjusted_p_val > significant] = None
        df.loc[:,col] = res.overlap
    
    ev_names = df.iloc[:,2:].columns
    def set_evidence(x):
        ev = ev_names[~df.iloc[x,2:].isnull()]
        if len(ev) == 0:
            if pvals[x] <= significant:
                ev = 'combined'
            else:
                ev = 'none'
        return ev
    evidence = list(map(lambda x: set_evidence(x), range(len(df.index))))
    df.loc[:,'evidence'] = evidence
    df.columns = ['term_id','evidence'] + list(map(lambda x: 'Genes_' + x, df.columns[2:]))  
    
    return df 
    
