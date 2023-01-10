import pandas as pd
import numpy as np

import scipy.stats as stats

# can accept: numpy array, pandas data frame, or list
def merge_p_values(scores, method='Fisher', scores_direction = None, expected_direction = None):
    # validate types of scores variable
    if type(scores) == list:
        scores = np.array(scores)
    if type(scores) != pd.DataFrame or type(scores) != np.ndarray:
        print("scores must be a pandas data frame, numpy array, or list")
        exit()
    if np.any(pd.isnull(scores)):
        print("Scores may not contain missing values")
        exit()
    # check data type float
    try:
        scores = scores.astype('float')
    except:
        print('Scores must be numeric')
        exit()
    # check range of data
    if np.any(scores < 0) or np.any(scores > 1):
        print('All values in scores must be in [0,1]')
        exit()
    # check merging method
    if method not in ["Fisher", "Brown", "Stouffer", "Strube"]:
        print("Only Fisher's, Brown's, Stouffer's and Strube's methods are currently supported")
        exit()
    # if only some p-values are directional, expected_direction must be a pandas dataframe or series
    if type(scores_direction) != pd.DataFrame:
        print("scores_direction must be a pandas Dataframe column labels matching the column names of the p-values dataframe or numpy array")
        exit()

    # convert zeroes to smallest available floats
    scores[scores == 0] = 1e-300

    # if scores is a numpy array
    if type(scores) == np.ndarray:
        if method in ['Brown', 'Strube']:
            print("Brown's or Strube's method cannot be used with a single list of p-values")
            exit()
        
        if method == 'Fisher':
            p_val = 1 - stats.chi2.cdf((fishersMethod(scores, scores_direction, expected_direction), 2*len(scores)))
            return p_val
        
        if method == 'Stouffer':
            p_val = 2 * (1 - stats.norm.cdf(abs(stouffersMethod(scores, scores_direction, expected_direction))))
            return p_val

    # if there is only one column, just return those p-values
    if len(scores.columns) == 1:
        return scores

    # if scores is a matrix with multiple columns, apply the following methods
    if method == 'Fisher':
        p_val = np.array(list(map(lambda x: 1 - stats.chi2.cdf((fishersMethod(scores.loc[x,:], scores_direction.loc[x,:], expected_direction), 2*len(scores))), scores.index.to_numpy())))

        # return dataframe
        return pd.DataFrame(p_val, index = scores.index)

    if method == 'Brown':
        cov_matrix = calculate_covariances(scores.transpose())
        p_val = brownsMethod(scores, scores_direction, expected_direction, cov_matrix = cov_matrix)

        return pd.DataFrame(p_val, index = scores.index)

    if method == 'Stouffer':
        # NOTE subset scores_direction and expected_direction
        p_val = np.array(list(map(lambda x: 2 * (1 - stats.norm.cdf(abs(stouffersMethod(scores.loc[x,:].to_numpy(), scores_direction.loc[x,:], expected_direction)))), scores.index.to_numpy())))

        # return dataframe
        return pd.DataFrame(p_val, index = scores.index)

    if method == 'Strube':
        p_val = brownsMethod(scores, scores_direction, expected_direction)

        return pd.DataFrame(p_val, index = scores.index)




def fishersMethod(p_values, scores_direction = None, expected_direction = None):
    if scores_direction != None and expected_direction != None:
        # apply directionality penalty where applicable
        directionality = expected_direction[direction_mask] * scores_direction[direction_mask] / np.abs(scores_direction[direction_mask])
        p_values_directional = p_values[direction_mask]

        chisq_directional = np.abs(-2 * sum(np.log(p_values_directional) * directionality))

        # calculate for non-directionality p-values
        if np.any(~direction_mask):
            p_values_nondirectional = p_values[~direction_mask]
            chisq_directional = np.concatenate([np.abs(-2 * sum(np.log(p_values_nondirectional))), chisq_directional])

        return np.sum(chisq_directional)

    else:
        return (-2 * np.log(p_values))


def brownsMethod(p_values, scores_direction = None, expected_direction = None, data_matrix=None, cov_matrix=None):
    # what has been provided
    if data_matrix == None and cov_matrix == None:
        print("Either data_matrix or cov_matrix must be supplied")
        exit()
    if data_matrix != None and cov_matrix != None:
        print("Both data_matrix and cov_matrix were supplied. Ignoring data_matrix")
    if cov_matrix == None:
        cov_matrix = calculateCovariances(data_matrix)

    n_datasets = len(cov_matrix.shape[1])
    expected = 2 * n_datasets
    cov_sum = 2 * sum(np.tril(cov_matrix, k=-1))
    var = 4 * n_datasets + cov_sum
    sf = var / (2 * expected)

    df = (2 * (expected ** 2)) / var
    if df > 2 * n_datasets:
        df = 2 * n_datasets
        sf = 1

    try: 
        p_val = np.array(list(map(lambda x: 1 - stats.chi2.cdf((fishersMethod(p_values.loc[x,:], scores_direction.loc[x,:], expected_direction), 2*len(p_values))), p_values.index.to_numpy())))
    except:
        p_val = np.array(list(map(lambda x: 1 - stats.chi2.cdf((fishersMethod(p_values.loc[x,:]), 2*len(p_values))), p_values.index.to_numpy())))

    p_brown = 1 - stats.chi2.cdf(p_val / sf, df)
    return p_brown

