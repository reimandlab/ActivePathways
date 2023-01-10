import pandas as pd
import numpy as np

import scipy.stats as stats
from statsmodels.distributions.empirical_distribution import ECDF

# can accept: numpy array, pandas data frame, or list
def merge_p_values(scores, method='Fisher', scores_direction = None, expected_direction = None):
    # validate types of scores variable
    if type(scores) == list:
        scores = np.array(scores)
    elif type(scores) == pd.Series:
        scores = scores.to_numpy()
    if type(scores) != pd.DataFrame or type(scores) != np.ndarray:
        print("scores must be a pandas data frame, series, numpy array, or list")
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
    
    # check scores_direction and expected_direction
    if np.logical_xor(scores_direction == None, expected_direction == None):
        print("Both scores_direction and expected_direction must be provided")
        exit()
    
    if scores_direction != None:
        # if only some p-values are directional, expected_direction must be a pandas dataframe or series
        if type(scores) != pd.Series or type(scores) != pd.DataFrame:
            print("Scores must be a pandas series or dataframe to apply directional penalty")
            exit()

        # check expected direction
        if type(expected_direction) == list:
            expected_direction = np.array(expected_direction)
        if type(expected_direction) != np.array or type(expected_direction) != pd.Series or np.any(np.abs(expected_direction) != 1):
            print("expected_direction must be a numeric type array composed of only -1's or 1's")
            exit()

        # check scores_direction
        if type(scores_direction) != pd.DataFrame:
            print("scores_direction must be a pandas Dataframe column labels matching the column names of the p-values dataframe or numpy array")
            exit()
        if np.any(scores_direction.isnull()):
            print("scores_direction must not contain any missing values. Fill NA's with 1 before processing")
            exit()
        try:
            scores_direction = scores_direction.astype('float')
        except:
            print("scores_direction must be numeric")
            exit()
        if np.any(~np.isin(scores_direction.index, scores.index)) or np.any(~np.isin(scores_direction.columns, scores.columns)):
            print("scores_direction index and columns must match scores index and columns")
            exit()
        if scores_direction.shape[1] < 2:
            print("A minimum of two datasets from the scores matrix should have corresponding directionality data in scores_direction. Ensure column names are identical")
            exit()

        # check both
        if scores_direction.shape[1] != len(expected_direction):
            print("expected_direction should have the same number of entries as columns in scores_direction")
        
    

    # convert zeroes to smallest available floats
    scores[scores == 0] = 1e-300

    # if scores is a numpy array
    if type(scores) == np.ndarray:
        if method in ['Brown', 'Strube']:
            print("Brown's or Strube's method cannot be used with a single list of p-values")
            exit()
        
        if method == 'Fisher':
            p_val = 1 - stats.chi2.cdf((fishers_method(scores, scores_direction, expected_direction), 2*len(scores)))
            return p_val
        
        if method == 'Stouffer':
            p_val = 2 * (1 - stats.norm.cdf(abs(stouffers_method(scores, scores_direction, expected_direction))))
            return p_val

    # if there is only one column, just return those p-values
    if len(scores.columns) == 1:
        return scores

    # if scores is a matrix with multiple columns, apply the following methods
    if method == 'Fisher':
        p_val = np.array(list(map(lambda x: 1 - stats.chi2.cdf((fishers_method(scores.loc[x,:], scores_direction.loc[x,:], expected_direction), 2*len(scores))), scores.index)))

        # return dataframe
        return pd.DataFrame(p_val, index = scores.index)

    if method == 'Brown':
        cov_matrix = calculate_covariances(scores.transpose())
        p_val = browns_method(scores, scores_direction, expected_direction, cov_matrix = cov_matrix)

        return pd.DataFrame(p_val, index = scores.index)

    if method == 'Stouffer':
        # NOTE subset scores_direction and expected_direction
        p_val = np.array(list(map(lambda x: 2 * (1 - stats.norm.cdf(abs(stouffers_method(scores.loc[x,:].to_numpy(), scores_direction.loc[x,:], expected_direction)))), scores.index)))

        # return dataframe
        return pd.DataFrame(p_val, index = scores.index)

    if method == 'Strube':
        p_val = strubes_method(scores, scores_direction, expected_direction)

        return pd.DataFrame(p_val, index = scores.index)


def fishers_method(p_values, scores_direction = None, expected_direction = None):
    if scores_direction != None:
        # apply directionality penalty where applicable
        directionality = expected_direction * scores_direction / np.abs(scores_direction)

        direction_mask = np.isin(p_values.index, scores_direction.index)
        p_values_directional = p_values[direction_mask]

        chisq_directional = np.abs(-2 * np.sum(np.log(p_values_directional) * directionality))

        # calculate for non-directionality p-values
        if np.any(~direction_mask):
            p_values_nondirectional = p_values[~direction_mask]
            chisq_directional = np.concatenate([np.abs(-2 * sum(np.log(p_values_nondirectional))), chisq_directional])

    else:
        chisq_directional = -2 * np.log(p_values)

    return np.sum(chisq_directional)


def browns_method(p_values, scores_direction = None, expected_direction = None, data_matrix=None, cov_matrix=None):
    # what has been provided
    if data_matrix == None and cov_matrix == None:
        print("Either data_matrix or cov_matrix must be supplied")
        exit()
    if data_matrix != None and cov_matrix != None:
        print("Both data_matrix and cov_matrix were supplied. Ignoring data_matrix")
    if cov_matrix == None:
        cov_matrix = calculate_covariances(data_matrix)

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
        p_val = np.array(list(map(lambda x: 1 - stats.chi2.cdf((fishers_method(p_values.loc[x,:], scores_direction.loc[x,:], expected_direction), 2*len(p_values))), p_values.index.to_numpy())))
    except:
        p_val = np.array(list(map(lambda x: 1 - stats.chi2.cdf((fishers_method(p_values.loc[x,:]), 2*len(p_values))), p_values.index.to_numpy())))

    p_brown = 1 - stats.chi2.cdf(p_val / sf, df)

    return p_brown


def stouffers_method(p_values, scores_direction = None, expected_direction = None):
    k = len(p_values)

    if scores_direction != None:
        # apply directionality penalty where applicable
        directionality = expected_direction * scores_direction / np.abs(scores_direction)

        direction_mask = np.isin(p_values.index, scores_direction.index)
        p_values_directional = p_values[direction_mask]

        z_directional = np.abs(sum(stats.norm.ppf(p_values_directional / 2) * directionality))

        # calculate for non-directionality p-values
        if np.any(~direction_mask):
            p_values_nondirectional = p_values[~direction_mask]
            z_directional = np.concatenate([np.abs(np.sum(stats.norm.ppf(p_values_nondirectional / 2))), z_directional])

    else:
        z_directional = stats.norm.ppf(p_values / 2)

    return np.sum(z_directional) / np.sqrt(k)


def strubes_method(p_values, scores_direction = None, expected_direction = None):
    # acquire the unadjusted z-value from stouffer's method
    try:
        stouffers_z_arr = np.array(list(map(lambda x: stouffers_method(p_values.loc[x,:], scores_direction.loc[x,:], expected_direction), p_values.index)))
    except:
        stouffers_z_arr = np.array(list(map(lambda x: stouffers_method(p_values.loc[x,:]), p_values.index)))

    k = p_values.shape[1]

    # correlation matrix
    cor_mtx = np.abs(p_values.corr())
    cor_mtx.loc[cor_mtx.isna()] = 0

    # adjust p-value
    adjusted_z = stouffers_z_arr * np.sqrt(k) / np.sqrt(np.sum(cor_mtx))
    p_strube = 2 * (1 - stats.norm.cdf(abs(adjusted_z)))

    return p_strube


def transform_data(data):
    # If all values in dat are the same (equal to y), return dat. The covariance matrix will be the zero matrix, and brown's method gives the p-value as y. Otherwise (dat - dmv) / dvsd is NaN and ecdf throws an error
    if np.all(data == data[0]):
        return data

    # NOTE no idea what's going on here copied directly from R
    def pop_sd(x):
        x = np.nanvar(x) * (len(x) - 1) / len(x)
        np.sqrt(x)
        return x 
    
    dvm = np.nanmean(data)
    dvsd = pop_sd(data)
    s = (data - dvm) / dvsd
    distr = ECDF(s)

    return np.array(list(map(lambda x: -2 * np.log(distr(x)), s)))

def calculate_covariances(data_matrix):
    transformed_matrix = pd.DataFrame(list(map(lambda x: transform_data(data_matrix.loc[x,:].to_numpy()), data_matrix.index))).transpose()
    
    return transformed_matrix.cov()


