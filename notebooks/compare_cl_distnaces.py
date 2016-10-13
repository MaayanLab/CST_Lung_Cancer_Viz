def main():
  '''
  This script will make cell-line by cell-line distnace vectors and using the
  gene expression data (with and withouth gene-zscoring) and PTM data. I'll
  then check how different PTM data processing methods (normalization/filtering)
  affect the distances between all cell line pairs in
  gene-expression space
  '''

  compare_cl_data_to_ptm(exp_type='exp_none')

  compare_cl_data_to_ptm(exp_type='exp-plex')

def compare_cl_data_to_ptm(exp_type='exp_none', mantel_method='pearson'):

  # compare distance matrices of cell lines based on different processed
  # versions of exp and ptm data
  #########################################################
  multiple_dist_metrics = ['euclidean','cosine']

  for inst_metric in multiple_dist_metrics:
    compare_exp_to_various_ptm_versions(exp_type=exp_type,
       dist_metric=inst_metric, mantel_method=mantel_method)

def compare_plex_to_exp():
  mantel_test('exp-plex', 'exp_none')

def compare_exp_to_various_ptm_versions(exp_type='exp_none',
      dist_metric='euclidean', mantel_method='pearson'):
  import numpy as np
  import pandas as pd
  from copy import deepcopy

  ptm_norms = [
  'none', 'row-zscore', 'col-qn', 'col-zscore',
  'col-qn_row-zscore', 'col-zscore_row-zscore'
  ]

  # only add filter for PTM data
  filter_before = ['filter_'+i for i in ptm_norms]
  filter_after = [i+'_filter' for i in ptm_norms]
  all_proc = ptm_norms + filter_before + filter_after
  all_proc = ['ptm_'+i for i in all_proc]

  cols = ['cor', 'pval', 'zscore']
  rows = deepcopy(all_proc)
  rows = [i.replace('ptm_','').replace('_',' ') for i in rows]

  mat = np.zeros((len(rows), len(cols)))

  # run mantel tests
  #####################################
  for i in range(len(all_proc)):

    ptm_proc = all_proc[i]

    results = mantel_test(exp_type, ptm_proc, dist_metric=dist_metric,
                          mantel_method=mantel_method)

    mat[i, 0] = results[0]
    mat[i, 1] = results[1]
    mat[i, 2] = results[2]

  # save as tsv
  df = pd.DataFrame(data=mat, columns=cols, index=rows)

  # use dash to encode expression data type
  exp_name = exp_type.replace('_none','').replace('_','-')
  filename = '../lung_cellline_3_1_16/lung_cl_all_ptm/compare_cl_dist/'+\
             'cl_'+exp_name+'_vs_ptm_'+dist_metric+'.txt'
  df.to_csv(filename, sep='\t')

def mantel_test(data_1, data_2, perms=10000, tail='upper',
                dist_metric='euclidean', mantel_method='pearson'):

  import Mantel

  print('compare ' + data_1 + ' to ' + data_2)

  # calculate distance matrices of both matrices
  dist_metric_1 = dist_metric
  dist_metric_2 = dist_metric
  if data_1 == 'exp-plex':
    dist_metric_1 = 'jaccard'
  if data_2 == 'exp-plex':
    dist_metric_2 = 'jaccard'

  dist_mat_1 = calc_cl_dist(data_type=data_1, dist_metric=dist_metric_1)
  dist_mat_2 = calc_cl_dist(data_type=data_2, dist_metric=dist_metric_2)

  # pearson or spearman
  results = Mantel.test(dist_mat_1, dist_mat_2, perms=perms, tail='upper', method=mantel_method)

  print(results)
  print('\n')

  return results

def calc_cl_dist(data_type='exp_none', dist_metric='euclidean'):
  '''
  calculate cell line distance based on data_type (e.g. expression) with
  optional filtering and normalization
  '''
  from scipy.spatial.distance import pdist, squareform
  from clustergrammer import Network
  from copy import deepcopy

  net = deepcopy(Network())
  filename = '../lung_cellline_3_1_16/lung_cl_all_ptm/precalc_processed/' + \
             data_type + '.txt'

  # load file and export dataframe
  net.load_file(filename)
  net.swap_nan_for_zero()
  tmp_df = net.dat_to_df()

  df = tmp_df['mat']

  print('data_type: ' + data_type + ' dist_metric: ' + dist_metric + ' shape'+ str(df.shape))

  # transpose to calc distance matrix of columns
  df = df.transpose()

  dist_mat = pdist(df, metric=dist_metric)

  return dist_mat

main()