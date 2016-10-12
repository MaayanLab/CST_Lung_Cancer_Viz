def main():
  '''
  This script will make cell-line by cell-line distnace vectors and using the
  gene expression data (with and withouth gene-zscoring) and PTM data. I'll
  then check how different PTM data processing methods (normalization/filtering)
  affect the distances (and/or similarities) between all cell line pairs in
  gene-expression space
  '''

  # compare similarity matrices of cell lines based on different processed
  # versions of exp and ptm data

  # no normalizations
  #######################
  results = mantel_test('exp_none', 'ptm_none')

  print('--- one process \n---------------------')
  results = mantel_test('exp_none', 'ptm_row-zscore')
  results = mantel_test('exp_none', 'ptm_col-qn')
  results = mantel_test('exp_none', 'ptm_filter_none')
  results = mantel_test('exp_none', 'ptm_col-zscore')

  print('--- two processes \n---------------------')
  results = mantel_test('exp_none', 'ptm_col-qn_row-zscore')
  results = mantel_test('exp_none', 'ptm_col-zscore_row-zscore')

  print('--- three processes \n---------------------')
  results = mantel_test('exp_none', 'ptm_filter_col-qn_row-zscore')
  results = mantel_test('exp_none', 'ptm_filter_col-zscore_row-zscore')
  results = mantel_test('exp_none', 'ptm_col-qn_row-zscore_filter')
  results = mantel_test('exp_none', 'ptm_col-zscore_row-zscore_filter')

def mantel_test(data_1, data_2, perms=10000, tail='upper'):

  import Mantel

  print('compare ' + data_1 + ' to ' + data_2)

  # calculate similarity matrices of both matrices
  sim_1 = calc_cl_sim(data_type=data_1)
  sim_2 = calc_cl_sim(data_type=data_2)

  results = Mantel.test(sim_1, sim_2, perms=perms, tail='upper')

  print(results)
  print('\n')

  return results

def calc_cl_sim(data_type='exp_none', dist_metric='euclidean'):
  '''
  calculate cell line similarity based on data_type (e.g. expression) with
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

  print('shape'+ str(df.shape))

  # transpose to calc distance matrix of columns
  df = df.transpose()

  # calculate the similarity of cell line data based on gene expression
  sim = 1 - pdist(df, metric=dist_metric)

  return sim

main()