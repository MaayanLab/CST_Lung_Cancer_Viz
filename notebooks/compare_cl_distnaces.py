def main():
  '''
  This script will make cell-line by cell-line distnace vectors and using the
  gene expression data (with and withouth gene-zscoring) and PTM data. I'll
  then check how different PTM data processing methods (normalization/filtering)
  affect the distances (and/or similarities) between all cell line pairs in
  gene-expression space
  '''

  # calculate similarities of cell lines
  ##########################################
  # # calculte similarity vector based on expression data
  sim_exp = calc_cl_sim(data_type='exp_none')

  # calculate ptm sim
  sim_ptm = calc_cl_sim(data_type='ptm_none')
  sim_ptm_norm = calc_cl_sim(data_type='ptm_col-qn_row-zscore')
  sim_ptm_filt = calc_cl_sim(data_type='ptm_filter_none')

  print('here')
  import Mantel

  results = Mantel.test(sim_exp, sim_ptm, perms=10000, tail='upper')
  print(results)

  results = Mantel.test(sim_exp, sim_ptm_norm, perms=10000, tail='upper')
  print(results)

  results = Mantel.test(sim_exp, sim_ptm_filt, perms=10000, tail='upper')
  print(results)

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

  print(df.shape)

  # transpose to calc distance matrix of columns
  df = df.transpose()

  # calculate the similarity of cell line data based on gene expression
  sim = 1 - pdist(df, metric=dist_metric)

  return sim

main()