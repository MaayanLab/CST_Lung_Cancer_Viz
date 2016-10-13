def main():
  print('I will set up a pairwise-compete matrix')

  inst_data_type = 'ptm_none'
  custom_dist_matrix(inst_data_type)

def custom_dist_matrix(data_type):
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

main()