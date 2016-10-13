def main():
  print('I will set up a pairwise-compete matrix')

  inst_data_type = 'ptm_none'
  compare_pdist_to_custom_sim_mat(inst_data_type)

def compare_pdist_to_custom_sim_mat(data_type='ptm_none', dist_metric='euclidean'):
  '''
  calculate cell line similarity based on data_type (e.g. expression) with
  optional filtering and normalization
  '''

  filename = '../lung_cellline_3_1_16/lung_cl_all_ptm/precalc_processed/' + \
             data_type + '.txt'

  sim_pdist = calc_pdist_sim(filename, data_type, dist_metric)

  sim_cust = calc_custom_sim(filename, data_type, dist_metric)


def calc_custom_sim(filename, data_type, dist_metric):
  df = get_df(filename)

  # return sim_

def calc_pdist_sim(filename, data_type, dist_metric):
  from scipy.spatial.distance import pdist, squareform

  df = get_df(filename)

  # transpose to calc distance matrix of columns
  df = df.transpose()

  sim_pdist = 1 - pdist(df, metric=dist_metric)

  return sim_pdist

def get_df(filename):
  from copy import deepcopy
  from clustergrammer import Network
  net = deepcopy(Network())

  # load file and export dataframe
  net.load_file(filename)
  net.swap_nan_for_zero()
  tmp_df = net.dat_to_df()

  df = tmp_df['mat']

  print(df.shape)
  return df

main()