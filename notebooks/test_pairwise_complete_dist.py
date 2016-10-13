def main():
  print('I will set up a pairwise-compete matrix')

  data_type = 'ptm_none'
  dist_metric = 'euclidean'
  compare_pdist_to_custom_sim_mat(data_type=data_type, dist_metric=dist_metric)

def compare_pdist_to_custom_sim_mat(data_type='ptm_none', dist_metric='euclidean'):
  '''
  calculate cell line similarity based on data_type (e.g. expression) with
  optional filtering and normalization
  '''

  filename = '../lung_cellline_3_1_16/lung_cl_all_ptm/precalc_processed/' + \
             data_type + '.txt'

  sim_pdist = calc_pdist_sim(filename, data_type, dist_metric)

  sim_custom = calc_custom_sim(filename, data_type, dist_metric)

def calc_custom_sim(filename, data_type, dist_metric):
  import numpy as np
  import scipy.spatial.distance as dist_fun

  df = get_df(filename)

  rows = df.index.tolist()
  cols = df.columns.tolist()

  dist_vector = np.zeros(666,)

  # write for-loop to calculate custom distance matrix and compare result
  # to pdist
  for col_1 in cols:

    vect_1 = df[col_1]
    vect_2 = vect_1

    if dist_metric == 'euclidean':
      inst_dist = dist_fun.euclidean(vect_1, vect_2)
    elif dist_metric == 'cosine':
      inst_dist = dist_fun.cosine(vect_1, vect_2)

    print(inst_dist)

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

  return df

main()