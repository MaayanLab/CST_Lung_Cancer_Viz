def main():
  '''
  calculate cell line distance based on data_type (e.g. expression) with
  optional filtering and normalization
  '''
  from scipy.spatial.distance import pdist, squareform
  import pandas as pd
  from clustergrammer import Network
  from copy import deepcopy

  # data_type = 'ptm45_col-qn_row-zscore'
  data_type = 'ptm45_none'

  filename = '../lung_cellline_3_1_16/lung_cl_all_ptm/precalc_processed/' + \
             data_type + '.txt'

  # load file and export dataframe
  net = deepcopy(Network())
  net.load_file(filename)
  net.swap_nan_for_zero()
  tmp_df = net.dat_to_df()
  df = tmp_df['mat']

  # get cell line names
  cols = df.columns.tolist()

  print(df.shape)

  # transpose to calculte cell line distance matrix

  df = df.transpose()

  dist_mat = pdist(df, metric='correlation')

  dist_mat = squareform(dist_mat)

  exp_df = pd.DataFrame(data=dist_mat, index=cols, columns=cols)

main()