from clustergrammer import Network

import pandas as pd
import numpy as np
from copy import deepcopy

def main():

  inst_data_type = 'ptm45_none'

  inst_df = load_data_as_df(inst_data_type)

  print(inst_df.shape)

  iqn_df(inst_df)

def load_data_as_df(data_type):


  filename = '../lung_cellline_3_1_16/lung_cl_all_ptm/precalc_processed/' + \
             data_type + '.txt'

  # load file and export dataframe
  net = deepcopy(Network())
  net.load_file(filename)
  # net.swap_nan_for_zero()
  tmp_df = net.dat_to_df()
  df = tmp_df['mat']

  return df

def iqn_df(df, axis='row', keep_orig=False):
  '''
  Calculate quantile normalization and remove the influence of missing data
  '''

  # get coumn names
  cols = df.columns.tolist()

  # initialize dataframe for common data
  com_dist = pd.DataFrame()

  df_meas = deepcopy(df)
  df_meas[np.isnan(df_meas)==False] = 1
  df_meas[np.isnan(df_meas)] = 0

  # get the length of the common distribution for common data
  com_dist_len = min(df_meas.sum())

  print(com_dist_len)


main()