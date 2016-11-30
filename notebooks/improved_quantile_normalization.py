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
  all_col = df.columns.tolist()

  # initialize dataframe for common data
  com_dist = pd.DataFrame()

  df_meas = deepcopy(df)
  df_meas[np.isnan(df_meas)==False] = 1
  df_meas[np.isnan(df_meas)] = 0

  # get the length of the common distribution for common data
  com_dist_len = int(min(df_meas.sum()))

  print('\nthe minimum number of measured PTMs in any cell line')
  print(com_dist_len)

  # gather the mesured cell line data into series
  meas_col = {}
  for inst_col in all_col:

    # get non-nan values and sort in place
    inst_series = df[inst_col]
    inst_series = inst_series.dropna()

    # save to dictionary
    meas_col[inst_col] = inst_series

    print(inst_col + ' column has ' + str(len(inst_series)) + ' measured values')

  # map measurements from each col onto com_dist_len number of measuremesnts so
  # that they can be averaged into a common distribution series with duplicate
  # indexes
  map_series_tmp = {}
  # series with unique indexes
  map_series = {}

  for inst_col in all_col:

    # get series
    inst_series = meas_col[inst_col]

    # loop through series and add to temporary series with duplicate indexes
    # in general
    inst_series_len = len(inst_series)

    # initialize series
    map_series[inst_col] = pd.Series()

    # make two lists to cnostruct series, list_values and list_indexes
    list_values = []
    list_indexes = []

    # gather sorted indexes and values

    # i is the original sorted index of the data
    for i in range(inst_series_len):

      # save values
      inst_value = inst_series[i]
      list_values.append(inst_value)

      # map i is the mapped sorted index - mapped so taht the current index
      # fits into the common index series
      map_i = int( round( (i/float(inst_series_len))*com_dist_len ) )

      # save indexes as strings
      # the 'names' are their mapped index
      list_indexes.append(str(map_i))

    # save values and indexes to series
    map_series_tmp[inst_col] = pd.Series(list_values, index=list_indexes)

    # average repeat values
    for i in range(com_dist_len):

      # define the current map index
      map_i = str(i)

      # grab the value(s) with this mapped index map_i
      inst_value = map_series_tmp[inst_col][map_i]

      # average over repeated values if necessary
      if type(inst_value) is pd.core.series.Series:
        inst_value = inst_value.mean()

      # add this value to the series
      map_series[inst_col][map_i] = inst_value

main()