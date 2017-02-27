from clustergrammer import Network
import pandas as pd
import numpy as np
from copy import deepcopy

def main():

  inst_data_type = 'ptm45_none'

  inst_df = load_data_as_df(inst_data_type)

  # perform 'improved' quantile normalization on inst_df
  iqn_df(inst_df)

def iqn_df(df, axis='row', keep_orig=False):

  # com_dist = calc_common_dist(df)

  com_dist_filename = 'intermediate_data/ptm45_common_distribution.txt'

  # read to dataframe, no header and row index
  tmp_com_dist = pd.read_table(com_dist_filename, header=None, index_col=0)

  # save to pandas series
  inst_data = tmp_com_dist.values.flatten()
  inst_index = tmp_com_dist.index.tolist()
  com_dist = pd.Series(data=inst_data, index=inst_index)

  qn_df = remap_data_using_common_dist(df, com_dist)

  filename = '../lung_cellline_3_1_16/lung_cl_all_ptm/precalc_processed/' + \
             'ptm45_col-iqn' + '.txt'

  qn_df.to_csv(filename, sep='\t', na_rep='nan')


def remap_data_using_common_dist(df, com_dist):
  '''
  Use the common distribution, loop through all measurements in a cell line and
  swap the values from the common distribution. Start by sorting the ptm
  measurements from a cell line, calculated the adjusted ranks for the ptm, then
  swap in the adjusted value
  '''

  # get the column names
  all_col = df.columns.tolist()

  # initialize qn dataframe
  qn_df = deepcopy(df)

  # the length of the common distribution

  com_dist_len = len(com_dist)

  # ptm_adj: hold adjusted ptm values
  ptm_adj = {}

  for inst_col in all_col:

    # get the measured data from a cell line in a series
    inst_series = df[inst_col].dropna()
    inst_series.sort_values(inplace=True)

    # list sorted ptms
    all_ptms = inst_series.index.tolist()
    inst_series_len = len(all_ptms)

    # initialize the cell line dictionary
    ptm_adj[inst_col] = {}

    for i in range(len(all_ptms)):

      # ptm name
      inst_ptm = all_ptms[i]

      # transform i to adj_rank
      adj_rank = int( round( (i/float(inst_series_len))*com_dist_len ) )

      # get the adjusted value using the adjusted rank
      adj_value = com_dist[adj_rank]

      # save the adjusted value in the dictionary
      ptm_adj[inst_col][inst_ptm] = adj_value

  print('transfering from dictionary to dataframe')

  # add the adjusted values back to the adjusted dataframe
  for inst_col in all_col:

    all_ptms = ptm_adj[inst_col].keys()

    print('-> '+ inst_col)

    # loop through the ptms in each cell line
    for inst_ptm in all_ptms:

      inst_adj_value = ptm_adj[inst_col][inst_ptm]

      # save to qn_df: quantile normalized dataframe
      qn_df[inst_col][inst_ptm] = inst_adj_value

  # print(qn_df.shape)
  # qn_df.to_csv('intermediate_data/qn_df.txt', sep='\t', na_rep='nan')

  return qn_df

def calc_common_dist(df):
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
    inst_series.sort_values(inplace=True)

    # print('there are ' + str(len(inst_series)) + ' values ')
    # print(inst_series)
    # print('\n\n')

    # save to dictionary
    meas_col[inst_col] = inst_series

    # print(inst_col + ' column has ' + str(len(inst_series)) + ' measured values')

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

    # assign mapped indexes to sorted values: map_series_tmp
    ############################################################
    # gather sorted indexes and values
    print('gather sorted indexes and values')
    # i is the original sorted index of the data
    for i in range(inst_series_len):

      # save values
      inst_value = inst_series[i]
      list_values.append(inst_value)

      # map i is the mapped sorted index - mapped so that the current index
      # fits into the common index series
      map_i = int( round( (i/float(inst_series_len))*com_dist_len ) )

      # save indexes as strings
      # the 'names' are their mapped index
      list_indexes.append(str(map_i))

    # save values and indexes to series
    map_series_tmp[inst_col] = pd.Series(list_values, index=list_indexes)

    # make real mapped series: map_series
    ############################################################
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

  # add map_series to dataframe com_dist
  for inst_col in all_col:
    # transfer from dictionary to dataframe
    com_dist[inst_col] = map_series[inst_col]

  # average the columns, axis = 1
  com_dist      = com_dist.mean(axis=1)
  com_dist_dict = com_dist.to_dict()

  com_dist.to_csv('intermediate_data/ptm45_common_distribution.txt', sep='\t', na_rep='nan')

  # return the common distribution
  return com_dist

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

main()