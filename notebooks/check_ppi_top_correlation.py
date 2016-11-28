from scipy.spatial.distance import pdist, squareform
from clustergrammer import Network
from copy import deepcopy
import pandas as pd
import itertools

def main():

  # load_ppi()

  inst_type = 'ptm45_filter_none'

  calc_top_corr(inst_type)

def calc_top_corr(data_type):
  ''' calc top correlations '''


  filename = '../lung_cellline_3_1_16/lung_cl_all_ptm/precalc_processed/' + \
             data_type + '.txt'

  # print('\n')
  # print(data_type)
  # print('-----------------')
  # load file and export dataframe
  net = deepcopy(Network())
  net.load_file(filename)
  net.swap_nan_for_zero()
  tmp_df = net.dat_to_df()
  df = tmp_df['mat']

  print(df.shape)

  dist_mat = pdist(df, metric='correlation')

  print(len(dist_mat))
  print(type(dist_mat))

  dist_series = pd.Series(data=dist_mat)

  dist_series.sort_values(inplace=True)

  print(dist_series[0:10])

  # generate name combinations
  combo_index = list(itertools.combinations(range(865),2))

  print('num of combo names')
  print(len(combo_index))


def load_ppi():
  ''' load PPI network '''
  filename = '../PPI/PPI.sig'

  f = open(filename, 'r')
  lines = f.readlines()
  f.close()

  print(len(lines))

  for i in range(10):

    if i > 0:
      inst_line = lines[i]

      inst_line = inst_line.strip().split()

      source = inst_line[0]

      target = inst_line[5]


      print('\n')
      print(inst_line)
      print(source)
      print(target)


main()