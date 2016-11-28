def main():

  load_ppi()

  inst_type = 'ptm45'

  calc_top_corr(inst_type)

def calc_top_corr(data_type):
  ''' calc top correlations '''
  from scipy.spatial.distance import pdist, squareform
  from clustergrammer import Network
  from copy import deepcopy


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