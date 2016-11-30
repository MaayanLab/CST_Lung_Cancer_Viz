from scipy.spatial.distance import pdist, squareform
from clustergrammer import Network
from copy import deepcopy
import pandas as pd
import itertools

def main():

  ppi_combos = load_ppi()

  inst_type = 'ptm45_filter_none'

  top_sorted = calc_top_corr(inst_type)

  calc_num_matches(ppi_combos, top_sorted)


def calc_num_matches(ppi_combos, top_sorted):

  print('ppi_combos')
  print(len(ppi_combos))

  print('\n')
  print(set(ppi_combos[0]))

  print('top_sorted')
  print(len(top_sorted))


  for ptm_pair in top_sorted:

    ptm_1 = ptm_pair[0].split('_')[0]
    ptm_2 = ptm_pair[0].split('_')[1]

    gene_pair = set([ptm_1, ptm_2])
    # print( 'checking gene pair ' + str(gene_pair))

    for inst_interaction in ppi_combos:

      inst_interaction = set(inst_interaction)

      if gene_pair == inst_interaction:
        print('found')
        print(gene_pair)
        print('\n\n')



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


  ptm_names = df.index.tolist()

  dist_mat = pdist(df, metric='correlation')

  # generate name combinations
  combo_index = list(itertools.combinations(range(865),2))
  combo_names = []
  for inst_index in combo_index:
    inst_source = ptm_names[inst_index[0]]
    inst_target = ptm_names[inst_index[1]]
    inst_tuple = (inst_source, inst_target)
    combo_names.append(inst_tuple)

  print('num of combo names')
  print(len(combo_index))

  dist_series = pd.Series(data=dist_mat, index=combo_names)

  # sort values
  dist_series.sort_values(inplace=True)

  sorted_names = dist_series.index.tolist()
  sorted_corrs = dist_series.data

  # print(dist_series[0:10])
  # print(sorted_names[0:10])

  # top_sorted = sorted_names[0:1000]
  top_sorted = sorted_names

  return top_sorted

def load_ppi():
  ''' load PPI network '''
  filename = '../PPI/PPI.sig'

  f = open(filename, 'r')
  lines = f.readlines()
  f.close()

  print(len(lines))

  ppi_combos = []

  for i in range(len(lines)):

    if i > 0:
      inst_line = lines[i]

      inst_line = inst_line.strip().split()

      source = inst_line[0]
      target = inst_line[5]

      inst_tuple = ( source, target )

      # print('\n')
      # print(inst_line)
      # print(source)
      # print(target)

      ppi_combos.append(inst_tuple)

  return ppi_combos


main()