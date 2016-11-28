import numpy as np
import pandas as pd
from scipy.stats import pearsonr

def main():
  data_type = ['ptm45_none', 'ptm45_col-qn', 'ptm45_col-qn_row-zscore','ptm45_filter_none', 'ptm45_filter_col-qn', 'ptm45_filter_col-qn_row-zscore']

  for inst_type in data_type:
    compare_duplicates_to_other(inst_type)

def compare_duplicates_to_other(data_type):
  '''
  calculate cell line distance based on data_type (e.g. expression) with
  optional filtering and normalization
  '''
  from scipy.spatial.distance import pdist, squareform
  from clustergrammer import Network
  from copy import deepcopy


  filename = '../lung_cellline_3_1_16/lung_cl_all_ptm/precalc_processed/' + \
             data_type + '.txt'

  print('\n')
  print(data_type)
  print('-----------------')
  # load file and export dataframe
  net = deepcopy(Network())
  net.load_file(filename)
  net.swap_nan_for_zero()
  tmp_df = net.dat_to_df()
  df = tmp_df['mat']

  print('-----------------')
  print(df.shape)

  # get cell line names
  cols = df.columns.tolist()

  # transpose to calculte cell line distance matrix

  calc_corr(df)

def calc_corr(df):

  cols = df.columns.tolist()

  dist_mat = np.zeros([len(cols), len(cols)])

  rep_corr = []
  other_corr = []

  rep_pval = []
  other_pval = []

  for x_index in range(len(cols)):
    for y_index in range(len(cols)):

      if x_index <= y_index:
        col_1 = cols[x_index]
        col_2 = cols[y_index]

        vect_1 = df[col_1]
        vect_2 = df[col_2]

        corr_info = pearsonr(vect_1, vect_2)

        inst_corr = corr_info[0]
        inst_pval = corr_info[1]

        dist_mat[x_index, y_index] = inst_corr

        if cols[x_index] != cols[y_index]:
          if cols[x_index].split('_')[0] == cols[y_index].split('_')[0]:
            rep_corr.append(inst_corr)
            rep_pval.append(inst_pval)
          else:
            other_corr.append(inst_corr)
            other_pval.append(inst_pval)

  exp_df = pd.DataFrame(data=dist_mat, index=cols, columns=cols)

  exp_df.to_csv('tmp_cl_dist_mat.txt', sep='\t')

  print('compare correlations')
  print(np.mean(other_corr))
  print(np.mean(other_pval))
  print(len(other_pval))
  print('\n')
  print(np.mean(rep_corr))
  print(np.mean(rep_pval))
  print(len(rep_pval))

def calc_pdist(df):
  df = df.transpose()

  dist_mat = pdist(df, metric='correlation')

  print(len(dist_mat))

  print(max(dist_mat))
  print(min(dist_mat))

  dist_mat = squareform(dist_mat)

  # dist_mat = 1 - dist_mat

  exp_df = pd.DataFrame(data=dist_mat, index=cols, columns=cols)

  exp_df.to_csv('tmp_cl_dist_mat.txt', sep='\t')

main()