def main():
  '''
  This script will make cell-line by cell-line distnace vectors and using the
  gene expression data (with and withouth gene-zscoring) and PTM data. I'll
  then check how different PTM data processing methods (normalization/filtering)
  affect the distances (and/or similarities) between all cell line pairs in
  gene-expression space
  '''

  # # calculte similarity vector based on expression data
  # sim_exp = calc_cl_sim(data_type='exp')

  # # sim_exp_filt = calc_cl_sim(data_type='exp', var_filter=1000)
  # sim_exp_filt = calc_cl_sim(data_type='exp', var_filter=1000)

  # calculate similarity vector based on ptm data
  save_gene_exp_compatible_ptm_data()

  # # compare similarity vectors based on expression and ptm data
  # sim_data = compare_sim_vectors(sim_exp, sim_exp_filt)
  # print(sim_data)

def save_gene_exp_compatible_ptm_data():
  '''
  This will make a PTM matrix that is compatible with (e.g. easily comparable)
  with the gene expression matrix: 1) average duplicate cell line measurements,
  2) only include the 37 cell lines that are found in CCLE, 3) put cell lines
  into the same order as the gene expression data.
  '''

  combine_and_save_ptm()

def combine_and_save_ptm():
  # combine all ptm data into single dataframe
  #################################################################
  # use simple col names

  from clustergrammer import Network
  import pandas as pd
  from copy import deepcopy

  ptm_data = {
  'phos': '../lung_cellline_3_1_16/lung_cellline_phospho/lung_cellline_TMT_phospho_combined_ratios.tsv',
  'act': '../lung_cellline_3_1_16/lung_cellline_Ack/lung_cellline_TMT_Ack_combined_ratios.tsv'
  }

  df_all = pd.DataFrame()

  for inst_type in ptm_data:

    net = deepcopy(Network())
    filename = ptm_data[inst_type]

    net.load_file(filename)

    tmp_df = net.dat_to_df()
    inst_df = tmp_df['mat']

    col_tuples = inst_df.columns.tolist()

    col_names = []
    for inst_tuple in col_tuples:
      col_names.append(inst_tuple[0])

    inst_df.columns = col_names

    print(inst_df.shape)

    df_all = pd.concat([df_all, inst_df], axis=0)

  filename_all_ptm = '../lung_cellline_3_1_16/lung_cellline_TMT_all_ptm_ratios.tsv'

  print(df_all.shape)

  df_all.to_csv(filename_all_ptm, sep='\t')

  all_rows = df_all.index.tolist()
  print(len(all_rows))
  all_rows = list(set(all_rows))
  print(len(all_rows))




  # sort names in place
  # col_names.sort()

  # check that the cell lines are in the same order in both exp and PTM data


def calc_cl_sim(data_type='exp', sum_filter=None, var_filter=None,
                    row_zscore=False, col_qn=False, col_zscore=False):
  '''
  calculate cell line similarity based on data_type (e.g. expression) with
  optional filtering and normalization
  '''

  from scipy.spatial.distance import pdist, squareform
  from clustergrammer import Network

  net = Network()

  all_data = {
    'exp':'../CCLE_gene_expression/CCLE_NSCLC_all_genes.txt'
  }

  filename = all_data[data_type]

  net.load_file(filename)

  # col normalization
  ######################

  # run col qn before anything
  if col_qn != False:
    print('column qn')
    net.normalize(axis='col', norm_type='qn')

  # run col zscore before anything
  if col_zscore != False:
    print('zscore the cols')
    net.normalize(axis='col', norm_type='zscore')

  # row normalization
  ######################

  # run row zscore after col qn
  if row_zscore != False:
    print('zscore the rows')
    net.normalize(axis='row', norm_type='zscore')

  # row filtering
  ######################

  # filter rows after normalization
  if sum_filter != None:
    print('filter top ' + str(sum_filter) + ' rows based on sum')
    net.filter_N_top('row', sum_filter, rank_type='sum')

  if var_filter != None:
    print('filter top ' + str(var_filter) + ' rows based on variance')
    net.filter_N_top('row', var_filter, rank_type='var')



  tmp_df = net.dat_to_df()

  df = tmp_df['mat']

  col_names = df.columns.tolist()

  print(df.shape)

  # transpose to calc distance matrix of columns
  df = df.transpose()

  # calculate the similarity of cell line data based on gene expression
  sim = 1 - pdist(df, metric='cosine')

  return sim

def compare_sim_vectors(sim_exp, sim_ptm):
  from scipy.spatial.distance import pdist, squareform
  import numpy as np

  # combine similarity vectors into matrix
  inst_mat = np.vstack((sim_exp, sim_ptm))

  sim_data = 1 - pdist(inst_mat, metric='cosine')

  return sim_data

  # calculate PTM sim mat of cell lines
  ########################################

  # I need to check whether a pairwise complete function exists
  # I will use clustergrammer to do normalization etc.



main()