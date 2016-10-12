def main():
  '''
  This script will make cell-line by cell-line distnace vectors and using the
  gene expression data (with and withouth gene-zscoring) and PTM data. I'll
  then check how different PTM data processing methods (normalization/filtering)
  affect the distances (and/or similarities) between all cell line pairs in
  gene-expression space
  '''

  # calculate similarities of cell lines
  ##########################################
  # # calculte similarity vector based on expression data
  sim_exp = calc_cl_sim(data_type='exp')

  # calculate ptm sim
  sim_ptm = calc_cl_sim(data_type='ptm')
  sim_ptm_norm = calc_cl_sim(data_type='ptm', col_qn=True, row_zscore=True)

  # # compare similarity vectors based on expression and ptm data
  # ###############################################################
  # sim_data = compare_sim_vectors(sim_exp, sim_ptm)
  # print(sim_data)

  print('here')
  import Mantel

  results = Mantel.test(sim_exp, sim_ptm, perms=10000, tail='upper')
  print(results)

  results = Mantel.test(sim_exp, sim_ptm_norm, perms=10000, tail='upper')
  print(results)

def calc_cl_sim(data_type='exp', sum_filter=None, var_filter=None,
                    row_zscore=False, col_qn=False, col_zscore=False.
                    dist_metric='euclidean'):
  '''
  calculate cell line similarity based on data_type (e.g. expression) with
  optional filtering and normalization
  '''

  from scipy.spatial.distance import pdist, squareform
  from clustergrammer import Network

  net = Network()

  all_data = {
    'exp':'../CCLE_gene_expression/CCLE_NSCLC_all_genes.txt',
    'ptm':'../lung_cellline_3_1_16/lung_cellline_TMT_all_ptm_ratios_CCLE_cl.tsv'
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

  net.swap_nan_for_zero()

  tmp_df = net.dat_to_df()

  df = tmp_df['mat']

  col_names = df.columns.tolist()

  print(df.shape)

  # transpose to calc distance matrix of columns
  df = df.transpose()

  # calculate the similarity of cell line data based on gene expression
  sim = 1 - pdist(df, metric=dist_metric)

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