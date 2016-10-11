def main():
  '''
  This script will make cell-line by cell-line distnace vectors and using the
  gene expression data (with and withouth gene-zscoring) and PTM data. I'll
  then check how different PTM data processing methods (normalization/filtering)
  affect the distances (and/or similarities) between all cell line pairs in
  gene-expression space

  '''

  # calculte similarity vector based on expression data
  sim_exp = calc_cl_exp_sim()

  # calculate similarity vector based on ptm data
  save_gene_exp_compatible_ptm_data()

  # compare similarity vectors based on expression and ptm data
  sim_data = compare_sim_vectors(sim_exp, sim_exp)

  print(sim_data)

def save_gene_exp_compatible_ptm_data():
  '''
  This will make a PTM matrix that is compatible with (e.g. easily comparable)
  with the gene expression matrix: 1) average duplicate cell line measurements,
  2) only include the 37 cell lines that are found in CCLE, 3) put cell lines
  into the same order as the gene expression data.
  '''

  # sort names in place
  # col_names.sort()

  # check that the cell lines are in the same order in both exp and PTM data

  pass

def calc_cl_exp_sim():

  from scipy.spatial.distance import pdist, squareform
  from clustergrammer import Network

  net = Network()

  filename = '../CCLE_gene_expression/CCLE_NSCLC_all_genes.txt'

  net.load_file(filename)

  tmp_df = net.dat_to_df()

  df = tmp_df['mat']

  col_names = df.columns.tolist()



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