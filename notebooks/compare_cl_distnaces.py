def main():
  '''
  This script will make cell-line by cell-line distnace vectors and using the
  gene expression data (with and withouth gene-zscoring) and PTM data. I'll
  then check how different PTM data processing methods (normalization/filtering)
  affect the distances (and/or similarities) between all cell line pairs in
  gene-expression space

  '''

  from clustergrammer import Network
  from scipy.spatial.distance import pdist, squareform
  import numpy as np

  net = Network()

  filename = '../CCLE_gene_expression/CCLE_NSCLC_all_genes.txt'

  net.load_file(filename)

  tmp_df = net.dat_to_df()

  df = tmp_df['mat']

  print(df.shape)

  col_names = df.columns.tolist()

  print(col_names)

  # transpose to calc distance matrix of columns
  df = df.transpose()

  # calculate the similarity of cell line data based on gene expression
  dm = 1 - pdist(df, metric='cosine')

  print(dm.shape)

  tmp = np.vstack((dm, dm))

  print(tmp.shape)
  print(tmp)

  # dm_mat = squareform(dm)

  compare_dist = 1 - pdist(tmp, metric='cosine')

  print(compare_dist)

  # convert to squareform matrix to check distance matrix size

  # calculate gene-exp sim mat of cell lines
  ###############################################
  # I am using similarity rather than distance because I am more interested in
  # optimizing cell-line similarity than distances

  # sanity check
  # make matrix of two gene-exp sim mat reduced matrices
  # calculate similarity

  # calculate PTM sim mat of cell lines
  ########################################

  # construct PTM matrix with
  # 1) averaged duplicate cell lines run in multiple plexes
  # 2) the same ordering for the cell lines in the columns (alpha) as gene-exp





main()