def main():
  '''
  This script will make cell-line by cell-line distnace vectors and using the
  gene expression data (with and withouth gene-zscoring) and PTM data. I'll
  then check how different PTM data processing methods (normalization/filtering)
  affect the distances (and/or similarities) between all cell line pairs in
  gene-expression space
  '''

  from clustergrammer import Network

  net = Network()

  filename = '../CCLE_gene_expression/CCLE_NSCLC_all_genes.txt'

  net.load_file(filename)

  tmp_df = net.dat_to_df()

  df = tmp_df['mat']

  print(df.shape)

main()