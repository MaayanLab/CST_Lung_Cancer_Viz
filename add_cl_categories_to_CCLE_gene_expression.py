def main():
  '''
  This will add cell line category information (including plexes and
  gene-expression groups to the gene expression data from CCLE)
  '''
  from clustergrammer import Network
  net = Network()

  # load original CCLE gene expression data for CST lung cancer cell lines
  filename = 'CCLE_gene_expression/CCLE_NSCLC_all_genes.txt'
  f = open(filename, 'r')
  lines = f.readlines()
  f.close()

  # load cell line info
  cl_info = net.load_json_to_dict('cell_line_info/cell_line_muts.json')

  # write to new file
  new_file = 'CCLE_gene_expression/CCLE_NSCLC_cats_all_genes.txt'
  fw = open(new_file, 'w')

  fw.close()

main()