def main():
  '''
  This will add cell line category information (including plexes and
  gene-expression groups to the gene expression data from CCLE)
  '''

  filename = 'CCLE_gene_expression/CCLE_NSCLC_all_genes.txt'
  f = open(filename, 'r')
  lines = f.readlines()
  f.close()

main()