def main():

  # make_tuple_version_of_lung_CL_phos()

  separate_ST_Y()

def make_tuple_version_of_lung_CL_phos():
  from make_clustergrammer import Network

  net = Network()

  net.load_file('txt/lung_CL_phos_ratios_non_tuple.tsv')

  net.write_matrix_to_tsv ('txt/lung_CL_phos_ratios.tsv')

def separate_ST_Y():
  # open lung cell line ratio phos
  f = open('txt/lung_CL_phos_ratios.tsv')
  lines = f.readlines()
  f.close()

  for i in range(len(lines)):

    inst_line = lines[i].strip().split('\t')

    if i == 0:
      print(inst_line)
      print(len(inst_line))



main()