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

  fw_st = open('txt/lung_CL_phos_ratios_ST.tsv', 'w')
  fw_y = open('txt/lung_CL_phos_ratios_Y.tsv', 'w')

  for i in range(len(lines)):

    inst_line = lines[i]

    if i == 0:
      fw_st.write(inst_line)
      fw_y.write(inst_line)

    else:

      split_line = inst_line.strip().split('\t')

      inst_ptm = split_line[0].split('_')

      # print(len(inst_ptm))

      if len(inst_ptm) < 2:
        print(inst_ptm)

      # print(split_line[0])
      # print(inst_ptm)


  fw_st.close()
  fw_y.close()

main()