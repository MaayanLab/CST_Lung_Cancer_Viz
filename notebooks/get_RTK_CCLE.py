def main():
  import numpy as np
  import pandas as pd
  from clustergrammer import Network

  rtk_list = load_rtks()

  net = Network()
  net.load_file('txt/tmp_cst_drug_treat_cl.txt')
  df_dict = net.dat_to_df()

  inst_df = df_dict['mat']

  inst_df = inst_df.ix[rtk_list]

  inst_df.to_csv('txt/RTK_exp_in_drug_treat_cl.txt', sep='\t')


def load_rtks():
  f = open('RTK_list.txt')
  lines = f.readlines()
  f.close()

  rtk_list = []

  for inst_line in lines:
    inst_line = inst_line.strip()
    rtk_list.append(inst_line)

  return rtk_list


main()