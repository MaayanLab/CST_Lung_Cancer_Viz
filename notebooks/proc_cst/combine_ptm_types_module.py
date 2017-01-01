def main():

  # combine all ptm data into single dataframe
  #################################################################
  # use simple col names

  from clustergrammer import Network
  import pandas as pd
  from copy import deepcopy

  ptm_data = {
  'phos': '../lung_cellline_3_1_16/lung_cellline_phospho/lung_cellline_TMT_phospho_combined_ratios.tsv',
  'act': '../lung_cellline_3_1_16/lung_cellline_Ack/lung_cellline_TMT_Ack_combined_ratios.tsv',
  'met_arg': '../lung_cellline_3_1_16/lung_cellline_Rme1/lung_cellline_TMT_Rme1_combined_ratios.tsv',
  'met_lys': '../lung_cellline_3_1_16/lung_cellline_Kme1/lung_cellline_TMT_Kme1_combined_ratios.tsv'
  }

  df_all = pd.DataFrame()

  for inst_type in ptm_data:

    net = deepcopy(Network())
    filename = ptm_data[inst_type]

    net.load_file(filename)

    tmp_df = net.dat_to_df()
    inst_df = tmp_df['mat']

    col_names = inst_df.columns.tolist()

    print('\ninst_type ' + inst_type)
    print(inst_df.shape)

    df_all = pd.concat([df_all, inst_df], axis=0)

  filename_all_ptm = '../lung_cellline_3_1_16/lung_cl_all_ptm/all_ptm_ratios.tsv'

  print('\nshape of mat')
  print(df_all.shape)

  df_all.to_csv(filename_all_ptm, sep='\t', na_rep='nan')

  print('\ncheck if rows are unique')
  all_rows = df_all.index.tolist()
  print(len(all_rows))
  all_rows = list(set(all_rows))
  print(len(all_rows))

  print('\nnumber of cell lines ')
  print(len(df_all.columns.tolist()))

  # I am manually removing trailing tabs need to improve this
  ################################################################
