def main():
  '''
  This will make a PTM matrix that is compatible with (e.g. easily comparable)
  with the gene expression matrix: 1) average duplicate cell line measurements,
  2) only include the 37 cell lines that are found in CCLE, 3) put cell lines
  into the same order as the gene expression data.
  '''
  import pandas as pd

  combine_and_save_ptm()

  average_plex_runs()

  keep_only_CCLE_cl()

def combine_and_save_ptm():
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

def average_plex_runs():
  import pandas as pd
  from clustergrammer import Network
  from copy import deepcopy

  print('averaging plex runs and saving to new tsv')

  # load all ptm ratios
  filename_all_ptm = '../lung_cellline_3_1_16/'+ \
                     'lung_cl_all_ptm/all_ptm_ratios.tsv'

  net_ptm = deepcopy(Network())
  net_ptm.load_file(filename_all_ptm)

  tmp_df = net_ptm.dat_to_df()
  df_ptm = tmp_df['mat']

  ptm_cols = df_ptm.columns.tolist()

  print('shape of all ptms')
  print(df_ptm.shape)
  print('\n')

  cl_with_duplicates = {}

  for inst_plex in ptm_cols:
    if '_plex' in inst_plex:

      inst_cl = inst_plex.split('_plex')[0]

      if inst_cl not in cl_with_duplicates:
        cl_with_duplicates[inst_cl] = []

      cl_with_duplicates[inst_cl].append(inst_plex)

  # print(cl_with_duplicates)

  df_add = deepcopy(df_ptm)

  # merge data
  for dup_cl in cl_with_duplicates:
    print(dup_cl)
    inst_plexes = cl_with_duplicates[dup_cl]

    df_dup = deepcopy(df_ptm[inst_plexes])

    print(df_dup.shape)

    # calc mean of col vectors
    df_mean = df_dup.mean(axis=1)

    print(df_mean.shape)

    # add series to df
    df_add[dup_cl] = df_mean

  print(df_add.shape)

  cl_add = df_add.columns.tolist()
  print(cl_add)

  cl_keep = []
  for inst_cl in cl_add:
    if 'plex' not in inst_cl:
      cl_keep.append(inst_cl)

  print('---------')
  print('remove duplicate plex cell lines from unique_cl version')
  cl_keep.sort()

  print(cl_keep)

  df_uni_cl = deepcopy(df_add[cl_keep])

  filename_unique_cl = '../lung_cellline_3_1_16/lung_cl_all_ptm/'+\
                      'all_ptm_ratios_uni_cl.tsv'

  df_uni_cl.to_csv(filename_unique_cl, sep='\t')

def keep_only_CCLE_cl():
  from clustergrammer import Network
  from copy import deepcopy
  # load all ptm ratios
  filename_all_ptm = '../lung_cellline_3_1_16/lung_cl_all_ptm/'+\
                     'all_ptm_ratios_uni_cl.tsv'

  net_ptm = deepcopy(Network())
  net_ptm.load_file(filename_all_ptm)

  tmp_df = net_ptm.dat_to_df()
  df_ptm_all_cl = tmp_df['mat']

  print('shape of all ptms')
  print(df_ptm_all_cl.shape)

  # get gene-exp cell lines
  net_exp = deepcopy(Network())
  net_exp.load_file('../CCLE_gene_expression/CCLE_NSCLC_all_genes.txt')

  tmp_df = net_exp.dat_to_df()
  df_exp = tmp_df['mat']

  print('shape of exp')
  print(df_exp.shape)

  # only keep gene expression cell lines
  #######################################
  cl_exp = df_exp.columns.tolist()
  df_ptm = df_ptm_all_cl[cl_exp]
  print(df_ptm.shape)

  filename_CCLE_cl = '../lung_cellline_3_1_16/lung_cl_all_ptm/'+\
                     'all_ptm_ratios_CCLE_cl.tsv'

  df_ptm.to_csv(filename_CCLE_cl, sep='\t')

main()