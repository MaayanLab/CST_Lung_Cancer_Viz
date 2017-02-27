def main():
  '''
  This will make a PTM matrix that is compatible with (e.g. easily comparable)
  with the gene expression matrix: 1) average duplicate cell line measurements,
  2) only include the 37 cell lines that are found in CCLE, 3) put cell lines
  into the same order as the gene expression data.
  '''
  import pandas as pd

  starting_data = 'all_ptm_ratios'
  average_plex_runs(starting_data)

  keep_only_CCLE_cl(starting_data)

def average_plex_runs(starting_data):
  import pandas as pd
  from clustergrammer import Network
  from copy import deepcopy

  print('averaging plex runs and saving to new tsv')

  # load all ptm ratios
  filename_all_ptm = '../lung_cellline_3_1_16/'+ \
                     'lung_cl_all_ptm/'+starting_data+'.tsv'

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
                      starting_data +'_uni_cl.tsv'

  df_uni_cl.to_csv(filename_unique_cl, sep='\t', na_rep='nan')

def keep_only_CCLE_cl(starting_data):
  from clustergrammer import Network
  from copy import deepcopy
  # load all ptm ratios
  filename_all_ptm = '../lung_cellline_3_1_16/lung_cl_all_ptm/'+\
                     starting_data +'_uni_cl.tsv'

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
                     starting_data +'_CCLE_cl.tsv'

  df_ptm.to_csv(filename_CCLE_cl, sep='\t', na_rep='nan')

main()