def main():

  # only run once (need to add pairwise-complete versions)
  make_processed_versions()

  # generate plex matrix (rows are plexes and columns are cell lines)
  # some cell lines have two plexes
  make_plex_matrix()

def make_processed_versions():
  '''
  This will pre-calculate different normalizations/filterings of the CCLE
  comparable PTM data and CCLE gene-expression data.
  '''

  data_types = ['exp', 'ptm', 'ptm45']

  for inst_type in data_types:
    precalc_processed_versions(inst_type)

def precalc_processed_versions(inst_type):
  from copy import deepcopy
  from clustergrammer import Network

  data_file_names = {
    'exp':'../CCLE_gene_expression/CCLE_NSCLC_all_genes.txt',
    'ptm':'../lung_cellline_3_1_16/lung_cl_all_ptm/all_ptm_ratios_CCLE_cl.tsv',
    'ptm45':'../lung_cellline_3_1_16/lung_cl_all_ptm/all_ptm_ratios.tsv',
  }

  filename = data_file_names[inst_type]

  norms = [
  'none', 'row-zscore', 'col-qn', 'col-zscore',
  'col-qn_row-zscore', 'col-zscore_row-zscore'
  ]

  # only add filter for PTM data
  if inst_type == 'ptm' or inst_type == 'ptm45':
    filter_before = ['filter_'+i for i in norms]
    filter_after = [i+'_filter' for i in norms]
    all_proc = norms + filter_before + filter_after
  else:
    all_proc = norms

  for inst_filt in all_proc:

    print('\n\n-- '+ inst_type +': all processes: ' + inst_filt)
    print('----------------------------')

    # load data into network so that norm/filtering can be easily done
    ####################################################################
    net = deepcopy(Network())
    net.load_file(filename)

    # perform normalizations and filters
    #######################################
    run_proc = inst_filt.split('_')

    for i in range(len(run_proc)):
      inst_proc = run_proc[i]
      proc_num = i + 1
      print( str(proc_num) + ': ' + inst_proc)

      print('**********')
      print(inst_type)
      print('**********')

      net = process_net(net, inst_proc, inst_type)

    # export dataframe (keep nans)
    ###############################
    tmp_df = net.dat_to_df()
    df = tmp_df['mat']

    print(df.shape)

    # write to file
    #################
    inst_filename = '../lung_cellline_3_1_16/lung_cl_all_ptm/'+\
                    'precalc_processed/' + inst_type + '_' + inst_filt + '.txt'
    df.to_csv(inst_filename, sep='\t', na_rep='nan')

def process_net(net, inst_proc, inst_type):


  print('processing network: ' + inst_proc  + ' ' + inst_type + '\n********************')

  if inst_proc == 'row-zscore':
    net.normalize(axis='row', norm_type='zscore')
  elif inst_proc == 'col-qn':
    net.normalize(axis='col', norm_type='qn')
  elif inst_proc == 'col-zscore':
    net.normalize(axis='col', norm_type='zscore')
  elif inst_proc == 'filter':
    if inst_type == 'ptm':
      # this removes ptms with missing data
      net.filter_threshold('row', threshold=0, num_occur=37)
    else:
      # this removes ptms with missing data
      net.filter_threshold('row', threshold=0, num_occur=45 )

  return net

def make_plex_matrix():
  '''
  Make a cell line matrix with plex rows and cell line columns.
  This will be used as a negative control that should show worsening correlation
  as data is normalized/filtered.
  '''
  import numpy as np
  import pandas as pd
  from clustergrammer import Network

  # load cl_info
  net = Network()
  cl_info = net.load_json_to_dict('../cell_line_info/cell_line_info_dict.json')

  # load cell line expression
  net.load_file('../CCLE_gene_expression/CCLE_NSCLC_all_genes.txt')
  tmp_df = net.dat_to_df()
  df = tmp_df['mat']

  cols = df.columns.tolist()

  rows = range(9)
  rows = [i+1 for i in rows]
  print(rows)

  mat = np.zeros((len(rows), len(cols)))

  for inst_col in cols:

    for inst_cl in cl_info:

      if inst_col in inst_cl:
        inst_plex = int(cl_info[inst_cl]['Plex'])

        if inst_plex != -1:
          # print(inst_col + ' in ' + inst_cl + ': ' + str(inst_plex))

          row_index = rows.index(inst_plex)
          col_index = cols.index(inst_col)

          mat[row_index, col_index] = 1


  df_plex = pd.DataFrame(data=mat, columns=cols, index=rows)

  filename = '../lung_cellline_3_1_16/lung_cl_all_ptm/precalc_processed/' + \
            'exp-plex.txt'
  df_plex.to_csv(filename, sep='\t')

main()