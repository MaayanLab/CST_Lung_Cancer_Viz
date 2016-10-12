def main():
  '''
  This will pre-calculate different normalizations/filterings of the CCLE
  comparable PTM data and CCLE gene-expression data.
  '''

  data_types = ['exp', 'ptm']

  for inst_type in data_types:
    precalc_processed_versions(inst_type)

def precalc_processed_versions(inst_type):
  from copy import deepcopy
  from clustergrammer import Network

  data_file_names = {
    'exp':'../CCLE_gene_expression/CCLE_NSCLC_all_genes.txt',
    'ptm':'../lung_cellline_3_1_16/lung_cl_all_ptm/all_ptm_ratios_CCLE_cl.tsv'
  }

  filename = data_file_names[inst_type]

  norms = [
  'none', 'row-zscore', 'col-qn', 'col-zscore',
  'col-qn_row-zscore', 'col-zscore_row-zscore'
  ]

  filters = ['filter_'+i for i in norms]

  all_proc = norms + filters

  for inst_filt in all_proc:

    print('process: ' + inst_filt)

    # load data into network so that norm/filtering can be easily done
    ######################################################################
    net = deepcopy(Network())
    net.load_file(filename)

    # perform normalizations and filters
    #######################################
    run_proc = inst_filt.split('_')

    print('\n-------------')
    for i in range(len(run_proc)):
      inst_proc = run_proc[i]
      proc_num = i + 1
      print( str(proc_num) + ': ' + inst_proc)


    # export dataframe
    ######################
    net.swap_nan_for_zero()
    tmp_df = net.dat_to_df()
    df = tmp_df['mat']

    print(df.shape)

    # write to file
    #################

main()