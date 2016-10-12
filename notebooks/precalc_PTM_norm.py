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

    print('\n\n-- all processes: ' + inst_filt)
    print('------------------------\n')

    # load data into network so that norm/filtering can be easily done
    ######################################################################
    net = deepcopy(Network())
    net.load_file(filename)

    # perform normalizations and filters
    #######################################
    run_proc = inst_filt.split('_')

    for i in range(len(run_proc)):
      inst_proc = run_proc[i]
      proc_num = i + 1
      print( str(proc_num) + ': ' + inst_proc)

      net = process_net(net, inst_proc)

    # export dataframe
    ######################
    net.swap_nan_for_zero()
    tmp_df = net.dat_to_df()
    df = tmp_df['mat']

    print(df.shape)

    # write to file
    #################
    inst_filename = '../lung_cellline_3_1_16/lung_cl_all_ptm/'+\
                    'precalc_processed/' + inst_type + '_' + inst_filt + '.txt'
    df.to_csv(inst_filename, sep='\t', na_rep='nan')

def process_net(net, inst_proc):

  print('processing network: ' + inst_proc + '\n********************')

  if inst_proc == 'row-zscore':
    net.normalize(axis='row', norm_type='zscore')
  elif inst_proc == 'col-qn':
    net.normalize(axis='col', norm_type='qn')
  elif inst_proc == 'col-zscore':
    net.normalize(axis='col', norm_type='zscore')
  elif inst_proc == 'filter':
    # this removes ptms with missing data
    net.filter_threshold('row', threshold=0, num_occur=37)

  return net

main()