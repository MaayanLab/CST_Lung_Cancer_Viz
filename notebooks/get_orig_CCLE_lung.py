def main():

  from clustergrammer import Network

  # load CCLE cell lines
  filename = '../lung_cellline_3_1_16/lung_cl_all_ptm/precalc_processed/ccle_cl_names.txt'
  f = open(filename, 'r')
  lines = f.readlines()
  f.close()

  cl_names = []
  for inst_line in lines:
    inst_line = inst_line.strip()
    cl_names.append(inst_line)


  filename = '../lung_cellline_3_1_16/lung_cl_all_ptm/precalc_processed/CCLE_lung.txt'

  net = Network()

  net.load_file(filename)

  ccle_lung = net.export_df()

  cols = ccle_lung.columns.tolist()

  # simplify cols, disguard meta-data
  ######################################
  simple_cols = []

  for inst_col in cols:
    proc_col = inst_col[0].split(': ')[1].replace('NCI','')

    if 'CALU' in proc_col:
      proc_col = proc_col.replace('CALU', 'Calu-')

    if 'LOU' in proc_col:
      proc_col = proc_col.replace('LOU', 'Lou-')

    if 'CAL' in proc_col:
      proc_col = proc_col.replace('CAL', 'CAL-')

    simple_cols.append(proc_col)

  ccle_lung.columns = simple_cols

  cols = ccle_lung.columns.tolist()

  found_cols = []

  for inst_col in cols:
    if inst_col in cl_names:
      found_cols.append(inst_col)

  # found all cell lines
  print('found ' + str(len(found_cols)))


  # save subset of cell lines that are also found in the CST PTM data
  ccle_cst_lung = ccle_lung[cl_names]

  save_filename = '../lung_cellline_3_1_16/lung_cl_all_ptm/precalc_processed/CCLE_CST_lung.txt'
  ccle_cst_lung.to_csv(save_filename, sep='\t')



main()