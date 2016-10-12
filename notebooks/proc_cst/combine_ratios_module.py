def main(ptm_type):

  all_cl, all_ptms = get_all_cl_and_ptm(ptm_type)

  combine_all_plexes(ptm_type, all_cl, all_ptms)

def combine_all_plexes(ptm_type, all_cl, all_ptms):

  print('\ncombine all plex data into one file')

  import scipy
  import numpy as np
  import pandas as pd
  from clustergrammer import Network

  mat = scipy.zeros([len(all_ptms), len(all_cl)])

  mat[:] = np.nan

  rows = sorted(all_ptms)
  cols = sorted(all_cl)

  exp_clusters = load_expression_clusters()

  print(exp_clusters)

  for file_num in [1,2,3,4,5,6,7,8,9]:

    inst_filename = '../lung_cellline_3_1_16/lung_cellline_'+ ptm_type +'/'+\
      'lung_cellline_TMT' + str(file_num) + '_'+ ptm_type +'_ratios.tsv'

    print(inst_filename)
    f = open(inst_filename, 'r')
    lines = f.readlines()
    f.close()

    plex_ptms = []

    for i in range(len(lines)):

      inst_line = lines[i].strip().split('\t')

      if i == 0:
        plex_cl = inst_line

      if i > 0:
        inst_ptm = inst_line[0]
        inst_data = inst_line[1:]

        if len(list(set(inst_data))) > 1:
          row_index = rows.index(inst_ptm)

          # save data to matrix
          for j in range(len(inst_data)):

            inst_cl = plex_cl[j] + '_' + str(file_num)

            col_index = cols.index(inst_cl)

            inst_val = float(inst_data[j])

            mat[row_index, col_index] = inst_val

  print(cols)

  # only keep plex info for duplicate cls
  duplicate_cl = ['H1437', 'H2073', 'H209']
  cols_plex = []
  for inst_col in cols:
    if inst_col.split('_')[0] in duplicate_cl:
      new_col = inst_col.split('_')[0] + '_plex_' + inst_col.split('_')[1]
    else:
      new_col = inst_col.split('_')[0]
    cols_plex.append(new_col)

  inst_df = pd.DataFrame(data=mat, columns=cols_plex, index=rows)

  new_filename = '../lung_cellline_3_1_16/lung_cellline_'+ ptm_type +'/'+\
    'lung_cellline_TMT_'+ptm_type+'_combined_ratios.tsv'
  inst_df.to_csv(new_filename, sep='\t', na_rep='nan')

  # # write matrix to new file
  # new_filename = '../lung_cellline_3_1_16/lung_cellline_'+ ptm_type +'/'+\
  #   'lung_cellline_TMT_'+ptm_type+'_combined_ratios.tsv'
  # fw = open(new_filename, 'w')

  # # write first line with col labels
  # ########################################
  # duplicate_cl = ['H1437', 'H2073', 'H209']
  # fw.write('\t')
  # for inst_col in cols:
  #   if inst_col.split('_')[0] in duplicate_cl:
  #     inst_name = inst_col.split('_')[0] + '_plex_' + inst_col.split('_')[1]
  #   else:
  #     inst_name = inst_col.split('_')[0]
  #   fw.write(inst_name+'\t')
  # fw.write('\n')

  # net = Network()

  # # load cell line info
  # cl_info = net.load_json_to_dict('../cell_line_info/cell_line_muts.json')

  # # histology
  # ########################################
  # fw.write('\t')
  # for inst_col in cols:
  #   tmp_col = inst_col.split('_')[0]

  #   if tmp_col in cl_info['hist']['NSCLC']:
  #     inst_hist = 'Histology: NSCLC'
  #   elif tmp_col in cl_info['hist']['SCLC']:
  #     inst_hist = 'Histology: SCLC'

  #   fw.write(inst_hist+'\t')
  # fw.write('\n')

  # # plex
  # ########################################
  # fw.write('\t')
  # for inst_col in cols:
  #   tmp_plex = inst_col.split('_')[1]

  #   inst_plex = 'Plex: Plex-'+ tmp_plex

  #   fw.write(inst_plex+'\t')
  # fw.write('\n')

  # # gender
  # ########################################
  # fw.write('\t')
  # for inst_col in cols:
  #   tmp_col = inst_col.split('_')[0]

  #   if tmp_col in cl_info['gender']['M']:
  #     inst_gender = 'Gender: M'
  #   elif tmp_col in cl_info['gender']['F']:
  #     inst_gender = 'Gender: F'
  #   else:
  #     inst_gender = 'Gender: N.A.'

  #   fw.write(inst_gender+'\t')
  # fw.write('\n')

  # # # mutations
  # # ########################################
  # # all_muts = cl_info['muts'].keys()
  # # for inst_mut in all_muts:

  # #   fw.write('\t')
  # #   for inst_col in cols:
  # #     tmp_col = inst_col.split('_')[0]

  # #     if tmp_col in cl_info['muts'][inst_mut]:
  # #       is_mut = inst_mut+': '+'mutated'
  # #     else:
  # #       is_mut = inst_mut +': '+'Not mutated'

  # #     fw.write(is_mut+'\t')
  # #   fw.write('\n')

  # # exp groups
  # ########################################
  # all_exp_groups = exp_clusters.keys()
  # for inst_group in all_exp_groups:

  #   fw.write('\t')
  #   for inst_col in cols:
  #     tmp_col = inst_col.split('_')[0]

  #     if tmp_col in exp_clusters[inst_group]:
  #       is_mut = inst_group+': '+'true'
  #     else:
  #       is_mut = inst_group +': '+'false'

  #     fw.write(is_mut+'\t')
  #   fw.write('\n')

  # # write rows
  # for i in range(len(rows)):
  #   inst_row = rows[i]
  #   fw.write(inst_row+'\t')

  #   inst_data = mat[i,:]

  #   for inst_val in inst_data:
  #     fw.write(str(inst_val)+'\t')

  #   fw.write('\n')

  # fw.close()

def get_all_cl_and_ptm(ptm_type):
  print('combine '+ ptm_type +' data from all plexes')

  import numpy as numpy

  all_cl = []

  all_ptms = []

  for file_num in [1,2,3,4,5,6,7,8,9]:

    inst_filename = '../lung_cellline_3_1_16/lung_cellline_'+ ptm_type +'/'+\
      'lung_cellline_TMT' + str(file_num) + '_'+ ptm_type +'_ratios.tsv'

    f = open(inst_filename, 'r')
    lines = f.readlines()
    f.close()

    plex_ptms = []

    for i in range(len(lines)):

      inst_line = lines[i].strip().split('\t')

      if i == 0:

        # no need to skip first col because of strip
        inst_cl = inst_line

        inst_cl = [x+'_'+str(file_num) for x in inst_cl]

        all_cl.extend(inst_cl)

      if i > 0:

        if 'H1417' in inst_line[0]:
          print('i '+str(i))
          print(i)
          print(inst_line)

        inst_data = inst_line[1:]

        check_data = list(set(inst_data))

        if len(check_data) > 1:

          plex_ptms.append(inst_line[0])

          if '_' not in inst_line[0]:
            print(inst_filename + '\t' + inst_line[0])

    plex_ptms = list(set(plex_ptms))

    all_ptms.extend(plex_ptms)

  all_ptms = list(set(all_ptms))
  all_cl = list(set(all_cl))

  print('ptms: '+str(len(all_ptms)))
  print('cl: '+str(len(all_cl)))

  if 'H1299' in all_ptms:
    print('********* found H1299 in all_ptms')

  return all_cl, all_ptms


def load_expression_clusters():

  f = open('../cell_line_info/expression_clusters.txt')
  lines = f.readlines()

  exp_clusters = {}

  for i in range(len(lines)):

    inst_group = 'exp-group-' + str(i + 1)

    inst_line = lines[i]

    inst_list = inst_line.strip().split('\t')

    exp_clusters[inst_group] = inst_list

  return exp_clusters

  f.close()