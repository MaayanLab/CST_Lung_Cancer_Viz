def main():

  cl_info = add_hist_plex_gender_exp()

  cl_info = add_mutations(cl_info)

  cl_info = make_non_plex_cl_version(cl_info)

  cl_info = add_corrected_expression_clusters(cl_info)

  save_dict_to_json(cl_info)

def add_corrected_expression_clusters(cl_info):
  '''
  add the gene expression clusters that have the correct cell line aliases
  (e.g. Calu-3 vs CALU3)
  '''

  # load gene expression clusters from gmteseque tsv file
  f = open('expression_clusters.txt')
  lines = f.readlines()
  exp_groups = {}
  for i in range(len(lines)):

    group_num = i + 1
    inst_line = lines[i]
    exp_groups[group_num] = inst_line.strip().split('\t')
  f.close()


  for inst_cl in cl_info:

    # initialize at -1 group
    cl_group = -1

    # remove plex name if necessary
    if '_plex_' in inst_cl:
      simple_cl = inst_cl.split('_')[0]
    else:
      simple_cl = inst_cl

    # check all expression groups and add number
    for inst_group in exp_groups:

      # check with simple cl name but save using full name (may include plex)
      if simple_cl in exp_groups[inst_group]:
        cl_group = inst_group

    # save group to cl_info
    cl_info[inst_cl]['Exp-group'] = cl_group


  return cl_info

def make_non_plex_cl_version(cl_info):
  from copy import deepcopy

  for inst_cl in ['H209_plex_1', 'H2073_plex_9', 'H1437_plex_8']:

    cl_name = inst_cl.split('_')[0]

    # make a non-plex copy
    cl_info[cl_name] = deepcopy(cl_info[inst_cl])
    # set the plex to -1 since it is present in two plexes
    cl_info[cl_name]['Plex'] =  -1

  return cl_info

def add_mutations(cl_info):
  print('add mutations\n')

  from clustergrammer import Network
  net = Network()
  old_cl_info = net.load_json_to_dict('cell_line_muts.json')

  cl_muts = old_cl_info['muts']

  for inst_cl in cl_info:

    # remove plex name if necessary
    if '_plex_' in inst_cl:
      simple_cl = inst_cl.split('_')[0]
    else:
      simple_cl = inst_cl

    for inst_mut in cl_muts:
      mutated_cls = cl_muts[inst_mut]

      if simple_cl in mutated_cls:
        has_mut = 'true'
      else:
        has_mut = 'false'

      mutation_title = 'mut-'+inst_mut

      # use the original long cell line name (with possible plex)
      cl_info[inst_cl][mutation_title] = has_mut

  return cl_info

def add_hist_plex_gender_exp():
  print('add hist, plex, gender, and expression-groups\n')
  filename = 'CST_cell_line_info.txt'
  f = open(filename, 'r')
  lines = f.readlines()
  f.close()

  cl_info = {}
  for i in range(len(lines)):
    inst_line = lines[i]
    inst_line = inst_line.strip().split('\t')

    # initialize cl_info dict
    if i == 0:
      cell_lines = inst_line

      for inst_cl in cell_lines:
        cl_info[inst_cl] = {}

    else:

      # populate with cell line information from CST_cell_line_info.txt
      for j in range(len(inst_line)):
        inst_cat = inst_line[j]

        cat_title = inst_cat.split(': ')[0]
        cat_state = inst_cat.split(': ')[1]

        cl = cell_lines[j]


        # handle expression groups differently (combine into single cat)
        if 'exp-group' not in cat_title:
          if 'Plex-' in cat_state:
            cat_state = cat_state.split('-')[1]
          cl_info[cl][cat_title] = cat_state

        else:


          # combine into single category
          if cat_state == 'true':
            group_num = cat_title.split('-')[2]
            cl_info[cl]['Exp-group'] = group_num

  return cl_info

def save_dict_to_json(inst_dict):
  print('save to cell_line_info_dict.json\n')
  from clustergrammer import Network
  net = Network()

  net.save_dict_to_json(inst_dict, 'cell_line_info_dict.json', indent='indent')

main()
