def main():

  cl_info = add_hist_plex_gender_exp()

  cl_info = add_mutations(cl_info)

  save_dict_to_json(cl_info)

def add_mutations(cl_info):
  print('add mutations\n')
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
