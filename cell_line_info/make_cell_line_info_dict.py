def main():
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
        cl = cell_lines[j]
        cl_info[cl] = {}

  save_dict_to_json(cl_info)


def save_dict_to_json(inst_dict):
  from clustergrammer import Network
  net = Network()

  net.save_dict_to_json(inst_dict, 'cell_line_info_dict.json', indent='indent')

main()
