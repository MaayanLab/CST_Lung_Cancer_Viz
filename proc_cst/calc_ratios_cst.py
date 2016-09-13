def main(ptm_type):
  '''
  calculate ratio of lung cancer cell lines to non-cancerous cell lines
  do this for all ptm types
  '''
  import numpy as np
  for file_num in [1,2,3,4,5,6,7,8,9]:

    inst_filename = 'lung_cellline_3_1_16/lung_cellline_'+ ptm_type +'/'+\
      'lung_cellline_TMT' + str(file_num) + '_'+ ptm_type +'_proc.tsv'

    f = open(inst_filename,'r')
    lines = f.readlines()
    f.close()

    new_filename = 'lung_cellline_3_1_16/lung_cellline_'+ ptm_type +'/'+\
      'lung_cellline_TMT' + str(file_num)+'_'+ ptm_type +'_ratios.tsv'

    fw = open(new_filename, 'w')

    cell_lines = []

    num_cl_per_plex = 5

    for i in range(len(lines)):

      inst_line = lines[i].strip().split('\t')

      # process col titles
      if i == 0:
        col_title_line = inst_line

        fw.write('\t')

        for col_num in range(len(col_title_line)):

          col_name = col_title_line[col_num]

          if col_num > 0:
            print(col_name)
            inst_cl = col_name.split(': ')[1]
            cell_lines.append(inst_cl)

          if col_num > 1:
            fw.write(inst_cl+'\t')

        fw.write('\n')

      # process data
      if i > 0:

        inst_ptm = inst_line[0]

        control_val = float(inst_line[1])

        ratios = []

        inst_data = inst_line[2:]

        # fw.write(inst_ptm+'_'+str(i)+'\t')
        fw.write(inst_ptm+'\t')

        for j in range(len(inst_data)):

          data_point = float(inst_data[j])

          if control_val !=0:
            inst_ratio = data_point/control_val

            if inst_ratio !='nan' and inst_ratio !=0:
              inst_ratio = np.log2(inst_ratio)

            if inst_ratio > 100:
              inst_ratio = 100
          else:
            inst_ratio = 'nan'

          inst_ratio = str(inst_ratio)

          fw.write(inst_ratio+'\t')

        fw.write('\n')

    fw.close()