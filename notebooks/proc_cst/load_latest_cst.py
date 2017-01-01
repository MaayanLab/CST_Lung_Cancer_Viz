def save_lung_data_to_tsv(ptm_type):
  import numpy as np

  for file_num in [1,2,3,4,5,6,7,8,9]:
    print( 'initial processing of file ' + str(file_num))
    inst_filename = '../lung_cellline_3_1_16/lung_cellline_'+ ptm_type +'/lung_cellline_TMT'+\
      str(file_num)+'_'+ ptm_type +'.txt'

    f = open(inst_filename,'r')
    lines = f.readlines()
    f.close()

    # 0 protein
    # 1 gene symbol
    # 2 residue
    # 3 (+/-)7 seq
    # 4 ACC#
    # 5 site group ID
    # 6 site ID

    # the number of expected data points
    # num_data = 36
    num_data = 6

    # save_filename = 'treated_cell_12_1_2015/treated_cl_'+inst_file.split('_')[2]+ '.tsv'
    save_filename = inst_filename.replace('.txt','_proc.tsv')

    fw = open(save_filename, 'w')

    col_start = 12
    col_end = 18

    for i in range(len(lines)):

      inst_line = lines[i].strip()
      inst_line = inst_line.split('\t')

      if i == 0:
        tmp_data_names = inst_line[col_start:]

        write_line = 'ptm\t' + '\t'.join(tmp_data_names) + '\n'

        fw.write(write_line)
      if i > 0:
        inst_gene = inst_line[1]
        inst_residue = inst_line[2]
        inst_ptm = inst_gene + '_' + inst_residue

        tmp_data = inst_line[col_start:]

        # add NaNs for missing data
        inst_data = []

        for j in range(num_data):

          # fill in missing values with NaN
          try:
            data_point = tmp_data[j]

            if data_point == '':
              data_point = 'NaN'
          except:
            data_point = 'NaN'

          inst_data.append(data_point)

        # save to file
        write_line = inst_ptm + '\t'  + '\t'.join(inst_data) + '\n'

        fw.write(write_line)

    fw.close()

def calc_1hr_ratios():
  '''
  This will calculate the treatment vs control ratios for 1hr treatment with all
  drugs. I'll start from the treated_cl_phospho.tsv file that has missing values
  filled in with NaNs.
  '''
  import numpy as np

  f = open('treated_cell_12_1_2015/treated_cl_phospho.tsv','r')
  lines = f.readlines()
  f.close()

  inst_line = lines[0].strip().split('\t')

  all_treatments = []
  treatment_index = []
  for i in range(len(inst_line)):
    inst_cl = inst_line[i]

    if 'ontrol' in inst_cl and 'Unstarved' not in inst_cl and 'ontrol2' not in inst_cl:
      print('\n'+str(i)+'\t'+inst_cl)
    elif '1hr' in inst_cl:

      print('\t'+str(i)+'\t'+inst_cl.replace(' 1uM','').replace(':1hr','').replace('NCI-','')+' 1hr')

      # make list of treatments
      inst_treatment = inst_cl.replace(' 1uM','').replace(':1hr',''.replace('NCI-',''))+' 1hr'
      all_treatments.append( inst_treatment )
      treatment_index.append(i)

  print('\n')
  print(all_treatments)
  print(treatment_index)

  treatment_dict = {}
  treatment_dict['2'] = 1
  treatment_dict['4'] = 1
  treatment_dict['6'] = 1
  treatment_dict['9'] = 7
  treatment_dict['14'] = 13
  treatment_dict['21'] = 19
  treatment_dict['27'] = 25
  treatment_dict['33'] = 31

  # layout of data
  #######################
  # 1 A549:control
  #   2 A549 iressa 1hr
  #   4 A549 gleevec 1hr
  #   6 A549 crizotinib 1hr

  # 7 MKN-45:control1
  #   9 MKN-45 crizotinib 1hr

  # 13  NCI-H1703:control
  #   14  H1703 gleevec 1hr

  # 19  NCI-H2228:Control Starved
  #   21  H2228 crizotinib 1hr

  # 25  NCI-H3122:Control Starved
  #   27  H3122 crizotinib 1hr

  # 31  NCI-H3255:control1
  #   33  H3255 iressa 1hr

  # fw = open('treated_cell_12_1_2015/ratios_1hr_phospho.tsv','w')
  # fw = open('treated_cell_12_1_2015/ratios_1hr_phospho_full.tsv','w')
  fw = open('treated_cell_12_1_2015/ratios_1hr_phospho_missing_one.tsv','w')

  max_ratio = 4

  for i in range(len(lines)):
    inst_line = lines[i].strip().split('\t')

    if i == 0:
      fw.write('\t')
      for j in range(len(inst_line)):
        inst_name = inst_line[j]
        inst_name = inst_name.replace(' 1uM','').replace(':1hr','').replace('NCI-','')+' 1hr'
        if j in treatment_index:
          fw.write(inst_name+'\t')

          control_index = treatment_dict[str(j)]
          print('use control '+ str(control_index) + ' for ' + str(j))

      fw.write('\n')

    # calc ratios for each ptm
    if i > 0:

      ptm_line = ''
      tmp_data = []

      inst_ptm = inst_line[0].replace('_','-') + '-' + str(i)
      ptm_line = ptm_line + inst_ptm+'\t'

      for j in range(len(inst_line)):
        if j in treatment_index:

          inst_data = float(inst_line[j])
          control_index = treatment_dict[str(j)]
          inst_control = float( inst_line[control_index] )

          if inst_control != 0:
            inst_ratio = inst_data/ float(inst_control)
          else:
            inst_ratio = 'nan'

          if inst_ratio != 'nan' and inst_ratio != 0:
            inst_ratio = np.log2(inst_ratio)

            # if inst_ratio < -max_ratio:
            #   inst_ratio=--3
            # elif inst_ratio > max_ratio:
            #   inst_ratio==3

          ptm_line = ptm_line + str(inst_ratio)+'\t'
          tmp_data.append(str(inst_ratio))

      ptm_line = ptm_line + '\n'

      tmp_data = list(set(tmp_data))

      if (len(tmp_data) > 6):

        # if 'nan' not in tmp_data:
        fw.write(ptm_line)

  fw.close()

def save_treatment_to_tsv():
  import numpy as np

  files = ['treated_cell_phospho_TMT_data.txt']

  for inst_file in files:

    inst_filename = 'treated_cell_12_1_2015/' + inst_file

    f = open(inst_filename,'r')
    lines = f.readlines()
    f.close()

    # 0 protein
    # 1 gene symbol
    # 2 residue
    # 3 (+/-)7 seq
    # 4 ACC#
    # 5 site group ID
    # 6 site ID

    # the number of expected data points
    num_data = 36
    # num_data = 6

    save_filename = 'treated_cell_12_1_2015/treated_cl_'+inst_file.split('_')[2]+ '.tsv'
    # inst_filename = 'lung_cellline_3_1_16/lung_cellline_phospho/lung_cellline_TMT1_phospho.txt'

    fw = open(save_filename, 'w')

    col_start = 7

    for i in range(len(lines)):

      inst_line = lines[i].strip()
      inst_line = inst_line.split('\t')

      if i == 0:
        tmp_data_names = inst_line[col_start:]

        print(tmp_data_names)

        write_line = 'ptm\t' + '\t'.join(tmp_data_names) + '\n'

        fw.write(write_line)
      if i > 0:
        inst_gene = inst_line[1]
        inst_residue = inst_line[2]
        inst_ptm = inst_gene + '_' + inst_residue

        tmp_data = inst_line[col_start:]

        # add NaNs for missing data
        inst_data = []

        for j in range(num_data):

          # fill in missing values with NaN
          try:
            data_point = tmp_data[j]

            if data_point == '':
              data_point = 'NaN'
          except:
            data_point = 'NaN'

          inst_data.append(data_point)

        # save to file
        write_line = inst_ptm + '\t'  + '\t'.join(inst_data) + '\n'

        fw.write(write_line)

    fw.close()

def calc_treatment_ratios():

  from clustergrammer import Network

  net = Network()

  net.load_tsv_to_net('treated_cell_12_1_2015/treated_cl_phospho.tsv')

def scraps():

  pass

  # # calc ratio of data
  # if inst_data[0] != 'NaN' and inst_data[0] !='0':

  #   inst_data = [float(i) for i in inst_data]

  #   inst_ctrl = inst_data[0]

  #   inst_data = [ np.log2( i/float(inst_ctrl) ) for i in inst_data]
  #   # inst_data = [ i/float(inst_ctrl) for i in inst_data]

  #   inst_data = inst_data[1:]

  #   inst_data = [str(i) for i in inst_data]

  #   inst_data = [i.replace('-inf','0') for i in inst_data]

