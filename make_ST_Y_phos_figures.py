
from clustergrammer import Network
from copy import deepcopy

all_filenames = ['lung_CL_phos_ratios_ST', 'lung_CL_phos_ratios_Y']


for inst_filename in all_filenames:

  if inst_filename == 'lung_CL_phos_ratios_ST':
    num_reach = 4
  else:
    num_reach = 2

  net = deepcopy(Network())

  load_filename = 'txt/' + inst_filename + '.tsv'
  net.load_file(load_filename)

  # quantile normalize to normalize cell lines
  net.normalize(axis='col', norm_type='qn')

  # only keep most differentially expressed genes
  net.filter_N_top('row', 250, 'sum')

  # take zscore
  net.normalize(axis='row', norm_type='zscore', keep_orig=True)

  net.swap_nan_for_zero()

  # # original
  # net.filter_threshold('row', 1.75, 4)
  # # modified version
  net.filter_threshold('row', 1.75, num_reach)

  # views = ['N_row_sum','N_row_var']
  views = []

  net.make_clust(dist_type='cos',views=views, dendro=True,
                 sim_mat=True)

  json_filename = 'json/' + inst_filename + '.json'
  net.write_json_to_file('viz', json_filename)
