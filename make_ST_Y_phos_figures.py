
from clustergrammer import Network
net = Network()

net.load_file('txt/lung_CL_phos_ratios_ST.tsv')

# quantile normalize to normalize cell lines
# net.normalize(axis='col', norm_type='qn')

# only keep most differentially expressed genes
net.filter_N_top('row', 500, 'sum')

# take zscore
# net.normalize(axis='row', norm_type='zscore', keep_orig=True)

net.swap_nan_for_zero()

# # original
# net.filter_threshold('row', 1.75, 4)
# # modified version
# net.filter_threshold('row', 1.75, 4)


views = ['N_row_sum','N_row_var']

net.make_clust(dist_type='cos',views=views, dendro=True,
               sim_mat=True)


net.write_json_to_file('viz', 'json/mult_view.json')

