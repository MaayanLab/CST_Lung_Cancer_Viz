import time
# import StringIO

start_time = time.time()

# import network class from Network.py

from clustergrammer import Network
net = Network()

print('here')
net.load_file('lung_cellline_3_1_16/lung_cellline_phospho/lung_cellline_TMT_phospho_combined_ratios.tsv')


# quantile normalize to normalize cell lines
net.normalize(axis='col', norm_type='qn')

# only keep most differentially expressed genes
net.filter_N_top('row', 250, 'sum')

# take zscore
net.normalize(axis='row', norm_type='zscore', keep_orig=True)

net.swap_nan_for_zero()

# original
# net.filter_threshold('row', 1.75, 4)
# modified version
net.filter_threshold('row', 1.75, 6)


views = ['N_row_sum','N_row_var']

# net.swap_nan_for_zero()

net.make_clust(dist_type='cos',views=views, dendro=True,
               sim_mat=True)


net.produce_view({'N_row_sum':10,'dist':'euclidean'})

net.write_json_to_file('viz', 'json/mult_view.json', 'indent')

# your code
elapsed_time = time.time() - start_time

print('\n\nelapsed time')
print(elapsed_time)
