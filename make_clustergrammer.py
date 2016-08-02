import time
start_time = time.time()

from clustergrammer import Network
net = Network()

# choose tsv file
####################
inst_name = 'ST'
net.load_file('txt/phos_ratios_all_treat_no_geld_ST.txt')
# net.load_file('txt/phos_ratios_all_treat_no_geld_Tyrosine.txt')


net.swap_nan_for_zero()

# net.normalize(axis='row', norm_type='zscore', keep_orig=True)

print(net.dat.keys())

views = ['N_row_sum', 'N_row_var']

net.make_clust(dist_type='cos',views=views , dendro=True,
               sim_mat=True, filter_sim=0.1, calc_cat_pval=False)
               # run_enrichr=['KEA_2015'])
               # run_enrichr=['ENCODE_TF_ChIP-seq_2014'])
               # run_enrichr=['GO_Biological_Process_2015'])

net.write_json_to_file('viz', 'json/'+inst_name+'.json', 'no-indent')
net.write_json_to_file('sim_row', 'json/'+inst_name+'_sim_row.json', 'no-indent')
net.write_json_to_file('sim_col', 'json/'+inst_name+'_sim_col.json', 'no-indent')

elapsed_time = time.time() - start_time
print('\n\nelapsed time: '+str(elapsed_time))
