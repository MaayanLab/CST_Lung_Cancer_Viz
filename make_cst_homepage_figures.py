def main():

  make_phos_homepage_viz()

  make_exp_homepage_viz()

def make_phos_homepage_viz():

  from clustergrammer import Network
  net = Network()

  net.load_file('lung_cellline_3_1_16/lung_cellline_phospho/lung_cellline_TMT_phospho_combined_ratios.tsv')

  # quantile normalize to normalize cell lines
  net.normalize(axis='col', norm_type='qn')

  # only keep most differentially regulated PTMs
  net.filter_N_top('row', 250, 'sum')

  # take zscore of rows
  net.normalize(axis='row', norm_type='zscore', keep_orig=True)

  net.swap_nan_for_zero()

  # threshold filter PTMs
  net.filter_threshold('row', threshold=1.75, num_occur=3)

  views = ['N_row_sum','N_row_var']
  net.make_clust(dist_type='cos',views=views, dendro=True,
                 sim_mat=True, calc_cat_pval=False)

  net.write_json_to_file('viz', 'json/homepage_phos.json', 'indent')

def make_exp_homepage_viz():

  from clustergrammer import Network
  net = Network()

  net.load_file('CCLE_gene_expression/CCLE_NSCLC_all_genes.txt')

  # threshold filter expression
  net.filter_threshold('row', threshold=3.0, num_occur=4)

  views = ['N_row_sum', 'N_row_var']
  net.make_clust(dist_type='cos',views=views, dendro=True,
                 sim_mat=True, calc_cat_pval=False)

  net.write_json_to_file('viz', 'json/homepage_exp.json', 'indent')

main()