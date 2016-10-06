
from copy import deepcopy
from clustergrammer import Network

net = deepcopy(Network())
cl_info = net.load_json_to_dict('../cell_line_info/cell_line_info_dict.json')

def make_simple_cl_names(tuple_cols):
    cl_names = []
    for inst_tuple in tuple_cols:
        cl_names.append(inst_tuple[0])
    return cl_names


def make_hist_cmap(cl_names):
    hist_cmap = []
    for inst_cl in cl_names:
        inst_hist = cl_info[inst_cl]['Histology']
        if inst_hist == 'NSCLC':
            hist_cmap.append(0)
        else:
            hist_cmap.append(1)
    return hist_cmap

def make_plex_cmap(cl_names):
    plex_cmap = []
    for inst_cl in cl_names:
        inst_plex = cl_info[inst_cl]['Plex']
        plex_cmap.append(inst_plex)
    return plex_cmap

def make_cl_tsne_hist_plex(mat, cmap_left=None, cmap_right=None):
    from matplotlib import pyplot as plt
    from tsne import bh_sne
    import numpy as np
    # the matrix needs to be transposed in order to cluster the numbers
    x_data = mat.transpose()

    # convert image data to float64 matrix. float64 is need for bh_sne
    x_data = np.asarray(x_data).astype('float64')

    # perform t-SNE embedding, lowered perplexity
    vis_data = bh_sne(x_data, perplexity=7)

    # plot the result
    vis_x = vis_data[:, 0]
    vis_y = vis_data[:, 1]

    fig, axarr = plt.subplots(ncols=2, figsize=(10,5))

    marker_size = 150

    if cmap_left == None:
        axarr[0].scatter(vis_x, vis_y, s=marker_size)
    else:
        axarr[0].scatter(vis_x, vis_y, c=cmap_left, cmap=plt.cm.get_cmap('prism',len(cmap_left)), s=marker_size)

    if cmap_right == None:
        axarr[1].scatter(vis_x, vis_y, marker_size)
    else:
        axarr[1].scatter(vis_x, vis_y, c=cmap_right, cmap=plt.cm.get_cmap('jet',len(cmap_right)), s=marker_size)

    plt.show()

def normalize_and_make_tsne(qn_col=False, zscore_row=False,
                            filter_missing=False):

    filename = '../lung_cellline_3_1_16/lung_cellline_phospho/' + \
    'lung_cellline_TMT_phospho_combined_ratios.tsv'

    # get matrix of data for tsne
    net = deepcopy(Network())
    net.load_file(filename)

    if qn_col == True:
        print('quantile normalize columns')
        net.normalize(axis='col', norm_type='qn')

    if zscore_row == True:
        print('zscore rows')
        net.normalize(axis='row', norm_type='zscore')

    if filter_missing == True:
        print('filter PTMs with missing values')
        # only include PTMS with no missing data
        net.filter_threshold('row', 0.0, 45)

    net.swap_nan_for_zero()
    inst_df = net.dat_to_df()
    mat = inst_df['mat'].values

    print(mat.shape)

    hist_cmap, plex_cmap = get_cmaps(inst_df['mat'])

    make_cl_tsne_hist_plex(mat, cmap_left=hist_cmap, cmap_right=plex_cmap)

def get_cmaps(inst_df):
    filename = '../lung_cellline_3_1_16/lung_cellline_phospho/' + \
        'lung_cellline_TMT_phospho_combined_ratios.tsv'

    # get matrix of data for tsne
    net = deepcopy(Network())
    net.load_file(filename)
    net.swap_nan_for_zero()
    inst_df = net.dat_to_df()

    # make colormaps based on histology and plexes
    tuple_cols = inst_df['mat'].columns.tolist()

    cl_names  = make_simple_cl_names(tuple_cols)
    hist_cmap = make_hist_cmap(cl_names)
    plex_cmap = make_plex_cmap(cl_names)

    return hist_cmap, plex_cmap
