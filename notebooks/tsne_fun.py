
from copy import deepcopy
from clustergrammer import Network

net = deepcopy(Network())
cl_info = net.load_json_to_dict('../cell_line_info/cell_line_info_dict.json')

def make_simple_cl_names(tuple_cols):
    cl_names = []

    if type(tuple_cols[0]) is tuple:
        for inst_tuple in tuple_cols:
            cl_names.append(inst_tuple[0])
    else:
        cl_names = tuple_cols

    return cl_names

def make_cmap(cl_names, cmap_type):
    inst_cmap = []
    color_dict = make_color_dict()

    if cmap_type == 'hist':
        for inst_cl in cl_names:
            inst_hist = cl_info[inst_cl]['Histology']
            if inst_hist == 'NSCLC':
                inst_cmap.append('red')
            else:
                inst_cmap.append('blue')
        return inst_cmap

    elif cmap_type == 'plex':
        for inst_cl in cl_names:
            inst_plex = cl_info[inst_cl]['Plex']

            inst_plex = int(inst_plex)

            inst_color = color_dict[inst_plex]

            inst_cmap.append(inst_color)

    elif cmap_type == 'exp-group':
        for inst_cl in cl_names:
            inst_plex = cl_info[inst_cl]['Exp-group']

            inst_plex = int(inst_plex)

            inst_color = color_dict[inst_plex]

            inst_cmap.append(inst_color)

    return inst_cmap

def make_color_dict():
    color_dict = {-1:'white', 1:'red', 2:'blue', 3:'brown',
    4:'lime', 5:'yellow', 6:'navy', 7:'gold', 8:'cyan', 9:'pink'}
    return color_dict

def make_multiple_cl_tsne(mat, cmap_left=None, cmap_right=None,
                           skl_version=True, random_state=0,
                           learning_rate=40):
    from matplotlib import pyplot as plt
    import numpy as np

    # the matrix needs to be transposed in order to cluster the numbers
    x_data = mat.transpose()

    # convert image data to float64 matrix. float64 is need for bh_sne
    x_data = np.asarray(x_data).astype('float64')

    if skl_version == False:
        from tsne import bh_sne
        # perform t-SNE embedding, lowered perplexity
        vis_data = bh_sne(x_data, perplexity=7)
        vis_x = vis_data[:, 0]
        vis_y = vis_data[:, 1]

    else:
        from sklearn import manifold
        # run tsne from sklearn
        ###########################
        tsne = manifold.TSNE(perplexity=7, n_iter=100000,
            random_state = random_state, method='exact', metric='correlation',
            learning_rate=learning_rate, verbose=0,
            n_iter_without_progress=1000, init='random', early_exaggeration=4)

        Y = tsne.fit_transform(x_data)
        vis_x = Y[:, 0]
        vis_y = Y[:, 1]

    fig, axarr = plt.subplots(ncols=2, figsize=(10,5))

    marker_size = 150

    # always require cmap
    axarr[0].scatter(vis_x, vis_y, c=cmap_left, \
        cmap=plt.cm.get_cmap('prism',len(cmap_left)), s=marker_size)

    axarr[1].scatter(vis_x, vis_y, c=cmap_right, \
        cmap=plt.cm.get_cmap('jet',len(cmap_right)), s=marker_size)

    plt.show()

def normalize_and_make_tsne(data_type= 'phospho', qn_col=False, zscore_row=False,
                            filter_missing=False, skl_version=True,
                            random_state=0, learning_rate=40):

    paths = {}
    paths['phospho'] = '../lung_cellline_3_1_16/lung_cellline_phospho/' + \
        'lung_cellline_TMT_phospho_combined_ratios.tsv'
    paths['exp'] = '../CCLE_gene_expression/CCLE_NSCLC_all_genes.txt'

    filename = paths[data_type]

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

    print('matrix shape '+str(mat.shape))

    cmap = get_cmaps(inst_df)

    make_multiple_cl_tsne(mat, cmap_left=cmap['hist'], cmap_right=cmap['plex'],
                           skl_version=skl_version, random_state=random_state,
                           learning_rate=learning_rate)

def get_cmaps(inst_df):

    # make colormaps based on histology and plexes
    tuple_cols = inst_df['mat'].columns.tolist()

    cl_names  = make_simple_cl_names(tuple_cols)

    cmap = {}
    cmap['hist'] = make_cmap(cl_names, 'hist')
    cmap['plex'] = make_cmap(cl_names, 'plex')
    cmap['exp-group'] = make_cmap(cl_names, 'exp-group')

    return cmap
