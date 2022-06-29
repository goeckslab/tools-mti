import os
import re
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import matplotlib.pyplot as plt
import plotly.graph_objects as go


def remove_unwanted_markers(adata, patterns = ['DAPI', 'AF']):

    """
    Removes unwanted markes from main dataframe (adata.X).
    adata.raw.X is unchanged, and the raw markers get saved in adata.uns['raw_var_names']
    """

    print(f"Removing markers matching {patterns}")

    # for every marker, see if it matches each pattern provided and keep track
    markers_to_remove = []
    for pattern in patterns:
        for marker in adata.var_names:

            if re.search(pattern, marker):

                markers_to_remove.append(marker)

    # subset adata.var_names by markers not identified above
    markers_to_keep = [x for x in adata.var_names if x not in markers_to_remove]

    # filter X by filtered var_names
    adata = adata[:,markers_to_keep]
    print("Remaining markers after removal: ", adata.var_names)

    return adata


def create_log_layer(adata):

    """
    Create a layer in anndata called 'log' that is np.log1p(adata.raw.X).
    Shape of adata.X is used (so any markers in raw but no longer in X will not be in log)
    """

    # adata.raw.X is immutable so make a temp df
    raw_df = pd.DataFrame(adata.raw.X, columns=adata.raw.var_names)
    
    # subset by markers in adata.X
    raw_df = raw_df.loc[:,adata.var_names]

    # create log layer
    adata.layers["log"] = raw_df
    adata.layers["log"] = np.log1p(adata.layers["log"])

    # also put the logged dataframe in obsm for sc.pp.neighbors
    adata.obsm['X_log'] = np.array(raw_df)

    return adata


def create_binary_layer(adata):

    """
    Creates a new layer in adata with positive (1) or
    negative (0) for every marker in every cell
    """

    # make sure X has been re-scaled between 0 and 1
    if np.all((adata.X >= 0) & (adata.X <= 1)):

        print("X normalized between 0 and 1")
        
        adata.layers['binary'] = np.where(adata.X > 0.5, 1, 0)

    else:

        raise ValueError("X data must be rescaled between 0-1 to call marker positivity")

    return adata


def calculate_distances(adata, seed, multisample = False):

    '''
    Calculates distance matrices and embeddings necessary for analysis.
    Returns nothing, each metric is automatically saved in the anndata object
    '''

    if multisample:
        print('Running neighbor graph, PCA, UMAP in multisample mode')
        rep = 'X'
        layer = adata.X
    else:
        print('Running neighbor graph, PCA, UMAP in single-sample mode')
        rep = 'X_log'
        layer = adata.layers["log"]
    
    # This defines neighbors based on biomarker expression (not to be mistaken w/ spatial neighbors)
    sc.pp.neighbors(adata, random_state = seed, use_rep = rep)

    # PCA dist matrix
    X_pca = sc.tl.pca(layer, random_state = seed, return_info = False)
    adata.obsm['X_pca'] = X_pca

    # UMAP 
    sc.tl.umap(adata, random_state = seed)

    print("The followings keys have been added to adata.obsm: ", adata.obsm)

    return adata


def concat_h5ads(metadata, input_dir):

    """
    Concatenates multiple anndata h5ad file using sample name from metadata as index,
    also concatenate all the spatial feature tables together with sample name
    """

    # create dict: {'sample_1' : AnndataObject_1,...'sample_n' : AnndataObject_n}
    adatas = {}
    spatial_tables = []
    for i in metadata.index:

        # read sample anndata, add it to dictionary
        sample_path = os.path.join(input_dir, metadata['file'][i])
        sample_anndata = ad.read_h5ad(sample_path)
        adatas[metadata['sample'][i]] = sample_anndata

        # extract spatial feature table from sample anndata, add to list of spatial tables
        sample_spatial_table = pd.DataFrame(sample_anndata.uns['spatial_feature_table']['spatial_x'])
        sample_spatial_table['sample'] = metadata['sample'][i]
        spatial_tables.append(sample_spatial_table)
        

    # concatenate into single multi-sample anndata object w/ sample name tracking
    multi_adata = ad.concat(adatas, label = "sample")

    # concatenate spatial feature tables
    multi_spatial = pd.concat(spatial_tables, ignore_index = True)

    return multi_adata, multi_spatial