import re
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc

import matplotlib.pyplot as plt
import plotly.graph_objects as go
import plotly.express as px

from sklearn.neighbors import KDTree


def count_neighbors_by_type(adata, radius, col = "phenotype"):

    ''' 
    Create a list of dictionaries that contains neighbor phenotype counts for all cells
    '''

    print(f"Finding {col} neighbors")

    phenotypes = adata.obs[col]

    # Create tree for queries.
    cell_coords = adata.obs[['X_centroid', 'Y_centroid']]
    kdt = KDTree(cell_coords)

    # Get neighbors.
    neighbors = kdt.query_radius(cell_coords, r=radius)

    # Remove the query cell from the neighbors list, which is the index of the neighbors list.
    neighbors = [n[n != i] for i, n in enumerate(neighbors)]

    # Count phenotypes in neighbors.
    def count_phenotypes(n):
        p, c = np.unique(phenotypes[n], return_counts=True)
        return dict(zip(p, c))

    neighbor_matrix = np.array(list(map(count_phenotypes, neighbors)))
    
    # Create a temp dataframe with all neighbors data and total neighbors
    neighbor_df = pd.DataFrame(list(neighbor_matrix)).fillna(0)
    neighbor_df.columns = ["%s_neighbors" % n for n in neighbor_df.columns]
    neighbor_df['total_neighbors'] = neighbor_df.sum(axis=1)
    
    # merge temp df w/ adata.obs
    for i,col in enumerate(neighbor_df.columns):
        adata.obs[col] = list(neighbor_df[col])

    return adata


def get_region(neighbors):

    """ 
    Assigns region for a single cell based on the proportion of surrounding cell's phenotypes
    """
    region = ""
    count_neighbor_types = neighbors[neighbors > 0].count()
    if count_neighbor_types > 1:
        # Potential border region: make it a border region if dominant cell type is not overwhelming % of neighbors.
        # To avoid any border regions and classify strictly into tumor/immune/stroma, set condition to
        # if dominant_type_percent < 0:
        dominant_type_percent = neighbors.max() / neighbors.sum()
        if dominant_type_percent < 0.75:
            region = '_'.join(list(neighbors[neighbors > 0].index)[0:2]) + " border"
        else:
            region = neighbors.idxmax()     
    elif count_neighbor_types == 1:
        # Homogenous region.
        region = neighbors.idxmax()
    else:
        # No neighbors.
        region = "Unassigned"
        
    return region


def assign_all_regions(adata, col = "phenotype"):

    """
    Apply get_region() to every cell in dataframe
    """
    print("Assigning regions")
    # make a temp df that filters down to cols relevant for region assignment
    # and removes the _neighbors labels from the column names
    cols_filt = list(adata.obs[col].unique())
    cols_filt = [c + "_neighbors" for c in cols_filt]
    if 'unknown_neighbors' in cols_filt:
        cols_filt.remove('unknown_neighbors')

    df = adata.obs[cols_filt]
    df.columns = [c[0:-10] for c in df.columns]

    # Compute cell regions.
    region = df.apply(get_region, axis=1)
    adata.obs = adata.obs.assign(region=region.values)

    # remove redundancy from borders ("stroma_tumor border" == "tumor_stroma border")
    unique_regions = adata.obs['region'].unique()
    renaming_dict = {}
    for label in unique_regions:
        if "border" in label:
            sorted_label = sorted(re.split("_", label.removesuffix(" border")))
            rebuilt_label = sorted_label[0] + "_" + sorted_label[1] + " border"
            renaming_dict[label] = rebuilt_label
        else:
            renaming_dict[label] = label

    adata.obs['region'] = [renaming_dict[x] for x in list(adata.obs['region'])]

    return adata


def calculate_positivity_ratio(marker_column):

    """
    Calculate the positivity ratio of a marker from a single column of anndata binary layer
    """

    # calculate ratio of positive cells to negative cells
    colsum = marker_column.sum()
    ratio = round((colsum / (len(marker_column) - colsum)),2)

    return ratio
    

def get_spatial_features(adata, resolution, col = 'phenotype'):

    """
    Make a new df w/ regions as rows and fill in regional features as columns
    Features generated: area, total_cells, imageid, {phenotype}_count, {phenotype}_proportion, {phenotype}_density
    """
    print("Calculating spatial features")
    # convert resolution to (microns/px)^2
    resolution = float(resolution)**2

    # make region df initialized with total cell counts
    regions_df = pd.DataFrame(adata.obs['region'].value_counts())
    regions_df = regions_df.reset_index()
    regions_df.columns = ['region', 'total_cells']

    # double check this is a single sample, and add the imageid to region df
    if adata.obs['imageid'].nunique() == 1:
        regions_df['imageid'] = list(adata.obs['imageid'])[0:len(regions_df.index)]

    # get region areas and phenotype counts
    new_cols = {'area' : []}
    for region in regions_df['region']:

        df_temp = adata.obs[adata.obs['region'] == region]

        # region areas, convert from px^2 to microns^2
        area = df_temp['Area'].sum()
        area = area * resolution
        new_cols['area'].append(area)

        # phenotype counts
        pheno_counts = df_temp[col].value_counts()
        for phenotype in adata.obs[col].unique():

            colname = f"{phenotype}_count"

            # handle missing phenotypes 
            if phenotype not in pheno_counts.index:
                phenotype_count = 0
            else:
                phenotype_count = pheno_counts[phenotype]
            
            # generate list of phenotype counts
            if colname not in new_cols:
                new_cols[colname] = [phenotype_count]
            else:
                new_cols[colname].append(phenotype_count)

    # populate regions dataframe with areas and phenotype counts
    for key in new_cols:
        regions_df[key] = new_cols[key]

    # calculate density and percent of each phenotype per region
    count_columns = [x for x in regions_df.columns if re.search("_count", x)]

    for pheno in count_columns:
        prefix = re.split("_", pheno)[0]
        # proportion
        regions_df[f"{prefix}_proportion"] = (regions_df[pheno] / regions_df['total_cells'])
        # density
        regions_df[f"{prefix}_density"] = (regions_df[pheno] / regions_df['area'])

    # save the spatial feature table in anndata
    adata.uns['spatial_feature_table'] = {
        'spatial_vars' : np.array(regions_df.columns),
        'spatial_regions' : np.array(regions_df['region']),
        'spatial_x' : regions_df
    }

    return adata


def plot_spatial(adata, columns = ['phenotype', 'region']): # need cols as an argument

    """
    Generates stacked barplots for phenotypes w/in each region (by count, proportion, and density)
    and spatial plots for phenotypes, leiden clusters, regions, and re-scaled marker intensities
    """

    regions_df = adata.uns['spatial_feature_table']['spatial_x']

    # stacked barplots
    for metric in ["count", "proportion", "density"]:

        # subsetting to metric and making the names pretty for the legend
        subset = [x for x in regions_df.columns if re.search(metric, x)]
        temp_df = regions_df.loc[:,subset]
        temp_df.columns = [re.split("_", x)[0] for x in temp_df.columns]
        temp_df['region'] = regions_df['region']

        fig = px.bar(temp_df, x= 'region', y= temp_df.columns, 
        labels = {'value': metric, 'variable':'phenotype'})
        fig.write_image(f"./figures/regional_stacked_{metric}.png", scale = 5)

    # Add Spatial Coordinates
    adata.obsm['spatial'] = np.stack([adata.obs['X_centroid'].values,adata.obs['Y_centroid'].values], axis=1)

    # spatial plots for various classifications
    # removed "leiden" from cols for moment
    for col in columns:
        sc.pl.spatial(adata, spot_size=20, color = col, save = f"_{col}.png", show = False)

    # spatial plots for re-scaled intensities of each marker
    for marker in adata.var_names:
        sc.pl.spatial(adata, spot_size=20, color = marker, use_raw=False, save = f"_{marker}.png", show = False)

    # sc.pl.spatial(adata, spot_size=20, color = adata.var_names, use_raw=False, save = f"_markers.png", show = False)

    return