import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc

import matplotlib.pyplot as plt
import plotly.graph_objects as go
import plotly.express as px


def parse_panel(adata, panel, col = 'phenotype'):

    """
    Reads scimap phenotyping panel and adds a new column ('broad_phenotype')
    to anndata for the broadest classification of each cell
    (usually tumor, stroma, immune, and unknown but robust to any phenotyping scheme)
    """

    # Read in phenotype panel
    panel = pd.read_csv(panel)

    # dictionary {phenotype : parent phenotype}, keep list of broadest classes
    phenotype_hierarchy = {"unknown" : "unknown"}
    broadest_phenotypes = ["unknown"]
    for idx, pheno in enumerate(panel.iloc[:,1]):
        pheno = pheno.lower()
        parent_phenotype = panel.iloc[idx,0].lower()
        if parent_phenotype == 'all':
            broadest_phenotypes.append(pheno)
            phenotype_hierarchy[pheno] = pheno
        else:
            phenotype_hierarchy[pheno] = parent_phenotype
    
    # create a new column assigning broadest phenotypes to each cell
    pheno_list = list(adata.obs[col])
    broad_phenotype = []
    for idx, pheno in enumerate(pheno_list):
        new_pheno = pheno.lower()
        counter = 0
        while new_pheno not in set(broadest_phenotypes):
            counter += 1 
            if counter == 20:
                # an arbitrary backboard so we don't run forever.
                # this should never take 20 cycles for a given phenotype
                raise ValueError(f"Can't find hierarchical phenotype for {pheno}")
            else:
                # grab the parent phenotype
                new_pheno = phenotype_hierarchy[new_pheno]

        broad_phenotype.append(new_pheno)

    adata.obs['broad_phenotype'] = broad_phenotype

    return adata


def sample_compositional_plots(adata, cols = ["phenotype", "leiden"]):

    '''
    Plots all the compositional plots for single-sample analysis
    '''

    for i,col in enumerate(cols):

        # histograms
        fig = px.histogram(adata.obs, x = col, histnorm='percent')
        fig.write_image(f'./figures/histogram_{col}.png', scale = 5)

        # heatmaps
        ax = sc.pl.heatmap(adata, adata.var_names, groupby = col, 
        cmap='viridis', dendrogram=True, use_raw=False, swap_axes=True,
        show = False, save = f"_{col}.png")

        # matrixplots
        mp = sc.pl.matrixplot(adata, adata.var_names, groupby = col,
        use_raw=False, return_fig=True, show = False)
        mp.add_totals().style(edge_color='black')
        mp.savefig(f"./figures/matrixplot_{col}.png")

        # umap
        sc.pl.umap(adata, color = col, show = False, save = f"_{col}.png")

        # ranked marker plots on log data
        sc.tl.rank_genes_groups(adata, groupby = col, layer="log", use_raw=False)
        sc.pl.rank_genes_groups(adata, show = False, save = f"_{col}.png")

    # co-occurrence matrix between columns
    if len(cols) == 2:
        co_mat = pd.crosstab(adata.obs[cols[0]], adata.obs[cols[1]])
        co_mat = np.log1p(co_mat)
        fig = px.imshow(co_mat, labels = dict(color = "log(cell counts)"))
        fig.write_image('./figures/cooccurrence_matrix.png', scale = 5)

    for marker in adata.var_names:

        sc.pl.umap(adata, color = marker, use_raw=False, show = False, save = f"_{marker}.png")

    return


def calculate_sample_stats(metadata, multi_adata, col = 'phenotype'):

    """
    Populates metadata table with the following metrics for each sample:
        - total cell count (column = count)
        - tissue area in microns^2 (column = area)
        - counts of cell phenotypes (column = phenotype_count, for ea. phenotype)*
        - densities of cell phenotypes (column = phenotype_density, for ea. phenotype)*
        - marker positivity ratios (column = marker_ratio, for ea. marker in X)

    * By default, this is broad_phenotype, but can use any categorical column specified as col arg
    """

    # initialize new cols
    new_cols = {'count' : [], 'area' : []}
    
    # populate area and counts columns
    for i, sample in enumerate(metadata['sample']):

        # subset obs to sample
        sample_obs = multi_adata.obs[multi_adata.obs['sample'] == sample]

        # cell counts
        new_cols['count'].append(len(sample_obs))

        # convert res to microns^2/px^2, sum cell areas in sample, convert to microns^2
        sq_resolution = metadata['resolution'][i]**2
        area = sample_obs['Area'].sum() * sq_resolution
        new_cols['area'].append(round(area, 2))

        # phenotype (default col arg) counts and densities per sample
        value_counts_table = sample_obs[col].value_counts()
        for phenotype in value_counts_table.index:

            phenotype_count = value_counts_table[phenotype]
            phenotype_density = phenotype_count / new_cols['area'][i]

            if f"{phenotype}_count" not in new_cols:
                new_cols[f"{phenotype}_count"] = [phenotype_count]
                new_cols[f"{phenotype}_density"] = [phenotype_density]
            else:
                new_cols[f"{phenotype}_count"].append(phenotype_count)
                new_cols[f"{phenotype}_density"].append(phenotype_density)

        # marker positivity ratios per sample based on binary layer
        sample_binary = multi_adata.layers['binary'][multi_adata.obs['sample'] == sample]
        sample_binary = pd.DataFrame(sample_binary, columns = multi_adata.var_names)

        for marker in sample_binary.columns:

            # calculate ratio of positive cells to negative cells
            colsum = sample_binary[marker].sum()
            ratio = round((colsum / (len(sample_binary) - colsum)),2)

            if f"{marker}_ratio" not in new_cols:
                new_cols[f"{marker}_ratio"] = [ratio]
            else:
                new_cols[f"{marker}_ratio"].append(ratio)
            
    # populate metadata table with new metrics
    for col in new_cols:
        metadata[col] = new_cols[col]

    return metadata