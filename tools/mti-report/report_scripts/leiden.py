import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import matplotlib.pyplot as plt

from sklearn.metrics import silhouette_score, silhouette_samples


def optimize_clustering(adata, resolutions, seed):

    '''
    Runs leiden clustering on adata.X across a range of resolutions.
    Returns dataframe of num clusters and silhouette score for each resolution tested.
    Also returns dict of the best scoring resolution 
    '''
    clusters = []
    sils = []
    for res in resolutions:
        res = float(res)
        # leiden clustering at each resolution
        sc.tl.leiden(adata, resolution = res, key_added = f'leiden_{res}', random_state = seed)
        # assess clustering quality
        n_clusts = adata.obs[f'leiden_{res}'].nunique()
        sil = silhouette_score(adata.obsm['X_pca'], adata.obs[f'leiden_{res}'])
        print(f'resolution of {res} produced {n_clusts} clusters with a silhoutte score of {sil}')
        clusters.append(n_clusts)
        sils.append(sil)

    # create a dataframe containing each resolution, number of clusters generated, and the avg. sil score
    sil_df = pd.DataFrame({'resolution' : resolutions,
                           'clusters' : clusters,
                           'score' : sils})

    # find the best avg. sil score, it's corresponding resolution, and number of clusters
    best = {'score' : sil_df['score'].max(),
            'clusts' : sil_df['clusters'][sil_df['score'].idxmax()],
            'resolution' : sil_df['resolution'][sil_df['score'].idxmax()]
            }

    # clean-up adata.obs: only retain the best clustering in the anndata file and update the leiden provenance
    for i,res in enumerate(sil_df['resolution']):
        res = float(res)
        if res == float(best['resolution']):
            adata.obs['leiden'] = adata.obs[f'leiden_{res}']
            adata.uns['leiden']['params']['resolution'] = res
            adata.obs = adata.obs.drop(columns = f'leiden_{res}')
        else:
            adata.obs = adata.obs.drop(columns = f'leiden_{res}')

    # see if a broad enough range of resolutions were tested
    # fig, ax1 = plt.subplots()
    # ax2 = ax1.twinx()
    # ax1.plot(sil_df['resolution'],sil_df['score'] , 'g-')
    # ax2.plot(sil_df['resolution'],sil_df['clusters'], 'b-')
    # ax1.grid(False)
    # ax2.grid(False)
    # ax1.set_xlabel('Resolution')
    # ax1.set_title('Leiden clustering optimization')
    # ax1.set_ylabel('Silhouette score', color='g')
    # ax2.set_ylabel('Number of clusters', color='b')
    # plt.tight_layout()
    # fig.savefig("./figures/sil_by_res.png", dpi = 300)

    return adata


def plot_sample_scores(adata, col):

    '''
    Calculates cell-specific silhouette scores for a clustering or phenotype (or any categorical var),
    saves to anndata, and generates silhouette plot
    '''
    # get list of unique categories (clusters, phenotypes...) in size order (how many cells belong to each)
    uniques = list(adata.obs[col].value_counts().index)

    # calculate the silhoutte score for each cell in the best clustering
    sample_sils = silhouette_samples(adata.obsm['X_pca'], adata.obs[col])
    adata.obs[f'silhouette_{col}'] = sample_sils

    # plot the grouped cell silhouette scores
    fig = plt.figure()
    ax = fig.add_axes([0.1,0.1,0.75,0.75])
    fig.set_size_inches(18, 7)

    # sil scores range between [-1,1], but we'll constrain the plot to min/max of data + padding
    ax.set_xlim([sample_sils.min() - 0.001, sample_sils.max() + 0.001])
    ax.set_ylim([0, len(adata.obs) + (len(uniques) + 1) * 10])

    y_lower = y_upper = 0
    for i,n in enumerate(uniques):
        # Aggregate the scores for samples belonging to cluster i, and sort them
        ith_cluster_silhouette_values = sample_sils[adata.obs[col] == n]
        ith_cluster_silhouette_values.sort()

        # set the upper bound of the cluster bar
        y_upper += len(ith_cluster_silhouette_values)

        # bar plot for each cluster, bounded along the y-axis by num cells in cluster
        ax.barh(range(y_lower,y_upper),ith_cluster_silhouette_values,height =1, label = n)

        # set the next cluster's lower bound as the previous cluster's upper bound
        y_lower += len(ith_cluster_silhouette_values)

    # formatting
    ax.set_title(f"Cell silhoutte scores grouped by {col}")
    ax.set_xlabel("Silhouette score")
    ax.set_yticks([])  

    # legend inherits labels from barh(label)
    ax.legend(title = col)

    fig.savefig(f"./figures/{col}_sils.png", dpi = 300)

    return