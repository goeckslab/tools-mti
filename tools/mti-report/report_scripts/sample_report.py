#!/usr/bin/env python

import sys
import os
import re
import argparse

import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc

import matplotlib.pyplot as plt
import plotly.graph_objects as go
import plotly.express as px

from sklearn.neighbors import KDTree
from sklearn.metrics import silhouette_score, silhouette_samples

import preprocess
import leiden
import compositional
import spatial

def get_args():

    # default resolution range
    default_leiden = list(np.linspace(.2, 1.5, 14))
    default_leiden = [round(r, 2) for r in default_leiden]

    parser = argparse.ArgumentParser(add_help = True,
    description = "")

    parser.add_argument("-f", "--file", 
    help = "anndata h5ad cell feature table",
    required = True)

    parser.add_argument("-p", "--panel", 
    help = "Scimap phenotyping panel in CSV format",
    required = True)

    parser.add_argument("-c", "--column", 
    help = "column name for phenotypes, default is 'phenotype",
    required = False,
    default = 'phenotype')

    parser.add_argument("--removeMarkers", 
    nargs = "+",
    help = "Patterns to remove markers by. Example: 'DAPI' will remove 'DAPI_1','DAPI_2', etc. \
            Markers are removed from adata.X but maintained in raw. default is 'DAPI' 'AF'",
    required = False,
    default = ['DAPI', 'AF'])

    parser.add_argument("--nuclearMarker", 
    help = "Name of nuclear marker in every cycle for QC purposes, default is 'DAPI'",
    required = False,
    default = 'DAPI')

    parser.add_argument("--resolution", 
    help = "Image resolution in microns/px for density calculation (default is 0.65 micron/px)",
    required = False,
    default = 0.65)

    parser.add_argument("--leidenRange", 
    nargs = "+",
    help = "list of resolutions to test (default is between (0.2,1.5) with 0.1 step)", 
    required = False,
    default = default_leiden)

    parser.add_argument("--radius", 
    help = "Radius (pixels) for neighborhood search (default is 30px)",
    required = False,
    default = 30)

    return parser.parse_args()


def main():

    SEED = 20

    # get user args and make output dir for figures
    args = get_args()
    os.mkdir(os.path.join(os.getcwd(),"figures"))

    # scanpy figure settings
    sc.set_figure_params(dpi_save = 500, transparent = True, figsize = (10,10))

    # Read anndata file
    print("Reading anndata file")
    adata = ad.read_h5ad(args.file)

    # cycle QC figures
    preprocess.plot_cycle_dropout(adata, pattern = args.nuclearMarker)

    # remove unwanted markers and make log and positivity layers
    adata = preprocess.remove_unwanted_markers(adata, patterns = args.removeMarkers)
    adata = preprocess.create_log_layer(adata)
    adata = preprocess.create_binary_layer(adata)

    # PCA, umap, expression neighborhood
    print("Calculating PCA, UMAP, and neighborhood graph")
    adata = preprocess.calculate_distances(adata, seed = SEED)

    # find optimal leiden clustering in range of resolutions, plot silhouette scores
    print(f"Running leiden clustering for each resolution in range {args.leidenRange}")
    adata = leiden.optimize_clustering(adata, resolutions = args.leidenRange, seed = SEED)
    leiden.plot_sample_scores(adata, "leiden")
    leiden.plot_sample_scores(adata, args.column)

    # parse the scimap phenotyping panel to get broader classifications
    adata = compositional.parse_panel(adata, panel = args.panel, col = args.column)

    # generate compositional plots 
    compositional.plot_histograms(adata)
    compositional.compositional_plots(adata)

    # create a spatial neighborhood graph and assign regions
    adata = spatial.count_neighbors_by_type(adata, args.radius)
    adata = spatial.assign_all_regions(adata, col = 'broad_phenotype')

    # analyze and plot regions
    adata = spatial.get_spatial_features(adata, args.resolution, col = 'broad_phenotype')
    spatial.plot_spatial(adata)

    # output h5ad naming
    bname = os.path.basename(args.file)
    bname = re.split('\.', bname)[0]
    out = bname + "_analyzed.h5ad"
    out_path = os.path.join(os.getcwd(),out)
    adata.write_h5ad(out_path)

    return 0


if __name__ == "__main__":
    sys.exit(main())