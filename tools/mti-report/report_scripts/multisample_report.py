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
import compositional
import spatial

def get_args():

    parser = argparse.ArgumentParser(add_help = True,
    description = "Multisample report generation -- assumes single sample analysis has been complete for all inputs")

    parser.add_argument("-i", "--input", 
    help = "Path to directory containing h5ad files for multiple samples",
    required = True)

    parser.add_argument("-m", "--metadata", 
    help = "CSV with columns 'file','sample','resolution' (w/ resolution in microns per px)",
    required = True)

    # parser.add_argument("-c", "--column", 
    # help = "column name for phenotypes, default is 'phenotype",
    # required = False,
    # default = 'phenotype')

    parser.add_argument("--resolution", 
    help = "Image resolution in microns/px for density calculation (default is 0.65 micron/px)",
    required = False,
    default = 0.65)

    return parser.parse_args()


def main():

    args = get_args()
    
    # concatenate sample h5ad files based on metadata table
    metadata = pd.read_csv(args.metadata)
    multi_adata, multi_spatial = preprocess.concat_h5ads(metadata, args.input)

    # marker positivity calls
    multi_adata = preprocess.create_binary_layer(multi_adata)

    # populate metadata table with sample metrics
    metadata = compositional.calculate_sample_stats(metadata, multi_adata, col = 'broad_phenotype')

    

    return 0


if __name__ == '__main__':
    sys.exit(main())