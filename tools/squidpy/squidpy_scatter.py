import argparse
import ast
import json
import warnings

import pandas as pd
import squidpy as sq
from anndata import read_h5ad


def main(inputs, anndata, output, output_plot):
    """
    Parameter
    ---------
    inputs : str
        File path to galaxy tool parameter.
    anndata : str
        File path to anndata containing phenotyping info.
    output : str
        File path to output.
    output_plot: str or None
        File path to save the plotting image.
    """
    warnings.simplefilter('ignore')

    with open(inputs, 'r') as param_handler:
        params = json.load(param_handler)

    adata = read_h5ad(anndata)

    if 'spatial' not in adata.obsm:
        try:
            adata.obsm['spatial'] = adata.obs[
                ['X_centroid', 'Y_centroid']].values
        except Exception as e:
            print(e)

    options = params['analyses']['options']

    for k, v in options.items():
        if not isinstance(v, str):
            continue

        if v in ('', 'none'):
            options[k] = None
            continue

        if k == 'genes':    # for spatial_autocorr and sepal
            options[k] = [e.strip() for e in v.split(',')]
        elif k == 'radius':    # for spatial_neighbors
            options[k] = ast.literal_eval(v)
        elif k == 'interactions':    # for ligrec
            options[k] = pd.read_csv(v, sep="\t")
        elif k == 'max_neighs':
            options[k] = int(v)      # for sepal

    if output_plot:
        plotting_options = params['analyses']['plotting_options']
        for k, v in plotting_options.items():
            if not isinstance(v, str):
                continue

            if v in ('', 'none'):
                plotting_options[k] = None
                continue

            if k == 'figsize':
                options[k] = ast.literal_eval(v)
            elif k in ('palette', 'score', 'source_groups', 'target_groups'):
                options[k] = [e.strip() for e in v.split(',')]
            elif k == 'means_range':        # ligrec
                v = v.strip()
                if v[0] == '(':
                    v = v[1:]
                if v[-1] == ')':
                    v = v[:-1]
                options[k] = tuple([float(e.strip()) for e in v.split(',', 1)])

    sq.pl.spatial_scatter(
        adata=adata,
        shape=options['shape'],
        color=options['color'],
        groups=options['groups'],
        use_raw=options['use_raw'],
        size=options['size'],
        # need to figure out scale and scalebar
    )


if __name__ == '__main__':
    aparser = argparse.ArgumentParser()
    aparser.add_argument("-i", "--inputs", dest="inputs", required=True)
    aparser.add_argument("-e", "--output", dest="output", required=True)
    aparser.add_argument("-a", "--anndata", dest="anndata", required=True)
    aparser.add_argument(
        "-p", "--output_plot", dest="output_plot", required=False)

    args = aparser.parse_args()

    main(args.inputs, args.anndata, args.output, args.output_plot)
