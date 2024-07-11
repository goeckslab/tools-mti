import argparse
import ast
import json
import warnings

import squidpy as sq
from anndata import read_h5ad


def main(inputs, output_plot):

    """
    inputs : str
        File path to galaxy tool JSON inputs config file
    output_plot: str
        File path to save the plotting image
    """
    warnings.simplefilter('ignore')

    # read inputs JSON
    with open(inputs, 'r') as param_handler:
        params = json.load(param_handler)

    # collapse param dict hierarchy, parse inputs
    plot_opts = params.pop('plot_opts')
    legend_opts = params.pop('legend_opts')
    aes_opts = params.pop('aesthetic_opts')
    options = {**params, **plot_opts, **legend_opts, **aes_opts}

    # read input anndata file
    adata_fh = options.pop('anndata')
    adata = read_h5ad(adata_fh)

    # ensure spatial coords in anndata.obsm['spatial']
    # if not, populate with user provided X/Y coord column names
    x, y = options.pop('x_coord'), options.pop('y_coord')
    if 'spatial' not in adata.obsm:
        try:
            adata.obsm['spatial'] = adata.obs[[x, y]].values
        except Exception as e:
            print(e)

    # scan thru tool params,
    # replace None values, and reformat specific parameters
    for k, v in options.items():
        if not isinstance(v, str):
            continue

        if v in ('', 'none'):
            options[k] = None
            continue

        if k == 'groups':
            options[k] = [e.strip() for e in v.split(',')]

        elif k == 'crop_coord':
            # split str on commas into tuple of coords
            # then nest in list (expected by squidpy function)
            options[k] = [tuple([int(e.strip()) for e in v.split(',')])]

        elif k == 'figsize':
            options[k] = ast.literal_eval(v)

    # not exposing this parameter for now. Only useful to change for ST data
    # and can otherwise just be problematic.
    # Explicitly setting to None is necessary to avoid an error
    options['shape'] = None

    # call squidpy spatial scatter function, unpack tool params
    sq.pl.spatial_scatter(
        adata=adata,
        save='image.png',
        **options
    )


if __name__ == '__main__':

    aparser = argparse.ArgumentParser()
    aparser.add_argument(
        "-i", "--inputs", dest="inputs", required=True)
    aparser.add_argument(
        "-p", "--output_plot", dest="output_plot", required=False)

    args = aparser.parse_args()

    main(args.inputs, args.output_plot)
