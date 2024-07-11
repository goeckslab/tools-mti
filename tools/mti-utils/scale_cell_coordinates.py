import argparse
import json
import warnings

import anndata as ad


def main(inputs, output):

    """
    inputs : str
        File path to galaxy tool JSON inputs config file
    output: str
        File path to save the output h5ad file
    """
    warnings.simplefilter('ignore')

    # read inputs JSON
    with open(inputs, 'r') as param_handler:
        params = json.load(param_handler)

    # read input anndata file
    adata = ad.read_h5ad(params['anndata'])

    # scale coords
    unit = params['unit']
    new_col_names = []
    for c in [params['x_coord'], params['y_coord']]:
        scaled_col_name = f'{c}_{unit}'
        adata.obs[scaled_col_name] = adata.obs[c] * params['resolution']
        new_col_names.append(scaled_col_name)

    # overwrite adata.obsm['spatial'] with scaled coordinates
    adata.obsm['spatial'] = adata.obs[scaled_col_name].values

    # write out anndata to h5ad file
    adata.write_h5ad(output)


if __name__ == '__main__':

    aparser = argparse.ArgumentParser()
    aparser.add_argument(
        "-i", "--inputs", dest="inputs", required=True)
    aparser.add_argument(
        "-o", "--output", dest="output", required=False)

    args = aparser.parse_args()

    main(args.inputs, args.output)
