import argparse
import json
import warnings

import scimap as sm
from anndata import read_h5ad


def main(inputs, outfile):
    """
    Parameters
    ---------
    inputs : str
        File path to galaxy tool parameter.
    anndata : str
        File path to anndata.
    output : str
        File path to output.
    """
    warnings.simplefilter('ignore')

    with open(inputs, 'r') as param_handler:
        params = json.load(param_handler)

    adata = read_h5ad(params['anndata'])

    if params['layer'] == 'x':
        params['layer'] = None

    df = sm.hl.scimap_to_csv(
        adata=adata,
        layer=params['layer'],
        CellID=params['cellid'],
    )

    df.to_csv(outfile, index=False)


if __name__ == '__main__':
    aparser = argparse.ArgumentParser()
    aparser.add_argument("-i", "--inputs", dest="inputs", required=True)
    aparser.add_argument("-e", "--outfile", dest="outfile", required=True)
    args = aparser.parse_args()

    main(args.inputs, args.outfile)
