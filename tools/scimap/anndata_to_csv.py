import argparse
import json
import warnings

import scimap as sm


def main(inputs, anndata, outfile):
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

    adata = read_h5ad(anndata)

    df = sm.hl.scimap_to_csv(
        adata = adata, 
        data_type = params['rep']['representation'],
        CellID = params['cellid'],
    )

    df.to_csv(outfile)


if __name__ == '__main__':
    aparser = argparse.ArgumentParser()
    aparser.add_argument("-i", "--inputs", dest="inputs", required=True)
    aparser.add_argument("-e", "--outfile", dest="outfile", required=True)
    args = aparser.parse_args()

    main(args.inputs, args.outfile)
