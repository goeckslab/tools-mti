import argparse
import json
import warnings

import scimap as sm


def main(inputs, outfile):
    """
    Parameter
    ---------
    inputs : str
        File path to galaxy tool parameter.

    outfile : str
        File path to estimator.
    """
    warnings.simplefilter('ignore')

    with open(inputs, 'r') as param_handler:
        params = json.load(param_handler)

    image_path = params['image_path']
    drop_markers = params['drop_markers']
    if not drop_markers:
        drop_markers = None
    else:
        drop_markers = [x.strip() for x in drop_markers.split(',')]
    options = params['options']
    for k, v in options.items():
        if v == '':
            options[k] = None

    adata = sm.pp.mcmicro_to_scimap(image_path, drop_markers=drop_markers, **options)

    adata.write(outfile)


if __name__ == '__main__':
    aparser = argparse.ArgumentParser()
    aparser.add_argument("-i", "--inputs", dest="inputs", required=True)
    aparser.add_argument("-e", "--outfile", dest="outfile", required=True)
    args = aparser.parse_args()

    main(args.inputs, args.outfile)
