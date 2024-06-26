import argparse
import json
import warnings

import pandas as pd
import scimap as sm
from anndata import read_h5ad


def main(inputs, anndata, output):
    """
    Parameter
    ---------
    inputs : str
        File path to galaxy tool parameter.
    anndata : str
        File path to anndata containing phenotyping info.
    output : str
        File path to output.
    """
    warnings.simplefilter('ignore')

    with open(inputs, 'r') as param_handler:
        params = json.load(param_handler)

    adata = read_h5ad(anndata)

    tool = params['analyses']['selected_tool']
    tool_func = getattr(sm.tl, tool)

    options = params['analyses']['options']

    # tool specific pre-processing
    if tool == 'cluster':
        options['method'] = params['analyses']['method']
        subset_genes = options.pop('subset_genes')
        if subset_genes:
            options['subset_genes'] = \
                [x.strip() for x in subset_genes.split(',')]
        sub_cluster_group = options.pop('sub_cluster_group')
        if sub_cluster_group:
            options['sub_cluster_group'] = \
                [x.strip() for x in sub_cluster_group.split(',')]
    elif tool == 'spatial_lda':
        max_weight_assignment = options.pop('max_weight_assignment')

    for k, v in options.items():
        if v == '':
            options[k] = None

    # tool execution
    tool_func(adata, **options)

    # spatial LDA post-processing
    if tool == 'spatial_lda':

        if max_weight_assignment:
            # assign cell to a motif based on maximum weight
            adata.uns['spatial_lda']['neighborhood_motif'] = \
                adata.uns['spatial_lda'].idxmax(axis=1)

            # merge motif assignment into adata.obs
            adata.obs = pd.merge(
                adata.obs,
                adata.uns['spatial_lda']['neighborhood_motif'],
                left_index=True,
                right_index=True
            )

        # write out LDA results as tabular files
        # so they're accessible to Galaxy users
        adata.uns['spatial_lda'].reset_index().to_csv(
            'lda_weights.txt', sep='\t', index=False)
        adata.uns['spatial_lda_probability'].T.reset_index(
            names='motif').to_csv(
                'lda_probabilities.txt', sep='\t', index=False)

        if 'spatial_lda_model' in adata.uns:
            adata.uns.pop('spatial_lda_model')

    adata.write(output)


if __name__ == '__main__':
    aparser = argparse.ArgumentParser()
    aparser.add_argument("-i", "--inputs", dest="inputs", required=True)
    aparser.add_argument("-e", "--output", dest="output", required=True)
    aparser.add_argument("-a", "--anndata", dest="anndata", required=True)

    args = aparser.parse_args()

    main(args.inputs, args.anndata, args.output)
