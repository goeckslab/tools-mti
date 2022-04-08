import argparse
import warnings

import pandas as pd
import scimap as sm
from anndata import read_h5ad


def main(
    adata,
    output,
    gating_workflow,
    gating_workflow_ext,
    manual_gates=None,
    manual_gates_ext=None,
    rescale_plots=False
):
    """
    Parameter
    ---------
    adata : str
        File path to the input AnnData.
    output : str
        File path to the output AnnData.
    gating_workflow : str
        File path to the gating workflow.
    gating_workflow_ext : str
        Datatype for gating workflow, either 'csv' or 'tabular'.
    manual_gates : str
        File path to the munual gating.
    manual_gates_ext : str
        Datatype for munual gate, either 'csv' or 'tabular'.
    rescale_plots : boolean
        Save plots from rescaling.
    """
    warnings.simplefilter('ignore')

    adata = read_h5ad(adata)
    # Rescale data
    if manual_gates:
        sep = ',' if manual_gates_ext == 'csv' else '\t'
        manual_gates = pd.read_csv(manual_gates, sep=sep)

    adata = sm.pp.rescale(adata, gate=manual_gates, save_fig=rescale_plots)

    # Phenotype cells
    # Load the gating workflow
    sep = ',' if gating_workflow_ext == 'csv' else '\t'
    phenotype = pd.read_csv(gating_workflow, sep=sep)
    adata = sm.tl.phenotype_cells(adata, phenotype=phenotype, label="phenotype")

    # Summary of the phenotyping
    print(adata.obs['phenotype'].value_counts())

    adata.write(output)


if __name__ == '__main__':
    aparser = argparse.ArgumentParser()
    aparser.add_argument("-a", "--adata", dest="adata", required=True)
    aparser.add_argument("-o", "--output", dest="output", required=True)
    aparser.add_argument("-g", "--gating_workflow", dest="gating_workflow", required=True)
    aparser.add_argument("-s", "--gating_workflow_ext", dest="gating_workflow_ext", required=True)
    aparser.add_argument("-m", "--manual_gates", dest="manual_gates", required=False)
    aparser.add_argument("-S", "--manual_gates_ext", dest="manual_gates_ext", required=False)
    aparser.add_argument("-p", "--rescale_plots", dest="rescale_plots", action="store_true",
                         default=False, required=False)

    args = aparser.parse_args()

    main(args.adata, args.output, args.gating_workflow,
         args.gating_workflow_ext, args.manual_gates,
         args.manual_gates_ext, args.rescale_plots)
