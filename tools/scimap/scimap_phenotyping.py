import argparse
import warnings

import pandas as pd
import scimap as sm
from anndata import read_h5ad


def main(
    adata,
    output,
    log,
    gating_workflow,
    gating_workflow_ext,
    manual_gates=None,
    manual_gates_ext=None,
    random_state=0
):
    """
    Parameter
    ---------
    adata : str
        File path to the input AnnData.
    output : str
        File path to the output AnnData.
    log: bool
        Boolean whether to log the input data prior to rescaling 
    gating_workflow : str
        File path to the gating workflow.
    gating_workflow_ext : str
        Datatype for gating workflow, either 'csv' or 'tabular'.
    manual_gates : str
        File path to the munual gating.
    manual_gates_ext : str
        Datatype for munual gate, either 'csv' or 'tabular'.
    random_state: int
        The seed used by the random number generator for GMM in sm.pp.rescale
    """
    warnings.simplefilter('ignore')

    adata = read_h5ad(adata)
    # Rescale data
    if manual_gates:
        sep = ',' if manual_gates_ext == 'csv' else '\t'
        manual_gates = pd.read_csv(manual_gates, sep=sep)

    adata = sm.pp.rescale(adata, gate=manual_gates, log = log, random_state = random_state)

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
    aparser.add_argument("-l", "--log", dest="log", action="store_true")
    aparser.add_argument("-g", "--gating_workflow", dest="gating_workflow", required=True)
    aparser.add_argument("-s", "--gating_workflow_ext", dest="gating_workflow_ext", required=True)
    aparser.add_argument("-m", "--manual_gates", dest="manual_gates", required=False)
    aparser.add_argument("-S", "--manual_gates_ext", dest="manual_gates_ext", required=False)
    aparser.add_argument("--random_state", dest="random_state", type=int, required=False)

    args = aparser.parse_args()

    if args.log:
        print("\n adata.raw.X will be log1p transformed \n")
  
    main(    
        adata = args.adata,
        output = args.output,
        log = args.log,
        gating_workflow = args.gating_workflow,
        gating_workflow_ext = args.gating_workflow_ext,
        manual_gates = args.manual_gates,
        manual_gates_ext = args.manual_gates_ext,
        random_state = args.random_state
    )