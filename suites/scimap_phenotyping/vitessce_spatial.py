import argparse
import json
import warnings

import scanpy as sc
from pathlib import Path
from anndata import read_h5ad
from vitessce import (
    VitessceConfig,
    Component as cm,
    AnnDataWrapper,
    OmeTiffWrapper,
    MultiImageWrapper,
)


def main(inputs, output, image, anndata=None, masks=None):
    """
    Parameter
    ---------
    inputs : str
        File path to galaxy tool parameter.
    output : str
        Output folder for saving web content.
    image : str
        File path to the OME Tiff image.
    anndata : str
        File path to anndata containing phenotyping info.
    masks : str
        File path to the image masks.
    """
    warnings.simplefilter('ignore')

    with open(inputs, 'r') as param_handler:
        params = json.load(param_handler)

    vc = VitessceConfig(name=None, description=None)
    dataset = vc.add_dataset()
    image_wrappers = [OmeTiffWrapper(img_path=image, name='OMETIFF')]
    if masks:
        image_wrappers.append(
            OmeTiffWrapper(img_path=masks, name='MASKS', is_bitmask=True)
        )
    dataset.add_object(MultiImageWrapper(image_wrappers))

    status = vc.add_view(dataset, cm.STATUS)
    spatial = vc.add_view(dataset, cm.SPATIAL)
    lc = vc.add_view(dataset, cm.LAYER_CONTROLLER)

    if not anndata:
        vc.layout(status / lc | spatial)
        config_dict = vc.export(to='files', base_url='http://localhost', out_dir=output)
        with open(Path(output).joinpath('config.json'), 'w') as f:
            json.dump(config_dict, f, indent=4)
        return

    adata = read_h5ad(anndata)

    params = params['do_phenotyping']
    embedding = params['scatterplot_embeddings']['embedding']
    embedding_options = params['scatterplot_embeddings']['options']
    if embedding == 'umap':
        sc.pp.neighbors(adata, **embedding_options)
        sc.tl.umap(adata)
        mappings_obsm = 'X_umap'
        mappings_obsm_name = "UMAP"
    elif embedding == 'tsne':
        sc.tl.tsne(adata, **embedding_options)
        mappings_obsm = 'X_tsne'
        mappings_obsm_name = "tSNE"
    else:         # pca
        sc.tl.pca(adata, **embedding_options)
        mappings_obsm = 'X_pca'
        mappings_obsm_name = "PCA"

    adata.obsm['XY_centroid'] = adata.obs[['X_centroid', 'Y_centroid']].values

    cell_set_obs = params['phenotype_factory']['phenotypes']
    if not isinstance(cell_set_obs, list):
        cell_set_obs = [x.strip() for x in cell_set_obs.split(',')]
    cell_set_obs_names = [obj[0].upper() + obj[1:] for obj in cell_set_obs]
    dataset.add_object(
        AnnDataWrapper(
            adata,
            mappings_obsm=[mappings_obsm],
            mappings_obsm_names=[mappings_obsm_name],
            spatial_centroid_obsm='XY_centroid',
            cell_set_obs=cell_set_obs,
            cell_set_obs_names=cell_set_obs_names,
            expression_matrix="X"
        )
    )

    cellsets = vc.add_view(dataset, cm.CELL_SETS)
    scattorplot = vc.add_view(dataset, cm.SCATTERPLOT, mapping=mappings_obsm_name)
    heatmap = vc.add_view(dataset, cm.HEATMAP)
    genes = vc.add_view(dataset, cm.GENES)
    cell_set_sizes = vc.add_view(dataset, cm.CELL_SET_SIZES)
    cell_set_expression = vc.add_view(dataset, cm.CELL_SET_EXPRESSION)
    vc.layout(
        (status / genes / cell_set_expression) |
        (cellsets / lc / scattorplot) |
        (cell_set_sizes / heatmap / spatial)
    )
    config_dict = vc.export(to='files', base_url='http://localhost', out_dir=output)

    with open(Path(output).joinpath('config.json'), 'w') as f:
        json.dump(config_dict, f, indent=4)


if __name__ == '__main__':
    aparser = argparse.ArgumentParser()
    aparser.add_argument("-i", "--inputs", dest="inputs", required=True)
    aparser.add_argument("-e", "--output", dest="output", required=True)
    aparser.add_argument("-g", "--image", dest="image", required=True)
    aparser.add_argument("-a", "--anndata", dest="anndata", required=False)
    aparser.add_argument("-m", "--masks", dest="masks", required=False)

    args = aparser.parse_args()

    main(args.inputs, args.output, args.image, args.anndata, args.masks)
