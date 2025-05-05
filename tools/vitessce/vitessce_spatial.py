import argparse
import json
import warnings
from os.path import isdir, join
from pathlib import Path
import urllib.parse

import scanpy as sc
from anndata import read_h5ad
from vitessce import (
    AnnDataWrapper,
    Component as cm,
    MultiImageWrapper,
    OmeTiffWrapper,
    VitessceConfig,
)
from vitessce.data_utils import (
    optimize_adata,
    VAR_CHUNK_SIZE,
)


def main(inputs, output, image, offsets=None, anndata=None, masks=None, config_path=None):
    """
    Parameter
    ---------
    inputs : str
        File path to galaxy inputs config file.
    output : str
        Output folder for saving web content.
    image : str
        File path to the OME Tiff image.
    anndata : str
        File path to anndata containing phenotyping info.
    masks : str
        File path to the image masks.
    config_path : str
        File path to the config file containing galaxy_url and dataset_id.
    """
    warnings.simplefilter('ignore')

    with open(inputs, 'r') as param_handler:
        params = json.load(param_handler)

    with open(config_path) as conf_fh:
        config = json.load(conf_fh)

    galaxy_url = config["galaxy_url"]
    dataset_id = config["dataset_id"]

    # initialize vitessce config and add OME-TIFF image, and masks if specified
    vc = VitessceConfig(schema_version="1.0.17", name=None, description=None)
    dataset = vc.add_dataset()

    # FIXME: grab offsets file for faster display. NEED TO TEST
    image_wrappers = [OmeTiffWrapper(img_path=image, offsets_path=offsets, name='OMETIFF')]
    if masks:
        image_wrappers.append(
            OmeTiffWrapper(img_path=masks, name='MASKS', is_bitmask=True)
        )
    dataset.add_object(MultiImageWrapper(image_wrappers))

    # set relative view sizes (w,h), full window dims are 12x12
    # if no anndata file, image and layer view can take up whole window
    if not anndata:
        spatial_dims = (9, 12)
        lc_dims = (3, 12)
    else:
        spatial_dims = (6, 6)
        lc_dims = (3, 6)

    # add views for the images, and the layer/channels selector
    spatial = vc.add_view(
        view_type=cm.SPATIAL,
        dataset=dataset,
        w=spatial_dims[0],
        h=spatial_dims[1])

    lc = vc.add_view(
        view_type=cm.LAYER_CONTROLLER,
        dataset=dataset,
        w=lc_dims[0],
        h=lc_dims[1])

    # Build the prefix that Vitessce should use
    display_prefix = (f"{galaxy_url}/api/datasets/{dataset_id}/display?filename=")
    
    # if no anndata file, export the config with these minimal components
    if not anndata:
        vc.layout(lc | spatial)
        config_dict = vc.export(
            to='files',
            base_url=display_prefix,
            out_dir=output)
        with open(Path(output).joinpath('config.json'), 'w') as f:
            json.dump(config_dict, f, indent=4)
        return

    # read anndata file, compute embeddings
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
    else:
        sc.tl.pca(adata, **embedding_options)
        mappings_obsm = 'X_pca'
        mappings_obsm_name = "PCA"

    # Add spatial coords to obsm, although uncertain if this is needed
    # FIXME: provide options for alternative coordinate colnames
    adata.obsm['spatial'] = adata.obs[['X_centroid', 'Y_centroid']].values

    # parse list of obs columns to use as cell type labels
    cell_set_obs = params['phenotype_factory']['phenotypes']
    if not isinstance(cell_set_obs, list):
        cell_set_obs = [x.strip() for x in cell_set_obs.split(',')]

    # write anndata out as zarr hierarchy
    zarr_filepath = join("data", "adata.zarr")
    if not isdir(zarr_filepath):
        adata = optimize_adata(
            adata,
            obs_cols=cell_set_obs,
            obsm_keys=[mappings_obsm, 'spatial'],
            optimize_X=True
        )
        adata.write_zarr(
            zarr_filepath,
            chunks=[adata.shape[0], VAR_CHUNK_SIZE]
        )

    # create a nicer label for the cell types to be displayed on the dashboard
    cell_set_obs_names = [obj[0].upper() + obj[1:] for obj in cell_set_obs]

    # add anndata zarr to vitessce config
    dataset.add_object(
        AnnDataWrapper(
            adata_path=zarr_filepath,
            obs_feature_matrix_path="X",    # FIXME: provide rep options
            obs_set_paths=['obs/' + x for x in cell_set_obs],
            obs_set_names=cell_set_obs_names,
            obs_locations_path='spatial',
            obs_embedding_paths=['obsm/' + mappings_obsm],
            obs_embedding_names=[mappings_obsm_name]
        )
    )

    # add views
    cellsets = vc.add_view(
        view_type=cm.OBS_SETS,
        dataset=dataset,
        w=3,
        h=3)

    scatterplot = vc.add_view(
        view_type=cm.SCATTERPLOT,
        dataset=dataset,
        mapping=mappings_obsm_name,
        w=3,
        h=6)

    heatmap = vc.add_view(
        view_type=cm.HEATMAP,
        dataset=dataset,
        w=3,
        h=3)

    genes = vc.add_view(
        view_type=cm.FEATURE_LIST,
        dataset=dataset,
        w=3,
        h=3)

    cell_set_sizes = vc.add_view(
        view_type=cm.OBS_SET_SIZES,
        dataset=dataset,
        w=3,
        h=3)

    cell_set_expression = vc.add_view(
        view_type=cm.OBS_SET_FEATURE_VALUE_DISTRIBUTION,
        dataset=dataset,
        w=3,
        h=6)

    # define the dashboard layout
    vc.layout(
        (cellsets / genes / cell_set_expression)
        | (lc / scatterplot)
        | (cell_set_sizes / heatmap / spatial)
    )

    # export the config file
    config_dict = vc.export(
        to='files',
        base_url=display_prefix,
        out_dir=output)

    with open(Path(output).joinpath('config.json'), 'w') as f:
        json.dump(config_dict, f, indent=4)


if __name__ == '__main__':
    aparser = argparse.ArgumentParser()
    aparser.add_argument("-i", "--inputs", dest="inputs", required=True)
    aparser.add_argument("-e", "--output", dest="output", required=True)
    aparser.add_argument("-g", "--image", dest="image", required=True)
    aparser.add_argument("-f", "--offsets", dest="offsets", required=False)
    aparser.add_argument("-a", "--anndata", dest="anndata", required=False)
    aparser.add_argument("-m", "--masks", dest="masks", required=False)
    aparser.add_argument("--galaxy_config", dest="config_path", required=True)

    args = aparser.parse_args()

    main(args.inputs, args.output, args.image, args.offsets, args.anndata, args.masks, args.config_path)
