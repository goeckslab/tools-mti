import argparse
import json
import warnings
from os.path import isdir, join
from pathlib import Path

import numpy as np
import pandas as pd
from anndata import read_h5ad
from sklearn.mixture import GaussianMixture
from sklearn.preprocessing import MinMaxScaler
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


# Generate binarized phenotype for a gate
def get_gate_phenotype(g, d):
    dd = d.copy()
    dd = np.where(dd < g, 0, dd)
    warnings.filterwarnings('ignore')
    dd = np.where(dd >= g, 1, dd)
    return dd


def get_gmm_phenotype(data):
    low = np.percentile(data, 0.01)
    high = np.percentile(data, 99.99)
    data = np.clip(data, low, high)

    sum = np.sum(data)
    median = np.median(data)
    data_med = data / sum * median

    data_log = np.log1p(data_med)
    data_log = data_log.reshape(-1, 1)

    scaler = MinMaxScaler(feature_range=(0, 1))
    data_norm = scaler.fit_transform(data_log)

    gmm = GaussianMixture(n_components=2)
    gmm.fit(data_norm)
    gate = np.mean(gmm.means_)

    return get_gate_phenotype(gate, np.ravel(data_norm))


def main(inputs, output, image, anndata, offsets=None, masks=None, config_path=None):
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
    config_path : str
        File path to the config containing galaxy_url and dataset_id.
    """
    warnings.simplefilter('ignore')

    with open(inputs, 'r') as param_handler:
        params = json.load(param_handler)

    with open(config_path) as conf_fh:
        config = json.load(conf_fh)

    marker = params['marker'].strip()
    from_gate = params['from_gate']
    to_gate = params['to_gate']
    increment = params['increment']
    x_coordinate = params['x_coordinate'].strip() or 'X_centroid'
    y_coordinate = params['y_coordinate'].strip() or 'Y_centroid'

    galaxy_url = config["galaxy_url"]
    dataset_id = config["dataset_id"]

    # Build the prefix that Vitessce should use
    display_prefix = (f"{galaxy_url}/api/datasets/{dataset_id}/display?filename=")

    adata = read_h5ad(anndata)

    # If no raw data is available make a copy
    if adata.raw is None:
        adata.raw = adata

    # Copy of the raw data if it exisits
    if adata.raw is not None:
        adata.X = adata.raw.X

    data = pd.DataFrame(
        adata.X,
        columns=adata.var.index,
        index=adata.obs.index
    )
    marker_values = data[[marker]].values
    marker_values_log = np.log1p(marker_values)

    # Identify the list of increments
    gate_names = []
    for num in np.arange(from_gate, to_gate, increment):
        num = round(num, 3)
        key = marker + '--' + str(num)
        adata.obs[key] = get_gate_phenotype(num, marker_values_log)
        gate_names.append(key)

    adata.obs['GMM_auto'] = get_gmm_phenotype(marker_values)
    gate_names.append('GMM_auto')

    adata.obsm['spatial'] = adata.obs[[x_coordinate, y_coordinate]].values

    # initialize vitessce config and add OME-TIFF image
    vc = VitessceConfig(schema_version="1.0.17", name=None, description=None)
    dataset = vc.add_dataset()
    image_wrappers = [OmeTiffWrapper(img_path=image, offsets_path=offsets, name='OMETIFF')]
    if masks:
        image_wrappers.append(
            OmeTiffWrapper(img_path=masks, name='MASKS', is_bitmask=True)
        )
    dataset.add_object(MultiImageWrapper(image_wrappers))

    # write anndata out as zarr hierarchy
    zarr_filepath = join("data", "adata.zarr")
    if not isdir(zarr_filepath):
        adata = optimize_adata(
            adata,
            obs_cols=gate_names,
            obsm_keys=['spatial'],
            optimize_X=True
        )
        adata.write_zarr(
            zarr_filepath,
            chunks=[adata.shape[0], VAR_CHUNK_SIZE]
        )

    # add anndata zarr to vitessce config
    dataset.add_object(
        AnnDataWrapper(
            adata_path=zarr_filepath,
            obs_feature_matrix_path="X",    # FIXME: provide rep options
            obs_set_paths=['obs/' + x for x in gate_names],
            obs_set_names=gate_names,
            obs_locations_path='spatial'
        )
    )

    # add views
    spatial = vc.add_view(
        view_type=cm.SPATIAL,
        dataset=dataset,
        w=6,
        h=12)

    cellsets = vc.add_view(
        view_type=cm.OBS_SETS,
        dataset=dataset,
        w=3,
        h=6)

    lc = vc.add_view(
        view_type=cm.LAYER_CONTROLLER,
        dataset=dataset,
        w=3,
        h=9)

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
        h=3)

    # define the dashboard layout
    vc.layout(
        (cellsets / genes / cell_set_expression)
        | (cell_set_sizes / lc)
        | (spatial)
    )

    # export config file
    config_dict = vc.export(to='files', base_url=display_prefix, out_dir=output)

    with open(Path(output).joinpath('config.json'), 'w') as f:
        json.dump(config_dict, f, indent=4)


if __name__ == '__main__':
    aparser = argparse.ArgumentParser()
    aparser.add_argument("-i", "--inputs", dest="inputs", required=True)
    aparser.add_argument("-e", "--output", dest="output", required=True)
    aparser.add_argument("-g", "--image", dest="image", required=True)
    aparser.add_argument("-a", "--anndata", dest="anndata", required=True)
    aparser.add_argument("-f", "--offsets", dest="offsets", required=False)
    aparser.add_argument("-m", "--masks", dest="masks", required=False)
    aparser.add_argument("--galaxy_config", dest="config_path", required=True)

    args = aparser.parse_args()

    main(args.inputs, args.output, args.image, args.anndata, args.offsets, args.masks, args.config_path)
