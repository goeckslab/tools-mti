import argparse
import json
import pickle
import os
import time
import warnings

import numpy as np
import skimage.io
import matplotlib.pyplot as plt

from cellpose import models, plot
from model_unpickler import SafeUnpickler


setattr(pickle, 'Unpickler', SafeUnpickler)


def main(inputs, img_path, output_dir):
    """
    Parameter
    ---------
    inputs : str
        File path to galaxy tool parameter
    img_path : str
        File path for the input image
    output_dir : str
        Folder to save the outputs.
    """
    warnings.simplefilter('ignore')

    with open(inputs, 'r') as param_handler:
        params = json.load(param_handler)

    options = params.pop('options')

    # handle mxnet option 
    if params['check_mkl']:
        mkl_enabled = models.check_mkl((not params['mxnet']))
    else:
        mkl_enabled = True

    if mkl_enabled and params['mkldnn']:
        os.environ["MXNET_SUBGRAPH_BACKEND"]="MKLDNN"
    else:
        os.environ["MXNET_SUBGRAPH_BACKEND"]=""

    channels = [params['chan'], params['chan2']]

    use_gpu = params['use_gpu']
    device, gpu = models.assign_device((not params['mxnet']), use_gpu)

    pretrained_model = params['pretrained_model']
    bacterial = 'bact' in pretrained_model

    net_avg=(not options['fast_mode'] and not options['no_net_avg'])
    
    omni = options['omni']
    if 'omni' in pretrained_model:
        omni = True
    
    tic = time.time()

    if params['mxnet']:
        if pretrained_model == 'cyto2':
            print('cyto2 model not available in mxnet, using cyto model')
            pretrained_model = 'cyto'
        if bacterial:
            print('bacterial models not available in mxnet, using pytorch')
            params['mxnet'] = False

    # omni changes not implemented for mxnet. Full parity for cpu/gpu in pytorch. 
    if omni and params['mxnet']:
        print('>>>> omni only implemented in pytorch.')
        print('>>>> omni set to false.')
        omni = False

    # For now, omni version is not compatible with 3D. WIP. 
    if omni and options['do_3D']:
        print('>>>> omni not yet compatible with 3D segmentation.')
        print('>>>> omni set to false.')
        omni = False

    if not bacterial:                
        model = models.Cellpose(
            gpu=gpu, device=device, model_type=pretrained_model,
            torch=(not params['mxnet']),omni=omni,
            net_avg=net_avg)
    else:
        cpmodel_path = models.model_path(pretrained_model, 0, True)
        model = models.CellposeModel(
            gpu=gpu, device=device, pretrained_model=cpmodel_path,
            torch=True, nclasses=3, omni=omni,
            net_avg=False)

    # handle diameters
    diameter = options['diameter']
    if diameter ==0:
        diameter = None

    img = skimage.io.imread(img_path)

    print(f"Image shape: {img.shape}")
    # transpose to Ly x Lx x nchann and reshape based on channels
    # if img_format.endswith('tiff') and params['channel_first']:
    #     img = np.transpose(img, (1, 2, 0))
    #     img = transforms.reshape(img, channels=channels)
    #     channels = [1, 2]
    masks, flows, styles, diams = model.eval(
        img, channels=channels, diameter=diameter,
        do_3D=options['do_3D'],
        net_avg=net_avg)
                                             )

    # save masks to tiff
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        skimage.io.imsave(os.path.join(output_dir, 'cp_masks.tif'),
                          masks.astype(np.uint16))

    # make segmentation show #
    if params['show_segmentation']:
        img = skimage.io.imread(img_path)
        # uniform image
        # if img_format.endswith('tiff') and params['channel_first']:
        #     img = np.transpose(img, (1, 2, 0))
        #     img = transforms.reshape(img, channels=channels)
        #     channels = [1, 2]

        maski = masks
        flowi = flows[0]
        fig = plt.figure(figsize=(12, 3))
        # can save images (set save_dir=None if not)
        plot.show_segmentation(fig, img, maski, flowi, channels=channels)
        fig.savefig(os.path.join(output_dir, 'segm_show.png'), dpi=300)
        plt.close(fig)


if __name__ == '__main__':
    aparser = argparse.ArgumentParser()
    aparser.add_argument("-i", "--inputs", dest="inputs", required=True)
    aparser.add_argument("-p", "--img_path", dest="img_path")
    aparser.add_argument("-O", "--output_dir", dest="output_dir")
    args = aparser.parse_args()

    main(args.inputs, args.img_path, args.output_dir)
