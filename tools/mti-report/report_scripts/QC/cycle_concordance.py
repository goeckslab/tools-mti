#!/usr/bin/env python

import os
import re
import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle
import tifffile
import zarr
from PIL import Image
from skimage.filters import threshold_otsu
from skimage.color import label2rgb

def get_args():

    parser = argparse.ArgumentParser(add_help = True,
    description = "Takes registered OME-TIFF and channel metadata, overlays DAPI channels in different colors to check cycle registration")

    parser.add_argument("-f", "--file", 
    help = "Registered ome.tiff file", 
    required = True)

    parser.add_argument("-m", "--markers", 
    help = "CSV with channel number and marker names", 
    required = True)

    parser.add_argument("-p", "--pattern", 
    help = "Pattern for extracting channels to align (default is DAPI)", 
    required = False,
    default = 'DAPI')

    parser.add_argument("-l", "--level", 
    help = "OME-tiff pyramid level (default is 0)",
    required = False,
    default = 0)

    parser.add_argument("-o", "--output",
    help = "output path (defaults to cwd/QC_registration.png)",
    required = False,
    default = os.path.join(os.getcwd(), "QC_registration.png"))

    return parser.parse_args()


def select_channels(markers, pattern):

    """
    Reads a CSV file containing 1-based channel indices and marker names, 
    return 0-based channel indices for channels that match pattern and marker names
    """

    markers = pd.read_csv(markers)
    
    marker_names = []
    channel_indices = []

    for i,marker in enumerate(markers.iloc[:,1]):

        if re.search(pattern, marker):

            marker_names.append(marker)
            channel_indices.append(i)

    print(f"Extracting channels from metadata: \n")

    for name in marker_names:
        print(name)

    if len(channel_indices) > 12:
        raise ValueError("Currently only supports up to 12 channels")

    return marker_names, channel_indices


def parse_image(image, level, channels):

    '''
    Reads registered ome.tiff file as array, subsets to a single pyramid level
    '''
    
    print("Parsing input image \n")

    level = int(level)

    with tifffile.TiffFile(image, is_ome=False) as tiff:

        level_array = zarr.open(tiff.aszarr(series = 0, level = level))
        channel_arrays = level_array.get_orthogonal_selection((channels, slice(None), slice(None)))

    print(f"Shape of pyramid level with selected channels: {channel_arrays.shape}")

    return channel_arrays


def remove_background(channel_array):

    """
    Simple background removal, doesn't account for small blobs 
    """

    threshold = threshold_otsu(channel_array)
    binary_mask = (channel_array > threshold).astype(int)
    
    return binary_mask


def convert_channels(channels, marker_names, output):

    """
    Convert each channel to RGB, remove bg, change color.
    Element-wise addition of channel arrays so that fully registered areas appear white. 
    """

    # manual color generator that works for now (lame, I know)
    if len(channels) == 3:
        palette = [(1,0,0), (0,1,0), (0,0,1)]
    else:
        palette = [(1,0,0), (1,0.5,0), (1,1,0),
                   (0.5,1,0), (0,1,0), (0,1,0.5),
                   (0,1,1), (0,0.5,1), (0,0,1),
                   (0.5,0,1), (1,0,1), (1,0,0.5)]

    # initialize an array w/ zeros in same shape as image X,Y dims, and Z=3 for RGB
    stacked_channels = np.zeros(shape = (channels.shape[1],channels.shape[2], 3))

    # iterate across channels
    for chan in range(0,channels.shape[0]):

        # create a background mask then convert to RGB, w/color incrementing by channel idx
        print(f"Removing channel background: {marker_names[chan]} \n")
        mask = remove_background(channels[chan,:,:])
        rgb = label2rgb(mask, channels[chan,:,:], colors = [palette[chan]], bg_label = 0)
    
        # px-wise addition of re-colored channels
        stacked_channels = stacked_channels + rgb

    # clip intensity values that are over the max after adding
    print("clipping intensities to (0,1)")
    max_pxs = stacked_channels[:,:] > 1
    stacked_channels[max_pxs] = 1

    # plotting surface
    # gs = GridSpec(1,6)
    # fig = plt.figure(figsize = (6,6))
    # ax1 = fig.add_subplot(gs[:,:-1]) 
    # ax2 = fig.add_subplot(gs[:,-1])
    fig,ax1 = plt.subplots()

    # add the image
    ax1.imshow(stacked_channels)
    ax1.set_axis_off()

    # populate legend and save figure
    handles = [Rectangle((0,0),1,1, color = n) for n in palette[0:len(marker_names)]]
    legend = ax1.legend(handles, marker_names, ncol=1, facecolor = 'black', 
        edgecolor = 'black', framealpha=1, fontsize = 6, loc = 'upper left') # need to adjust automatic location 
    for i,text in enumerate(legend.get_texts()):
        text.set_color(palette[i])
    # ax2.set_axis_off()
    fig.tight_layout()
    fig.savefig(output, dpi = 500)

    return


def main():

    args = get_args()

    marker_names, channel_indices = select_channels(args.markers, args.pattern)

    channel_arrays = parse_image(args.file, args.level, channel_indices)

    convert_channels(channel_arrays, marker_names, args.output)

    return 0


if __name__ == "__main__":
    sys.exit(main())