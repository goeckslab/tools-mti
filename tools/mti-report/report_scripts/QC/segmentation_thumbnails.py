#!/usr/bin/env python

import os
import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
import tifffile
import zarr
from PIL import Image

from skimage.filters import threshold_otsu
from skimage.measure import regionprops
from skimage.color import label2rgb
from skimage.draw import rectangle_perimeter
from skimage.morphology import remove_small_objects

def get_args():

    parser = argparse.ArgumentParser(add_help = True,
    description = "Produces a thumbnail of the DAPI channel from a registered ome.tiff and a full-frame image for context")

    parser.add_argument("-f", "--file", 
    help = "Registered ome.tiff file", 
    required = True)

    parser.add_argument("-s", "--segmentation",
    help = "Cell segmentation mask in tiff format",
    required = True)

    parser.add_argument("-l", "--level", 
    help = "OME-tiff pyramid level (default is 0)", 
    required = False,
    default = 0)

    parser.add_argument("-c", "--channel", 
    help = "Channel number (default is 1)", 
    required = False,
    default = 1)

    return parser.parse_args()


def parse_channel(image, level, channel):

    '''
    Reads registered ome.tiff file as array, subsets to pyramid level and channel of interest
    '''
    
    print("Parsing input image \n")

    level = int(level)
    channel = int(channel - 1)
    
    with tifffile.TiffFile(image, is_ome=False) as tiff:

        level_array = zarr.open(tiff.aszarr(series = 0, level = level))
        channel_array = level_array.get_orthogonal_selection((channel, slice(None), slice(None)))

    print(f"Level {level} shape: ", level_array.shape)
    print(f"Level {level} channel {channel} shape: ", channel_array.shape, "\n")

    return channel_array


def remove_background(channel_array):

    print("Removing the image background \n")

    threshold = threshold_otsu(channel_array)
    binary_mask = (channel_array > threshold)
    binary_mask = remove_small_objects(binary_mask).astype(int)
    
    return binary_mask


def gather_tissue_properties(channel_array, binary_mask):

    '''
    Find the tissue center and area. Determines cropping radius from tissue center 
    based on tissue area
    '''

    properties = regionprops(binary_mask, channel_array)
    tissue_center = properties[0].centroid
    tissue_area = properties[0].area
    coords = properties[0].coords

    print("Tissue center is: ", tissue_center)
    print("Tissue area (number of pixels): ", tissue_area, "\n")

    # calculate a reasonable crop diameter
    # going w/ 40% of tissue len if it was a box for absolutely no good reason
    # this should probably be a tweakable param
    tissue_len = np.sqrt(tissue_area)
    crop_radius = (tissue_len * 0.40) // 2

    return tissue_center, crop_radius, coords


def get_crops(tissue_coords, crop_radius):

    """
    Takes tissue coordinates array and a cropping radius as inputs,
    calculates 6 tiles by which to split the tissue, finds each tile's center,
    and computes coordinates of rectangle corners to crop tiles from images and masks by 
    """

    # y,x to convert origin from upper left to lower left 
    y,x = np.moveaxis(tissue_coords, 1, 0)
    
    # Make quantile cuts by both x and y coordinates, this will result in 6 tiles
    df = pd.DataFrame({'x':x,'y':y})
    df = df.assign(
    X_bin = pd.qcut(df.x, 3, labels = ["x1", "x2", "x3"]),
    Y_bin = pd.qcut(df.y, 2, labels = ["y1", "y2"])
    )

    # merge x and y cuts into 2D tiles
    df = df.assign(tile = pd.Categorical(df.filter(regex='_bin').apply(tuple, 1)))

    # get cropping coordinates by projecting a radius around the centroid of each tile
    groups = df.groupby("tile")
    coordinate_dict = {}
    for i, (name, group) in enumerate(groups):

        centroid = (sum(group["x"]) / len(group), sum(group["y"]) / len(group))

        corners = [
            int(centroid[0] - crop_radius),
            int(centroid[0] + crop_radius),
            int(centroid[1] - crop_radius),
            int(centroid[1] + crop_radius)
        ]

        coordinate_dict[i] = [corners, centroid]

    return coordinate_dict


def process_segmentation_mask(segmentation_mask, crop_coords):

    '''
    Reads segmentation mask, crops according to coordinates given and outputs as PIL object
    '''

    with tifffile.TiffFile(segmentation_mask, is_ome=False) as tiff:

        # read cropped mask array
        mask_cropped = tiff.asarray()[crop_coords[2]:crop_coords[3], crop_coords[0]:crop_coords[1]]

    seg = Image.fromarray(mask_cropped)

    return seg

def crop_and_render(channel_array, binary_mask, crop_dict, segmentation_mask, output):
    
    '''
    Turns single-channel array into RGB image, calculates cropping coordinates based on 
    tissue center and cropping radius, calls function to crop and process segmentation mask,
    and finally renders the cropped thumbnail with and without the segmentation mask overlaid.
    Also produces a full-frame image overlaid with a rectangle that shows the cropped region for context.
    '''

    # convert channel array to RGB w/ background removed
    rgb_image = label2rgb(binary_mask, channel_array, colors = ['black', 'white'], alpha = 0.1, bg_label = 0)
    
    # hacky way to filter out tiles that have a small proportion of tissue 
    # currently tile must contain more than 5% tissue pixels. This should probably be a user param
    # maybe consider making crop_dict a dataframe so after purging bad tiles, can reset index
    drop_tiles = []
    for key in crop_dict:
        rect = crop_dict[key][0]
        mask_tile = binary_mask[rect[2]:rect[3],rect[0]:rect[1]]
        proportion_tissue = np.sum(mask_tile) / (len(mask_tile)**2)
        print(key, proportion_tissue)
        if proportion_tissue <= 0.05:
            drop_tiles.append(key)

    for tile in drop_tiles:
        del crop_dict[tile]

    print(f"Rendering and saving images:")

    for key in crop_dict:

        rect = crop_dict[key][0]
        cent = crop_dict[key][1]

        # cropped image w/ segmetation mask
        seg = process_segmentation_mask(segmentation_mask, crop_coords = rect)

        # name each segmentation tile w/ tile index
        seg_out = os.path.join(output, f"segmentation_{key}.png")
        print(f"\t {seg_out}")

        fig,ax = plt.subplots()
        ax.imshow(rgb_image[rect[2]:rect[3],rect[0]:rect[1]])
        ax.contour(seg, [0.5], linewidths=0.5, colors='y')
        at = AnchoredText(key, prop=dict(size=15,color='yellow'), frameon=True, loc='upper left')
        at.patch.set_edgecolor('yellow')
        at.patch.set_facecolor('black')
        ax.add_artist(at)
        fig.savefig(seg_out, dpi = 500)

    # uncropped image with rectangle of cropped region for context
    context_out = os.path.join(output, "context_image.png")
    print(f"\t {context_out}")

    fig,ax = plt.subplots()
    ax.imshow(rgb_image)
    for key in crop_dict:
        rect = crop_dict[key][0]
        cent = crop_dict[key][1]
        row, col = rectangle_perimeter(start=(rect[0],rect[2]), end=(rect[1], rect[3]))
        ax.plot(row, col, "--y")
        ax.annotate(key, cent, color = 'yellow')
    fig.savefig(context_out, dpi = 300)

    print("done! \n")

    return 0


def main():

    args = get_args()

    channel_array = parse_channel(args.file, level = args.level, channel = args.channel)

    binary_mask = remove_background(channel_array)

    tissue_center, crop_radius, coords = gather_tissue_properties(channel_array, binary_mask)

    crops = get_crops(coords, crop_radius)

    return crop_and_render(channel_array = channel_array, 
                           binary_mask = binary_mask,
                           crop_dict = crops,
                           segmentation_mask = args.segmentation,
                           output = os.getcwd())

if __name__ == '__main__':

    sys.exit(main())