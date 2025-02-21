# ------------------------------------------------------------------------------------
# Stripped down and modified from:
# https://github.com/ohsu-comp-bio/ashlar/blob/master/pyramid_upgrade.py
# ------------------------------------------------------------------------------------

import os
import sys 
import argparse
import csv
import xml.etree.ElementTree
from tifffile import tiffcomment


def fix_attrib_namespace(elt):
    """Prefix un-namespaced XML attributes with the tag's namespace."""
    # This fixes ElementTree's inability to round-trip XML with a default
    # namespace ("cannot use non-qualified names with default_namespace option"
    # error). 7-year-old BPO issue here: https://bugs.python.org/issue17088
    # Code inspired by https://gist.github.com/provegard/1381912 .
    if elt.tag[0] == "{":
        uri, _ = elt.tag[1:].rsplit("}", 1)
        new_attrib = {}
        for name, value in elt.attrib.items():
            if name[0] != "{":
                # For un-namespaced attributes, copy namespace from element.
                name = f"{{{uri}}}{name}"
            new_attrib[name] = value
        elt.attrib = new_attrib
    for child in elt:
        fix_attrib_namespace(child)


def main(image_fh, marker_file):
    """
    Parameters
    ---------
    image_fh : str
        File path to the OME Tiff image.
    marker_file : str
        File path to CSV containing marker name information.
    """

    # parse marker file, create list of new marker names
    new_channel_names = []
    with open(marker_file, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            new_channel_names.append(row['marker_name'])

    # read OME-XML metadata and parse down to channels
    xml_ns = {"ome": "http://www.openmicroscopy.org/Schemas/OME/2016-06"}
    root = xml.etree.ElementTree.fromstring(tiffcomment(image_fh))
    image = root.find("ome:Image", xml_ns)
    pixels = image.find("ome:Pixels", xml_ns)
    channels = pixels.findall("ome:Channel", xml_ns)

    # name channels
    for channel, name in zip(channels, new_channel_names):
        channel.attrib["Name"] = name

    # encode new xml and set image metadata
    fix_attrib_namespace(root)
    new_ome_xml = xml.etree.ElementTree.tostring(
        root, 
        encoding='utf-8',
        xml_declaration=True,
        default_namespace=xml_ns["ome"])

    tiffcomment(image_fh, comment=new_ome_xml)

    print("Updated OME-TIFF metadata:")
    print(tiffcomment(image_fh))

if __name__ == '__main__':
    aparser = argparse.ArgumentParser()
    aparser.add_argument("-i", "--image", dest="image", required=True)
    aparser.add_argument("-m", "--markers", dest="markers", required=True)

    args = aparser.parse_args()

    main(args.image, args.markers)