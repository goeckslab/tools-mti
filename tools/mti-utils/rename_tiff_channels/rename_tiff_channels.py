import sys
import pandas
import click
import tifffile
import bioformats
import javabridge

@click.command()
@click.argument('input_tiff', type=str)
@click.argument('channels_csv', type=str)
@click.argument('output_name', type=str)
def rename_tiff_channels(input_tiff, channels_csv, output_name):
    '''
    Tool for renaming the channel XML elements in an OME-TIFF image
    '''
    # Process channels_csv
    channels_df = pandas.read_csv(channels_csv)

    # The JavaBridge approaches..
    javabridge.start_vm(class_path=bioformats.JARS)

    # Get the file metadata and pixel values
    with tifffile.TiffFile(input_tiff) as tif:
        md = tif.ome_metadata

    # Create the OMEXML object
    xml_obj = bioformats.OMEXML(md)

    # Iterate through the channels and change their names
    for channel_index in range(xml_obj.image().Pixels.get_channel_count()):
        new_name = channels_df.loc[channel_index, 'marker_name']
        xml_obj.image().Pixels.channel(channel_index).set_Name(new_name)

    # This is a header for writing the XML out
    header = '<?xml version="1.0" encoding="UTF-8"?><!-- Warning: this comment is an OME-XML metadata block, which contains crucial dimensional parameters and other important metadata. Please edit cautiously (if at all), and back up the original data before doing so. For more information, see the OME-TIFF web site: https://docs.openmicroscopy.org/latest/ome-model/ome-tiff/. -->'
    # Set up the xml to write and do that
    xml_str = xml_obj.to_xml()
    removed_ome = '<'.join(xml_str.split('<ome:'))
    removed_ome = '/'.join(removed_ome.split('/ome:'))
    with open(output_name, "w") as out_tiff:
        out_tiff.write(header+removed_ome)

    # Slay the evil JavaBridge
    javabridge.kill_vm()


if __name__ == "__main__":
    rename_tiff_channels()
