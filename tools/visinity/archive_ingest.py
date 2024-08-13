import json
import sys

from os.path import isfile

from zipfile import ZipFile


with open(sys.argv[1], 'r') as f:
    arg_dict = json.load(f)

arch_zip = arg_dict['vis_archive']
og_channels = arg_dict['channel_files']
og_segs = arg_dict['label_files']


data_path = '/visinity/minerva_analysis/data'
conf_path = f'{data_path}/config.json'


with ZipFile(arch_zip, 'r') as zippy:
    zippy.extractall(f'{data_path}')

with open(f'{data_path}/names.csv', 'r') as f:
    name_list = f.read().split(',')

with open(conf_path, 'r') as f:
    config = json.load(f)


def reconfigure(channel, seg, old_config):
    subconfig = dict(old_config)
    if not subconfig['segmentation'].startswith(data_path) or not isfile(subconfig['segmentation']):
        subconfig['segmentation'] = seg
    subconfig['channelFile'] = channel
    print(f'\tnew: {subconfig}')
    return subconfig


new_config = {}

for i, name in enumerate(name_list):
    print(f'reconfiguring {name}')
    name = name_list[i]
    channel = og_channels[i]
    seg = og_segs[i]
    # modify config
    print(f'changing config:\n\told: {config[name]}')
    new_config[name] = reconfigure(channel, seg, config[name])

with open(conf_path, 'w') as f:
    f.write(json.dumps(new_config))
