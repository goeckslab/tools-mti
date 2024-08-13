import csv
import json
import os
import re
import shutil
import sys
import time


from zipfile import ZipFile, ZIP_BZIP2

from minerva_analysis import app


app.config['IS_DOCKER'] = False

# create a non-serving test client to replace http request
requests = app.test_client()

args = sys.argv[1:]
print(f"args: {args}")

input_dict = {}
with open(args[0], 'r') as f:
    input_dict = json.load(f)

out_file = args[1]

name_list = []
with open(input_dict['dataset_names'], 'r') as f:
    _reader = csv.reader(f)
    header = next(_reader)
    assert header[0].lower() == 'names'
    for row_list in _reader:
        row = row_list[0]
        if row in name_list:
            raise ValueError(f'dataset name {row} is a duplicate!'
                             ' please ensure names are unique')
        s = re.sub('[^A-z0-9.~]+', '-', row)
        if s != row:
            row = s
        name_list.append(row)

channel_files = input_dict['channel_files']
label_files = input_dict['label_files']

assert len(name_list) == len(channel_files) == len(label_files)

radius = input_dict['n_radius']
celltype_csv = input_dict['celltype_csv']

# ensure env
data_path = f"{os.getcwd()}/minerva_analysis/data"
config_path = f"{data_path}/config.json"

if not os.path.isdir(data_path):
    os.mkdir(data_path)

if not os.path.isfile(config_path):
    with open(config_path, 'w') as f:
        f.write('{}')

url = 'http://127.0.0.1:8000'

with open(f'{data_path}/names.csv', 'w') as name_csv_file:
    name_csv_file.write(
        ','.join(name_list)
    )

shutil.copyfile(celltype_csv, 'celltype.csv')

payload_data = {
    "action": "Upload",
    "celltype_file": open('celltype.csv', 'rb'),
    "neighborhood_radius": int(radius),
}

for i, name in enumerate(name_list):
    new_quant_csv = f'quant{i}.csv'
    shutil.copyfile(input_dict['quant_csvs'][i], new_quant_csv)
    _sub_d = {
        f'name-{i+1}': name,
        f'channel_file-{i+1}': channel_files[i],
        f'label_file-{i+1}': label_files[i],
        f'csv_file-{i+1}': open(new_quant_csv, 'rb'),
    }
    payload_data.update(_sub_d)


# payload = {
#   'data': {
#     "action": "Upload",
#     "name-1": dataset_name,
#     "neighborhood_radius": int(n_radius),
#     "channel_file-1": channel_file,
#     "label_file-1": label_file,
#     "csv_file-1": open(quant_target, 'rb'),
#     "celltype_file": open(celltype_target, 'rb')
#   }
# }

resp = requests.post(f"{url}/upload", data=payload_data)

resp_regex = 'passVariablesToFrontend\((.*)\);'
json_string = re.search(resp_regex, resp.text).group(1)
combinedChannelData = json.loads(json_string)


# def prep_header(_pre_sorted_headers: List[dict]):
#     header_list = []
#     # x_pos_terms = ['x_centroid', 'cellposition_x', 'x_position', 'x']
#     # y_pos_terms = ['y_centroid', 'cellposition_y', 'y_position', 'y']
#     # phenotype_terms = ['celltype', 'phenotype']

#     _header_prefix = [None, None, None]
#     _sorted_headers = []
#     for header in _pre_sorted_headers[1:]:
#         if header['fullName'].lower() in x_pos_terms:
#             assert (_header_prefix[0] is None)
#             _header_prefix[0] = header
#         elif header['fullName'].lower() in y_pos_terms:
#             assert (_header_prefix[1] is None)
#             _header_prefix[1] = header
#         elif header['fullName'].lower() in phenotype_terms:
#             assert (_header_prefix[2] is None)
#             _header_prefix[2] = header
#         else:
#             _sorted_headers.append(header)
#     _i = 1
#     all_headers = _header_prefix + _sorted_headers
#     all_headers: dict
#     for header in all_headers:
#         header_list.extend(
#             chunkify(header['fullName'], _i)
#         )
#     return header_list


def chunkify(fullname, idx):
    non_norm_markers = [
        'majoraxislength', 'number_neighbors', 'minoraxislength', 'neighbor_1',
        'solidity', 'eccentricity', 'y_position', 'x_position', 'neighbor_2',
        'percent_touching', 'orientation', 'neighbor_4', 'extent', 'cellid',
        'field_col', 'eulernumber', 'neighbor_3', 'neighbor_5', 'perimeter',
        'field_row']

    return [
        {
          "name": f"fullName{idx}",
          "value": fullname
        },
        {
          "name": f"name{idx}",
          "value": fullname
        },
        {
          "name": f"normalize{idx}",
          "value": "off" if fullname.lower() in non_norm_markers else "on"
        },
    ]


x_pos_terms = ['x_centroid', 'cellposition_x', 'x_position', 'x']
y_pos_terms = ['y_centroid', 'cellposition_y', 'y_position', 'y']
phenotype_terms = ['celltype', 'phenotype']

# headers should all be the same, so using [0]
pre_sorted_headers = combinedChannelData[0]["csvHeader"]

# fix issue with filename "*_<digits>" being mistakenly assumed to be 
# non-segmentation files; weirdly channel names dont seem to matter,
# only if it has an underscore and number
for channel in combinedChannelData:
    channel["labelName"] = 'segmentation'

idfield = chunkify(pre_sorted_headers[0]['fullName'], 0)

header_prefix = [None, None, None]
sorted_headers = []

for header in pre_sorted_headers[1:]:
    if header['fullName'].lower() in x_pos_terms:
        assert (header_prefix[0] is None)
        header_prefix[0] = header
    elif header['fullName'].lower() in y_pos_terms:
        assert (header_prefix[1] is None)
        header_prefix[1] = header
    elif header['fullName'].lower() in phenotype_terms:
        assert (header_prefix[2] is None)
        header_prefix[2] = header
    else:
        sorted_headers.append(header)


all_headers = header_prefix + sorted_headers
all_headers: dict

headerlist = []
i = 1
for header in all_headers:
    headerlist.extend(
        chunkify(header['fullName'], i)
    )
    i += 1

payload = {
  'originalData': combinedChannelData,
  'headerList': headerlist,
  'normalizeCsv': False,
  'idField': idfield
}

# execute request
start = time.time()
resp = requests.post(f"{url}/save_config", json=payload)
duration = time.time() - start
print(f"processing duration:\n\t{duration} sec ({duration/60} mins)")

time.sleep(5)

# extract relevant config.json bit
with open(config_path, 'r') as innie:
    config = json.load(innie)


# dataset_config = config[dataset_name]
# dataset_config['linkedDatasets'] = []
# segmentation_filename = label_file.split('/')[-1]

# dataset_dir = f"{data_path}/{dataset_name}"
# if os.path.isfile(f"{dataset_dir}/{segmentation_filename}"):
#     print('its there!')
#     seg_path = f"{dataset_dir}/{segmentation_filename}"
# else:
#     print('it aint there!')
#     seg_path = label_file
# dataset_config['segmentation'] = seg_path


# # write it to file
# with open(f"{dataset_dir}/subconfig.json", "w") as outtie:
#     json.dump({dataset_name: dataset_config}, outtie)

# zip it n ship it
with ZipFile(f'{out_file}', 'w', ZIP_BZIP2) as zipped:
    path_len = len(data_path) + 1
    for base, dirs, files in os.walk(data_path):
        for file in files:
            file_name = os.path.join(base, file)
            zipped.write(file_name, file_name[path_len:])
