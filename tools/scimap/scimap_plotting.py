import argparse
import json
import os
import warnings

import matplotlib.pyplot as plt
import numpy as np
import scimap as sm
import seaborn as sns
from anndata import read_h5ad
from matplotlib import colormaps

sns.set(color_codes=True)


def main(inputs, anndata, output):
    """
    Parameter
    ---------
    inputs : str
        File path to galaxy tool parameter.
    anndata : str
        File path to anndata containing phenotyping info.
    output : str
        File path to output.
    """
    warnings.simplefilter('ignore')

    with open(inputs, 'r') as param_handler:
        params = json.load(param_handler)

    adata = read_h5ad(anndata)

    tool = params['analyses']['selected_tool']
    options = params['analyses']['options']

    if tool == 'stacked_barplot':

        # parse list text arguments
        for o in options.copy():
            opt = options.pop(o)
            if o == 'matplotlib_cmap':
                matplotlib_cmap = opt
            elif opt != "":
                options[o] = [x.strip() for x in opt.split(',')]

        # add base args into options dict to pass to tool
        options['x_axis'] = params['analyses']['x_axis']
        options['y_axis'] = params['analyses']['y_axis']
        options['method'] = params['analyses']['method']

        options['return_data'] = True

        df = sm.pl.stacked_barplot(adata, **options)

        # Pick cmap to use
        if matplotlib_cmap not in list(colormaps):
            print('Using default color options')
            print('For different y-axis category colors,')
            print('enter one of the following strings for the colormap param:')
            print(list(colormaps))
            num_phenotypes = len(df.columns) - 1
            if num_phenotypes <= 9:
                matplotlib_cmap = "Set1"
            elif num_phenotypes > 9 and num_phenotypes <= 20:
                matplotlib_cmap = plt.cm.tab20
            else:
                matplotlib_cmap = plt.cm.gist_ncar

        # Plotting
        sns.set_theme(style="white")
        ax = df.plot.bar(stacked=True, cmap=matplotlib_cmap)
        fig = ax.get_figure()
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(
            reversed(handles),
            reversed(labels),
            bbox_to_anchor=(1, 1),
            loc='upper left'
        )
        plt.tight_layout()

        # # save and close
        fig.savefig('out.png', dpi=300)
        plt.close(fig)

    if tool == 'voronoi':

        plt.style.use('fast')

        tool_func = getattr(sm.pl, tool)

        # x_lim/y_lim need to be parsed from comma-sep str to integer tuples
        for lim in ['x_lim', 'y_lim']:
            opt_list = options.pop(lim)
            if opt_list:
                options[lim] = [int(x.strip()) for x in opt_list.split(',')]

        # parse list text arguments
        for cat in ['overlay_points_categories', 'overlay_drop_categories']:
            opt_list = options.pop(cat)
            if opt_list:
                options[cat] = [x.strip() for x in opt_list.split(',')]

        # add base args into options dict to pass to tool
        options['color_by'] = params['analyses']['color_by']
        options['x_coordinate'] = params['analyses']['x_coordinate']
        options['y_coordinate'] = params['analyses']['y_coordinate']

        # fill any missing params with None as tool expects
        for k, v in options.items():
            if v == '':
                options[k] = None

        options['saveDir'] = os.getcwd()
        options['fileName'] = 'out.png'

        if options['size_max'] is None:
            options['size_max'] = np.inf

        # call the tool and unpack all options
        tool_func(adata, **options)


if __name__ == '__main__':
    aparser = argparse.ArgumentParser()
    aparser.add_argument("-i", "--inputs", dest="inputs", required=True)
    aparser.add_argument("-a", "--anndata", dest="anndata", required=True)
    aparser.add_argument("-e", "--output", dest="output", required=True)

    args = aparser.parse_args()

    main(args.inputs, args.anndata, args.output)
