import os
import sys
import re

import numpy as np
import pandas as pd
import anndata as ad

import plotly.graph_objects as go
import matplotlib.pyplot as plt


def plot_cycle_dropout(adata, pattern = 'DAPI'):

    """
    Produces QC plots to show degredation between imaging cycles
    - scatter plot of DAPI/DNA channels between cycle 1 and every subsequent cycle
    - line plot of R-squared values between cycle 1 and every subsequent cycle DAPI/DNA
    - violin plot of DAPI/DNA raw intensity distribution for every cycle
    """

    df = pd.DataFrame(adata.raw.X, columns = adata.raw.var_names)
    QC_cols = [x for x in df.columns if re.search(pattern, x)]

    # iterate thru DAPI/DNA cols, using the first as the reference channel (assuming segmentation was run on this channel)
    scores = []
    for i,marker in enumerate(QC_cols):
        if i == 0:
            ref = marker
        else:

            # line of best fit
            theta = np.polyfit(df[ref], df[marker], 1)
            y_line = theta[1] + theta[0] * df[ref]

            # R-squared value
            corr_matrix = np.corrcoef(df[ref], df[marker])
            corr = corr_matrix[0,1]
            R_sq = round(corr**2,2)
            scores.append(R_sq)
            R_sq = str(R_sq)
            annotation = r"$R^{2} = $" + R_sq

            # scatter plot w/line of best fit and R^2 annotation between current DAPI/DNA and reference DAPI/DNA
            fig,ax = plt.subplots()
            ax.scatter(x = df[ref], y = df[marker], s = 0.5)
            ax.plot(df[ref], y_line, ls = '--', color = 'grey', alpha = 0.5)
            ax.annotate(annotation, (0,df[marker].max()-2000))
            ax.set_xlabel(ref)
            ax.set_ylabel(marker)
            ax.set_aspect('equal')
            fig.savefig(f"./figures/{ref}_{marker}_scatter.png", dpi = 500)
        
    # R^2 values were saved in scores list, plot R^2 across cycles
    cycles = list(range(2,len(QC_cols)+1))
    fig,ax = plt.subplots()
    ax.plot(cycles,scores)
    ax.set_xlabel("Cycle number")
    ax.set_ylabel(r"$R^{2}$")
    ax.set_title("Goodness of fit between cycle 1 and every subsequent cycle")
    fig.savefig("./figures/cycle_Rsquared.png", dpi = 500)

    # finally, a violin plot showing the distribution of DAPI/DNA intensities in each cycle
    fig = go.Figure()
    for col in QC_cols:
        fig.add_trace(go.Violin(
            y= df[col],
            name=col,
        ))

    fig.update_layout(title_text="DAPI intensity distribution by cycle", 
    xaxis_title = "Cycle", yaxis_title = "Intensity", title_x = 0.5)
    fig.write_image("./figures/cycle_violin.png", scale = 5)

    return