from typing import Optional
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib.backends.backend_pdf # required for compilation

def plot_prediction(
    prediction_df: pd.DataFrame,
    output_file: str,
    color_dict: Optional[dict] = None,
):

    SMALL_SIZE = 8
    MEDIUM_SIZE = 12
    BIGGER_SIZE = 18

    plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

    non_label_columns = ['number_probes', 'timestamp']
    label_columns = []
    for c in prediction_df.columns:
        if c not in non_label_columns:
            label_columns.append(c)

    title = Path(output_file).stem + '\n' + 'Measured probes: {}'.format(prediction_df['number_probes'].item())

    
    plt.figure(figsize = (len(label_columns)//4, 7))
    for i, c in enumerate(label_columns):
        v = prediction_df[c].item()

        if color_dict is not None:
            color = color_dict[c]
        else:
            color = 'grey'

        plt.axhline(0.8, color = 'grey', linestyle = '--', zorder = 1)
        plt.axhline(0.95, color = 'grey', linestyle = '--', zorder = 1)
        plt.bar(
            x = i, 
            height = v, 
            edgecolor = 'black', 
            color = color, 
            zorder = 2
        )
        
    plt.xlim(-1, len(label_columns))
    plt.ylim(0, 1.05)
    plt.xticks(np.arange(len(label_columns)), label_columns, rotation = 90)
    plt.ylabel('Score')
    plt.title(title)
    plt.tight_layout()
    plt.savefig(output_file, bbox_inches="tight", dpi = 300)
    plt.close()


def plot_prediction_over_time(
    prediction_df: pd.DataFrame,
    output_file: str,
    color_dict: Optional[dict] = None,
    average_score_threshold: Optional[float] = 0.1,
):

    SMALL_SIZE = 12
    MEDIUM_SIZE = 16
    BIGGER_SIZE = 18

    plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=BIGGER_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

    non_label_columns = ['number_probes', 'timestamp']
    label_columns = []
    for c in prediction_df.columns:
        if c not in non_label_columns:
            label_columns.append(c)

    title = Path(output_file).stem

    classes_to_plot = list()
    for lc in label_columns:
        if np.mean(prediction_df[lc]) >= average_score_threshold:
            classes_to_plot.append(lc)

    plt.figure(figsize = (20, 8))
    for _, c in enumerate(classes_to_plot):
        y = np.array(prediction_df[c])
        x = np.arange(len(y))

        if color_dict is not None:
            color = color_dict[c]
        else:
            color = None

        plt.axhline(0.8, color = 'grey', linestyle = '--', zorder = 1)
        plt.axhline(0.95, color = 'grey', linestyle = '--', zorder = 1)
        plt.plot(x, y, color = 'black', zorder = 2)
        plt.scatter(
            x = x, 
            y = y, 
            c = color, 
            s = 60, 
            zorder = 3, 
            label = c, 
            edgecolor='black'
        )
        
        
    plt.ylim(0, 1.05)
    plt.xticks(
        np.arange(len(prediction_df['number_probes'])), 
        prediction_df['number_probes'], rotation = 45
    )
    plt.xlabel('Number of measured probes')
    plt.ylabel('Score')
    plt.title(title)
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.tight_layout()
    plt.savefig(output_file, bbox_inches="tight")
    plt.close()

    