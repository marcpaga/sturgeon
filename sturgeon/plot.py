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

    non_label_columns = ['number_probes']
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

        plt.bar(x = i, height = v, edgecolor = 'black', color = color)
        
    plt.xlim(-1, len(label_columns))
    plt.ylim(0, 1.05)
    plt.xticks(np.arange(len(label_columns)), label_columns, rotation = 90)
    plt.ylabel('Score')
    plt.title(title)
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()



def plot_prediction_over_time(
    prediction_df: pd.DataFrame,
    output_file: str,
    color_dict: Optional[dict] = None,
):
    raise NotImplementedError