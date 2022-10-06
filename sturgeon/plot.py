from typing import Optional

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

def plot_prediction(
    prediction_df: pd.DataFrame,
    output_file: str,
    color_dict: Optional[dict] = None,
):
    raise NotImplementedError


def plot_prediction_over_time(
    prediction_df: pd.DataFrame,
    output_file: str,
    color_dict: Optional[dict] = None,
):
    raise NotImplementedError