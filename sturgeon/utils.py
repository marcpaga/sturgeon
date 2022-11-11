import zipfile
import logging
import os
from pathlib import Path

import pandas as pd
import numpy as np

from sturgeon.constants import (
    UNMETHYL_VALUE,
    NOMEASURE_VALUE,
)

def validate_model_file(zip_file: str):
    """Validate the contents of a zip file

    Args:
        zip_file (str): path to the model zip file to be validated

    Raises an error if critical components are missing. Produces warnings if
    non-critical components are missing.
    """
    files_in_zip = zipfile.ZipFile(zip_file).namelist()

    mandatory_files = [
        'model.onnx',
        'decoding.json',
        'probes.csv'
    ]

    for mf in mandatory_files:
        if mf not in files_in_zip:
            err_msg = 'Mandatory file: {mf}, not found in {zf}'.format(
                mf = mf,
                zf = zip_file,
            )
            logging.error(err_msg)
            return False
            

    non_mandatory_files = {
        'calibration.npy': 'score calibration will not be possible',
        'colors.json': 'default colors will be used',
        'weight_scores.npz': 'weighted average score will not be possible'
    }

    for nmf, err in non_mandatory_files.items():
        if nmf not in files_in_zip:
            wrn_msg = 'Mandatory file: {mf}, not found in {zf}: {err}'.format(
                mf = mf,
                zf = zip_file,
                err = err,
            )
            logging.warning(wrn_msg)

    return True

def load_bed_file(bed_file: str):
    """Read the contents of a bed file

    Args:
        bed_file (str): path to methylation bed file to be read

    Returns a pandas.DataFrame with the contents of the bed file
    """

    bed_df =  pd.read_csv(
        bed_file, 
        header = 0, 
        index_col = None, 
        delim_whitespace=True
    )
    return bed_df

def validate_bed_file(bed_df: pd.DataFrame, probes_df: pd.DataFrame):
    """Validate the contents of a bed file

    Args:
        bed_file (pd.DataFrame): loaded dataframe with methylation status
        probes_df (pd.DataFrame): dataframe with the probes information

    Raises an error if critical components are missing. Produces warnings if
    non-critical components are missing.
    """

    mandatory_columns = ["methylation_call", "probe_id"]
    for mc in mandatory_columns:
        err_msg = "{} column missing in bed file".format(mc)
        assert mc in bed_df.columns, err_msg

    invalid_vals_idx = np.where(~bed_df['methylation_call'].isin([0, 1]))[0]
    if len(invalid_vals_idx) > 0:
        err_msg = """
        Valid methylation values in 'methylation_call' column are 0 or 1.
        Found the following invalid values:\n
        """
        for i in invalid_vals_idx:
            v = bed_df['methylation_call'].tolist()[i]
            err_msg += "Line {i}: {v}\n".format(i = i+1, v = v)
            raise ValueError(err_msg)

    return True

def get_available_models(print_str = False):

    model_dir = os.path.join(os.path.dirname(__file__), 'include/models')

    available_models = list()
    for model in os.listdir(model_dir):
        if not model.endswith('.zip'):
            continue
        available_models.append(Path(model).stem)
    
    if print_str:
        available_models = "\n"+"\n".join(available_models)

    return available_models

def get_model_path(model_name):

    if os.path.isfile(model_name):
        return model_name

    model_dir = os.path.join(os.path.dirname(__file__), 'include/models')
    if model_name not in get_available_models():
        err_msg = """
        The following model could not be found in {}
        """.format(model_dir)
        raise ValueError(err_msg)

    return os.path.join(model_dir, model_name + '.zip')

