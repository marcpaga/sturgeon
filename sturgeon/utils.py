import zipfile
import logging
import os
from pathlib import Path
import platform

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

def read_probes_file(probes_file: str):

    samtools_bed_file = False
    with open(probes_file, 'r') as handle:
        for line in handle:
            break
        if line.startswith('#'):
            samtools_bed_file = True
            column_header = line.strip('#').strip('\n').split('\t')

    if samtools_bed_file:
        probes_df = pd.read_csv(
            probes_file,
            header = None,
            index_col = None,
            sep = '\t',
            comment='#',
            names=column_header,
        )
    # normal probes file
    else:
        probes_df = pd.read_csv(
            probes_file,
            header = 0,
            index_col = None,
            sep = ' ',
        )

    return probes_df

def validate_megalodon_file(megalodon_file):
    """Very simple check of the first line in the file to check
    for correct column names.

    Args:
        megalodon_file (str): path to the file to be validated

    Returns:
        (bool, str): True, None if it passes the validation; False, error message
        if it failes.
    """

    mandatory_columns = [
        "read_id", 
        "chrm",
        "strand",
        "pos",
        "mod_log_prob",
        "can_log_prob",
        "mod_base",
    ]

    with open(megalodon_file, 'r') as f:
        for line in f:
            break

    column_names = line.strip('\n').split('\t')
    for mc in mandatory_columns:
        if mc not in column_names:
            err_msg = "{} column missing in megalodon file".format(mc)
            return False, err_msg

    return True, None


def validate_modkit_file(modkit_file):
    """Very simple check of the first line in the file to check
    for correct column names.

    Args:
        modkit_file (str): path to the file to be validated

    Returns:
        (bool, str): True, None if it passes the validation; False, error message
        if it failes.
    """

    mandatory_columns = [
        "read_id", 
        "chrom",
        "ref_position",
        "mod_qual",
        "mod_code",
        "canonical_base",
        "modified_primary_base",
    ]

    with open(modkit_file, 'r') as f:
        for line in f:
            break

    column_names = line.strip('\n').split('\t')
    for mc in mandatory_columns:
        if mc not in column_names:
            err_msg = "{} column missing in modkit file".format(mc)
            return False, err_msg

    return True, None

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

def creation_date(path_to_file):
    """
    Try to get the date that a file was created, falling back to when it was
    last modified if that isn't possible.
    See http://stackoverflow.com/a/39501288/1709587 for explanation.
    """
    if platform.system() == 'Windows':
        return os.path.getctime(path_to_file)
    else:
        stat = os.stat(path_to_file)
        try:
            return stat.st_birthtime
        except AttributeError:
            # We're probably on Linux. No easy way to get creation dates here,
            # so we'll settle for when its content was last modified.
            return stat.st_mtime

def softmax(x):
    c = x.max()
    logsumexp = np.log(np.exp(x - c).sum())
    return x - c - logsumexp

def merge_predictions(prediction_df, decoding_dict, merge_dict):

    old_class_columns = list(decoding_dict.values())
    non_class_columns = np.array(prediction_df.columns[~np.isin(prediction_df.columns, old_class_columns)])
   
    rev_merge_dict = dict()
    for k, v in merge_dict.items():
        for c in v:
            rev_merge_dict[c] = k
    for v in decoding_dict.values():
        try:
            rev_merge_dict[v]
        except KeyError:
            rev_merge_dict[v] = v
    
    for k, v in merge_dict.items():
        prediction_df[k] = 0
        for c in v:
            prediction_df[k] += prediction_df[c]
            prediction_df = prediction_df.drop(c, axis=1)

    new_class_columns = np.array(prediction_df.columns[~np.isin(prediction_df.columns, non_class_columns)])

    new_class_columns = np.sort(new_class_columns)

    new_decoding_dict = dict()
    for i, c in enumerate(new_class_columns):
        new_decoding_dict[i] = c

    final_column_order = non_class_columns.tolist() + new_class_columns.tolist()

    return prediction_df[final_column_order], new_decoding_dict