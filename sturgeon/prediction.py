from pathlib import Path
import json
import zipfile
from typing import List
import logging
from copy import deepcopy

import numpy as np
import pandas as pd
import onnxruntime

from sturgeon.utils import load_bed_file, softmax, merge_predictions
from sturgeon.constants import METHYL_VALUE, UNMETHYL_VALUE, NOMEASURE_VALUE

def bed_to_numpy(
    bed_df: pd.DataFrame, 
    probes_df: pd.DataFrame,
) -> np.ndarray:
    """Convert the bed dataframe into a numpy array ready for prediction

    Args:
        bed_file (pd.DataFrame): loaded dataframe with methylation status
        probes_df (pd.DataFrame): dataframe with the probes information

    Returns a `OrtValue` with shape [number_probes]
    """
    methyl_col = 'methylation_call'

    bed_df.loc[bed_df[methyl_col] == 0, methyl_col] = UNMETHYL_VALUE
    bed_df = bed_df.set_index('probe_id')
    bed_df = bed_df.reindex(index=probes_df['ID_REF'])
    bed_df = bed_df.reset_index()

    x = np.transpose(np.array(bed_df[methyl_col]))
    x = x.astype(np.float32)
    x[np.isnan(x)] = NOMEASURE_VALUE

    t = 'Total amount of measurable probes:'
    logging.info('{0:45s} {1:6d}'.format(t, len(x)))
    
    t = 'Number of not measured probes:'
    i = np.sum(x == NOMEASURE_VALUE)
    p = (i/len(x))*100
    logging.info('{0:45s} {1:6d} ({2:3.2f}%)'.format(t, i, p))

    t = 'Number of measured probes:'
    i = np.sum(x != NOMEASURE_VALUE)
    p = (i/len(x))*100
    logging.info('{0:45s} {1:6d} ({2:3.2f}%)'.format(t, i, p))

    t = 'Number of measured methylated probes:'
    m = np.sum(x == METHYL_VALUE)
    p = (m/i)*100
    logging.info('{0:45s} {1:6d} ({2:3.2f}%)'.format(t, m, p))

    t = 'Number of measured non-methylated probes:'
    m = np.sum(x == UNMETHYL_VALUE)
    p = (m/i)*100
    logging.info('{0:45s} {1:6d} ({2:3.2f}%)'.format(t, m, p))

    x = np.expand_dims(x, 0)
    
    return x

def model_forward(
    x: np.ndarray, 
    inference_session: onnxruntime.InferenceSession
) -> np.ndarray:


    # compute ONNX Runtime output prediction
    x = onnxruntime.OrtValue.ortvalue_from_numpy(x)
    ort_inputs = {inference_session.get_inputs()[0].name: x}
    ort_outs = inference_session.run(
        ['predictions'], 
        ort_inputs,
    )[0]

    uncalibrated_scores = ort_outs[0]

    return uncalibrated_scores

def load_model(model_file):

    with zipfile.ZipFile(model_file, 'r') as zipf:

        try:
            merge_dict = json.load(zipf.open('merge.json'))
        except KeyError:
            merge_dict = None

        logging.debug("Loading probes information")
        probes_df = pd.read_csv(
            zipf.open('probes.csv'), 
            header = 0, 
            index_col = None,
        )

        logging.debug("Loading the decoding dict")
        decoding_dict = json.load(zipf.open('decoding.json'))

        logging.debug("Loading calibration matrix")
        try:
            temperatures = np.load(zipf.open('calibration.npy'))
            temperatures = temperatures.flatten()
        except KeyError:
            logging.debug("No calibration matrix in zip file")
            temperatures = None

        logging.info("Starting inference session")
        so = onnxruntime.SessionOptions()
        so.inter_op_num_threads = 1
        so.intra_op_num_threads = 1

        inference_session = onnxruntime.InferenceSession(
            zipf.read('model.onnx'), 
            providers = ['CPUExecutionProvider'],
            sess_options = so,
        )

    return inference_session, probes_df, decoding_dict, temperatures, merge_dict


def predict_sample(
    inference_session: str, 
    bed_file: str,
    decoding_dict: dict,
    probes_df: pd.DataFrame,
    temperatures: np.ndarray,
    merge_dict: dict, 
):

    requires_merging = False
    requires_calibration = False

    if merge_dict is not None:
        requires_merging = True

    if temperatures is not None:
        requires_calibration = True

    logging.info("Loading bed file: {}".format(bed_file))
    x = bed_to_numpy(
        bed_df = load_bed_file(bed_file), 
        probes_df = probes_df
    )
    scores = model_forward(x, inference_session)

    calibrated_scores = np.empty_like(scores)
    for i in range(scores.shape[0]):
        if requires_calibration:
            calibrated_scores[i, :]  = np.exp(softmax(scores[i, :]/temperatures[i]))
        else:
            calibrated_scores[i, :]  = np.exp(softmax(scores[i, :]))

    calibrated_df = dict()
    for k, v in decoding_dict.items():
        calibrated_df[v] = calibrated_scores[:, int(k)]
    calibrated_df = pd.DataFrame(calibrated_df)

    if requires_merging:
        calibrated_df, decoding_dict = merge_predictions(calibrated_df, decoding_dict, merge_dict)
        
    n = np.sum(x != NOMEASURE_VALUE)

    avg_scores = {
        'number_probes': [n],
    }

    final_scores = list()
    arr = np.array(calibrated_df[calibrated_df.columns])
    best_m = np.where(arr == np.max(arr))[0]
    if len(best_m) > 0:
        best_m = best_m[0]
    for colname in calibrated_df.columns:
        score = np.array(calibrated_df[colname])[best_m].item()
        avg_scores[colname] = [score]
        final_scores.append(score)
    final_scores = np.array(final_scores)
    
    prediction_df = pd.DataFrame(avg_scores, index = [0])

    top3 = final_scores.argsort(-1)[::-1][:3]
    for i, t in enumerate(top3):
        logging.info('Top {0}: {1:30s} ({2:4.3f})'.format(
            i+1, calibrated_df.columns[t], final_scores[t]
        ))

    return prediction_df

            