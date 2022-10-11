import os
from pathlib import Path
import json
import zipfile
from typing import Optional, List

import numpy as np
import pandas as pd
import onnxruntime

from sturgeon.utils import load_bed_file
from sturgeon.calibration import HistogramCalibration
from sturgeon.plot import plot_prediction
from sturgeon.constants import UNMETHYL_VALUE, NOMEASURE_VALUE

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
    x = np.expand_dims(x, 0)
    
    return x

def predict_sample(
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

    uncalibrated_scores = np.exp(ort_outs[0])

    return uncalibrated_scores

def predict_samples(
    model_file: str, 
    bed_files: List[str], 
    output_dir: str,
    plot_results: Optional[bool] = False,
):

    model_name = Path(model_file).stem
    with zipfile.ZipFile(model_file, 'r') as zipf:

        decoding_dict = json.load(zipf.open('decoding.json'))
        probes_df = pd.read_csv(
            zipf.open('probes.csv'), 
            header = 0, 
            index_col = None,
        )

        try:
            calibrators = list()
            calibration_matrix = np.load(zipf.open('calibration.npy'))
            for i in range(calibration_matrix.shape[-1]):
                calibrator = HistogramCalibration(
                    num_classes = len(decoding_dict)
                )
                calibrator.load_matrix(calibration_matrix[:, :, i])
                calibrators.append(calibrator)
        except FileNotFoundError:
            calibrators = None

        try:
            color_dict = json.load(zipf.open('colors.json'))
        except FileNotFoundError:
            color_dict = None

        try:
            weight_matrix = np.load(zipf.open('weight_scores.npz'))
            mean_probes_per_timepoint = weight_matrix['avgsites']
            accuracy_per_timepoint_per_model = weight_matrix['performance']
            calculate_weighted_mean = True
        except FileNotFoundError:
            calculate_weighted_mean = False

        so = onnxruntime.SessionOptions()
        so.inter_op_num_threads = 1
        so.intra_op_num_threads = 1

        inference_session = onnxruntime.InferenceSession(
            zipf.read('model.onnx'), 
            providers = ['CPUExecutionProvider'],
            sess_options = so,
        )

        for bed_file in bed_files:

            file_name = Path(bed_file).stem
            output_stem = os.path.join(
                output_dir, 
                file_name+'_{}'.format(model_name)
            )
            
            x = bed_to_numpy(
                bed_df = load_bed_file(bed_file), 
                probes_df = probes_df
            )
            scores = predict_sample(x, inference_session)

            n = np.sum(x != NOMEASURE_VALUE)
            if calculate_weighted_mean:
                # take the weighted average of all models
                calculated_weights = np.ones(scores.shape, dtype=float)

                for m in range(calculated_weights.shape[0]):

                    weights = accuracy_per_timepoint_per_model[m]
                    n_probes = mean_probes_per_timepoint[m]
                    t = n_probes.searchsorted(n)
                    t = int(t)
                    if t == weights.shape[0]:
                        calculated_weights[m, :] = weights[t-1]
                    elif t == 0:
                        calculated_weights[m, :] = weights[t]
                    else:
                        weights = weights[t-1:t+1]
                        x = [n_probes[t-1], n_probes[t]]
                        for i in range(weights.shape[1]):
                            y = weights[:, i]
                            calculated_weights[m, i] = np.interp(n, x, y)

                final_scores = np.zeros(scores.shape[1])
                for i in range(scores.shape[1]):
                    final_scores[i] = np.average(
                        a = scores[:, i], 
                        weights = calculated_weights[:, i]
                    )

            else:
                final_scores = scores.mean(0)

            
            prediction_df = {'number_probes': n}
            for i in range(final_scores.shape[0]):
                prediction_df[decoding_dict[str(i)]] = final_scores[i]
            prediction_df = pd.DataFrame(prediction_df, index = [0])
            prediction_df.to_csv(
                output_stem + '.csv', 
                header = True, 
                index = False,
            )

            if plot_results:
                plot_prediction(
                    prediction_df = prediction_df,
                    color_dict = color_dict,
                    output_file = output_stem + '.pdf'
                )

            