import os
from typing import Optional, List
import logging
import zipfile
import json
from pathlib import Path
from copy import deepcopy

from sturgeon.utils import validate_model_file, get_model_path
from sturgeon.prediction import load_model, predict_sample
from sturgeon.plot import plot_prediction

def predict(
    input_path: List[str],
    model_files: List[str],
    output_path: str,
    plot_results: Optional[bool] = False,
):

    logging.info("Sturgeon start up")
    logging.info("Prediction program")
    
    bed_files = list()
    if os.path.isfile(input_path):
        bed_files.append(input_path)
    elif os.path.isdir(input_path):
        for f in os.listdir(input_path):
            if not f.endswith('.bed'):
                continue
            bed_files.append(os.path.join(input_path, f)) 
    else:
        err_msg = '''
        --input-path must be a directory or file, given: {}
        '''.format(input_path)
        logging.error(err_msg)
        raise ValueError(err_msg)

    if not os.path.isdir(output_path):
        os.makedirs(output_path)
        logging.info('''
        Output path does not exist, creating: {}
        '''.format(output_path))
    else:
        if len(os.listdir(output_path)):    
            err_msg = '''
            --output-path {} contains files. Delete them or provide another path
            '''.format(output_path)
            logging.error(err_msg)
            raise ValueError(err_msg)

    logging.info("Found a total of {} bed files".format(len(bed_files)))
    logging.info("Found a total of {} model files".format(len(model_files)))
    logging.info("Results will be saved in: {}".format(output_path))

    for model in model_files:

        model = get_model_path(model)

        logging.info("Validating model: {}".format(model))
        valid_model = validate_model_file(model)

        if valid_model:
            logging.info("Successful model validation")
        else:
            logging.error("Model did no pass validation, it will be skipped")
            continue

        logging.info("Starting prediction")

        inference_session, probes_df, decoding_dict, temperatures, weight_matrix, merge_dict = load_model(model)                

        for bed_file in bed_files:

            bed_name = Path(bed_file).stem

            prediction_df = predict_sample(
                inference_session = inference_session,
                bed_file = bed_file,
                decoding_dict = deepcopy(decoding_dict),
                probes_df = probes_df,
                weight_matrix = weight_matrix,
                temperatures = temperatures,
                merge_dict = merge_dict,
            )

            output_csv = os.path.join(
                output_path, 
                bed_name + '_{}.csv'.format(Path(model).stem)
            )
            logging.info('Saving results to: {}'.format(output_csv))
            prediction_df.to_csv(
                output_csv, 
                header = True, 
                index = False,
            )

            if plot_results:
                output_pdf = os.path.join(
                    output_path, 
                    bed_name + '_{}.pdf'.format(Path(model).stem)
                )
                logging.info('Plotting results to: {}'.format(output_pdf))

                with zipfile.ZipFile(model, 'r') as zipf:
                    logging.debug("Loading colors dict")
                    try:
                        color_dict = json.load(zipf.open('colors.json'))
                    except FileNotFoundError:
                        logging.debug("No colors dict in zip file")
                        color_dict = None
                
                plot_prediction(
                    prediction_df = prediction_df,
                    color_dict = color_dict,
                    output_file = output_pdf
                )
            else:
                logging.info('Skipping plotting results')

    