import os
from typing import Optional, List
import logging

from sturgeon.utils import validate_model_file, get_model_path
from sturgeon.prediction import predict_samples

def predict(
    input_path: List[str],
    model_files: List[str],
    output_path: str,
    plot_results: Optional[bool] = False,
):
    
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
        raise ValueError(err_msg)

    if not os.path.isdir(output_path):
        err_msg = '''
        --output-path must be a directory, given: {}
        '''.format(output_path)
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
        predict_samples(
            bed_files = bed_files,
            model_file = model,
            output_dir = output_path,
            plot_results = plot_results,
        )




    





    


    