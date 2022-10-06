import os
from typing import Optional, List

from sturgeon.utils import validate_model_file
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

    for model in model_files:
        validate_model_file(model)

        predict_samples(
            bed_files = bed_files,
            model_file = model,
            output_dir = output_path,
            plot_results = plot_results,
        )




    





    


    