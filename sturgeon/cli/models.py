import os
import shutil
from pathlib import Path
import logging

from sturgeon.utils import get_available_models, validate_model_file


def actions_models(action, model_files):

    if action == 'list':
        list_models()
    elif action == 'add':
        add_models(model_files)
    elif action == 'delete':
        del_models(model_files)

def list_models():

    available_models = get_available_models(print_str=True)

    msg = '''
    The following models are available: {}
    '''.format(available_models)

    print(msg)

    return None


def add_models(model_files):

    for model in model_files:

        model_name = Path(model).stem
        output_file = os.path.join(
            os.path.dirname(__file__), 
            '../include/models', 
            model_name + '.zip'
        )

        if os.path.isfile(output_file):
            logging.warning('''
            Model file with the same name already exists, delete that file first
            or use a different name. Skipping: {}
            '''.format(model)
            )
            continue

        if not model.endswith('.zip'):
            logging.warning("Given model file is not a zip file, skipping: {}".format(model))
            continue

        if not validate_model_file(model):
            logging.error("Given model file is not a valid model file, skipping: {}".format(model))
            continue

        logging.info("Given model file is valid, copying to model path")
        shutil.copyfile(
            model,
            output_file,
        )
        logging.info("Model added: {}".format(output_file))


    list_models()



def del_models(model_files):

    for model in model_files:

        if not os.path.isfile(model):
            logging.warning('''
            Model file nout found. Skipping: {}
            '''.format(model)
            )
            continue

        logging.info("Deleting model file: {}".format(model))
        os.remove(model)


    list_models()