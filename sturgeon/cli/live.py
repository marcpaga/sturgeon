import os
import logging
from typing import List
import time
from pathlib import Path
from copy import deepcopy
import shutil
import zipfile
import json

import numpy as np
import pandas as pd
import pysam

from sturgeon.callmapping import (
    bam_to_calls, 
    merge_probes_methyl_calls,
    probes_methyl_calls_to_bed,
    mega_file_to_bed,
)
from sturgeon.utils import (
    validate_model_file, 
    get_model_path, 
    creation_date, 
    validate_megalodon_file,
)

from sturgeon.prediction import predict_sample, load_model
from sturgeon.plot import plot_prediction, plot_prediction_over_time
from sturgeon.utils import read_probes_file


def live(
    input_path: str,
    output_path: str,
    model_files: List[str],
    source: str,
    probes_file: str,
    reference_genome: str,
    margin: int,
    neg_threshold: float,
    pos_threshold: float,
    plot_results: bool,
    cooldown: int,
):
    """
    """

    logging.info("Sturgeon start up")
    logging.info("Live prediction program")

    logging.info('Watching the following folder: {}'.format(input_path))

    if not os.path.exists(output_path):
        os.makedirs(output_path)

    if reference_genome is not None:
        probes_file = os.path.join(
            os.path.dirname(__file__), 
            '../include/static', 
            'probes_{}.bed'.format(reference_genome)
        )

    if source == 'guppy':

        probes_df = read_probes_file(probes_file)
        
        live_guppy(
            input_path = input_path,
            output_path = output_path,
            model_files = model_files,
            probes_df = probes_df,
            margin = margin,
            neg_threshold = neg_threshold,
            pos_threshold = pos_threshold,
            plot_results = plot_results,
            cooldown = cooldown,
        )

    elif source == 'megalodon':

        live_megalodon(
            input_path = input_path,
            output_path = output_path,
            model_files = model_files,
            probes_file = probes_file,
            margin = margin,
            neg_threshold = neg_threshold,
            pos_threshold = pos_threshold,
            plot_results = plot_results,
            cooldown = cooldown,
        )


def live_guppy(
    input_path: str,
    output_path: str,
    model_files: List[str],
    probes_df: pd.DataFrame,
    margin: int,
    neg_threshold: float,
    pos_threshold: float,
    plot_results: bool,
    cooldown: int,
):

    # keep track of processed bam files
    bam_files = dict()
    logging.info('Starting live prediction from bam files')

    while True:

        logging.info(
            '''
            Looking for new bam files, so far found {}
            '''.format(len(bam_files))
        )

        # get all available bam files and sort them by timestamp so that we
        # process them in an orderly manner
        available_bam_files = list()
        available_bam_timestamps = list()
        for file in os.listdir(input_path):
            if not file.endswith('.bam'):
                continue
            f = os.path.join(input_path, file)
            available_bam_files.append(f)
            available_bam_timestamps.append(creation_date(f))

        available_bam_files = np.array(available_bam_files)
        available_bam_timestamps = np.array(available_bam_timestamps)
        creation_order = np.argsort(available_bam_timestamps)

        available_bam_files = available_bam_files[creation_order]
        available_bam_timestamps = available_bam_timestamps[creation_order]

        for file_path, timestamp in zip(available_bam_files, available_bam_timestamps):
            
            # check if we have processed this bam file already
            file_name = Path(file_path).stem
            try:
                # already processed bam file, skip
                bam_files[file_name]
                continue
            except KeyError:
                logging.info('New bam file found: {}'.format(file_path))
                
            # generate index file if not existing
            bai_file = file_path + '.bai'
            if not os.path.exists(bai_file):
                logging.info(
                    '''
                    Index file not found for bam file: {}
                    '''.format(file_path)
                )
                logging.info(
                    '''
                    Generating index file: {}
                    '''.format(bai_file)
                )
                # sometimes we catch the bam file mid write up, try again
                # in the next iteration
                try:
                    pysam.index(file_path)
                except pysam.utils.SamtoolsError:
                    logging.warning(
                        '''
                        Making index file failed, bam file might not be 
                        complete yet, will try again later
                        '''
                    )
                    continue

                logging.info(
                    '''
                    Generated index file: {}
                    '''.format(bai_file)
                )

            # extract methylation calls from the bam file
            logging.info('Getting methylation calls from file')
            calls_per_probe_file = os.path.join(
                output_path, 
                file_name + '_probes_methyl_calls.txt',
            )
            calls_per_read_file = os.path.join(
                output_path, 
                file_name + '_read_methyl_calls.txt'
            )
            if not os.path.isfile(calls_per_probe_file) or not os.path.isfile(calls_per_read_file):
                
                probes_methyl_df = deepcopy(probes_df)
                calls_per_probe, calls_per_read = bam_to_calls(
                    bam_file = file_path,
                    probes_df = probes_methyl_df,
                    margin = margin,
                    neg_threshold = neg_threshold,
                    pos_threshold = pos_threshold,
                )
                calls_per_probe.to_csv(
                    calls_per_probe_file,
                    header = True, index = False, sep = '\t'
                )
                calls_per_read.to_csv(
                    calls_per_read_file, 
                    header = True, index = False, sep = '\t'
                )

            merged_output_file = os.path.join(
                output_path, 
                'merged_probes_methyl_calls_{}.txt'.format(len(bam_files))
            )

            # if this is the first processed file, there's no need to merge
            if len(bam_files) == 0:
                shutil.copyfile(calls_per_probe_file, merged_output_file)
            # here we merge the current probe file with the previous merge file
            else:
                merge_probes_methyl_calls(
                    [
                        calls_per_probe_file,
                        os.path.join(
                            output_path, 
                            'merged_probes_methyl_calls_{}.txt'.format(len(bam_files)-1)
                        )
                    ], 
                    merged_output_file
                )

            # conver the probe file to bed, so that we can predict
            bed_output_file = os.path.join(
                output_path,
                'merged_probes_methyl_calls_{}.bed'.format(len(bam_files))
            )
            probes_methyl_calls_to_bed(
                merged_output_file,
                bed_output_file
            )

            # add to bam_files as we consider this file processed
            bam_files[file_name] = file_path

            # make a prediction with each model
            for model in model_files:

                model = get_model_path(model)

                logging.info("Validating model: {}".format(model))
                valid_model = validate_model_file(model)

                if valid_model:
                    logging.info("Successful model validation")
                else:
                    logging.error(
                        '''
                        Model did no pass validation, it will be skipped
                        '''
                    )
                    continue

                inference_session, array_probes_df, decoding_dict, temperatures, merge_dict = load_model(model)    

                logging.info("Starting prediction")
                prediction_df = predict_sample(
                    inference_session = inference_session,
                    bed_file = bed_output_file,
                    decoding_dict = deepcopy(decoding_dict),
                    probes_df = array_probes_df,
                    temperatures = temperatures,
                    merge_dict = merge_dict,
                )
                prediction_df['timestamp'] = timestamp

                output_csv = os.path.join(
                    output_path, 
                    'predictions_{}.csv'.format(Path(model).stem)
                )
                prediction_df.to_csv(
                    output_csv,
                    mode = 'a',
                    index = False,
                    header = not os.path.exists(output_csv),
                )

                if plot_results:

                    # plot the last prediction
                    output_pdf = os.path.join(
                        output_path, 
                        'predictions_{}_{}.pdf'.format(
                            len(bam_files)-1,
                            Path(model).stem,
                        )
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

                    output_pdf = os.path.join(
                        output_path, 
                        'predictions_overtime_{}.pdf'.format(
                            Path(model).stem,
                        )
                    )
                    predictions_time = pd.read_csv(
                        output_csv,
                        header = 0,
                        index_col = None,
                    )
                    plot_prediction_over_time(
                        prediction_df = predictions_time,
                        color_dict = color_dict,
                        output_file = output_pdf,
                    )

                else:
                    logging.info('Skipping plotting results')


        logging.info(
            '''
            No new bam files found, sleeping for {} seconds
            '''.format(cooldown)
        )
        time.sleep(cooldown)

def live_megalodon(
    input_path: str,
    output_path: str,
    model_files: List[str],
    probes_file: str,
    margin: int,
    neg_threshold: float,
    pos_threshold: float,
    plot_results: bool,
    cooldown: int,
):

    # keep track of processed bam files
    meg_files = dict()
    logging.info('Starting live prediction from megalodon output files')
    
    while True:

        logging.info(
            '''
            Looking for new megalodon files, so far found {}
            '''.format(len(meg_files))
        )

        # get all available bam files and sort them by timestamp so that we
        # process them in an orderly manner
        available_meg_files = list()
        available_meg_timestamps = list()
        for file in os.listdir(input_path):
            if not file.endswith('.txt'):
                continue
            f = os.path.join(input_path, file)

            success, msg = validate_megalodon_file(f)
                
            if not success:
                logging.info(
                    '''
                    File {}, did not pass validation either not a megalodon file or
                    an invalid megalodon file. Reason: {}.
                    '''.format(f, msg)
                )
                continue

            available_meg_files.append(f)
            available_meg_timestamps.append(creation_date(f))

        available_meg_files = np.array(available_meg_files)
        available_meg_timestamps = np.array(available_meg_timestamps)
        creation_order = np.argsort(available_meg_timestamps)

        available_meg_files = available_meg_files[creation_order]
        available_meg_timestamps = available_meg_timestamps[creation_order]

        for file_path, timestamp in zip(available_meg_files, available_meg_timestamps):
            
            # check if we have processed this bam file already
            file_name = Path(file_path).stem
            try:
                # already processed bam file, skip
                meg_files[file_name]
                continue
            except KeyError:
                logging.info('New megalodon file found: {}'.format(file_path))
                

            # extract methylation calls from the bam file
            logging.info('Getting methylation calls from file')
            calls_per_probe_file = os.path.join(
                output_path, 
                file_name + '_probes_methyl_calls.txt',
            )

            if not os.path.isfile(calls_per_probe_file):
                
                calls_per_probe = mega_file_to_bed(
                    input_file = file_path,
                    probes_file = probes_file,
                    margin = margin,
                    neg_threshold = neg_threshold,
                    pos_threshold = pos_threshold,
                )
                calls_per_probe.to_csv(
                    calls_per_probe_file,
                    header = True, index = False, sep = '\t'
                )


            merged_output_file = os.path.join(
                output_path, 
                'merged_probes_methyl_calls_{}.txt'.format(len(meg_files))
            )

            # if this is the first processed file, there's no need to merge
            if len(meg_files) == 0:
                shutil.copyfile(calls_per_probe_file, merged_output_file)
            # here we merge the current probe file with the previous merge file
            else:
                merge_probes_methyl_calls(
                    [
                        calls_per_probe_file,
                        os.path.join(
                            output_path, 
                            'merged_probes_methyl_calls_{}.txt'.format(len(meg_files)-1)
                        )
                    ], 
                    merged_output_file
                )

            # conver the probe file to bed, so that we can predict
            bed_output_file = os.path.join(
                output_path,
                'merged_probes_methyl_calls_{}.bed'.format(len(meg_files))
            )
            probes_methyl_calls_to_bed(
                merged_output_file,
                bed_output_file
            )

            # add to bam_files as we consider this file processed
            meg_files[file_name] = file_path

            # make a prediction with each model
            for model in model_files:

                model = get_model_path(model)

                logging.info("Validating model: {}".format(model))
                valid_model = validate_model_file(model)

                if valid_model:
                    logging.info("Successful model validation")
                else:
                    logging.error(
                        '''
                        Model did no pass validation, it will be skipped
                        '''
                    )
                    continue

                inference_session, array_probes_df, decoding_dict, temperatures, merge_dict = load_model(model)    

                logging.info("Starting prediction")
                prediction_df = predict_sample(
                    inference_session = inference_session,
                    bed_file = bed_output_file,
                    decoding_dict = deepcopy(decoding_dict),
                    probes_df = array_probes_df,
                    temperatures = temperatures,
                    merge_dict = merge_dict,
                )
                prediction_df['timestamp'] = timestamp

                output_csv = os.path.join(
                    output_path, 
                    'predictions_{}.csv'.format(Path(model).stem)
                )
                prediction_df.to_csv(
                    output_csv,
                    mode = 'a',
                    index = False,
                    header = not os.path.exists(output_csv),
                )

                if plot_results:

                    # plot the last prediction
                    output_pdf = os.path.join(
                        output_path, 
                        'predictions_{}_{}.pdf'.format(
                            len(meg_files)-1,
                            Path(model).stem,
                        )
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

                    output_pdf = os.path.join(
                        output_path, 
                        'predictions_overtime_{}.pdf'.format(
                            Path(model).stem,
                        )
                    )
                    predictions_time = pd.read_csv(
                        output_csv,
                        header = 0,
                        index_col = None,
                    )
                    plot_prediction_over_time(
                        prediction_df = predictions_time,
                        color_dict = color_dict,
                        output_file = output_pdf,
                    )

                else:
                    logging.info('Skipping plotting results')


        logging.info(
            '''
            No new bam files found, sleeping for {} seconds
            '''.format(cooldown)
        )
        time.sleep(cooldown)