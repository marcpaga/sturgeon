import os
import logging
from typing import List
import time
from pathlib import Path
from copy import deepcopy

import pandas as pd
import pysam

from sturgeon.bam import (
    bam_to_calls, 
    merge_probes_methyl_calls,
    probes_methyl_calls_to_bed,
)
from sturgeon.utils import validate_model_file, get_model_path
from sturgeon.prediction import predict_samples


def livebam(
    input_path: str,
    output_path: str,
    model_files: List[str],
    probes_file: str,
    margin: int,
    neg_threshold: float,
    pos_threshold: float,
    save_methyl_read_calls: bool,
    plot_results: bool,
    cooldown: int,
):
    """
    """

    logging.info("Sturgeon start up")
    logging.info("Live bam prediction program")

    logging.info('Starting live prediction from bam files')
    logging.info('Watching the following folder: {}'.format(input_path))

    bam_files = dict()

    probes_df = pd.read_csv(
        probes_file, 
        header = 0, 
        index_col = None, 
        sep = ' ',
    )

    probes_df['methylation_calls'] = 0
    probes_df['unmethylation_calls'] = 0
    probes_df['total_calls'] = 0

    while True:

        logging.info('Looking for new bam files, so far found {}'.format(len(bam_files)))
        for file in os.listdir(input_path):
            if not file.endswith('.bam'):
                continue
            
            file_name = Path(file).stem
            file_path = os.path.join(input_path, file)
            try:
                bam_files[file_name]
                continue
            except KeyError:
                logging.info('New bam file found: {}'.format(file_path))
                
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
                try:
                    pysam.index(file_path)
                except pysam.utils.SamtoolsError:
                    logging.warning('Making index file failed, bam file might not be complete yet, will try again later')
                    continue

                logging.info(
                    '''
                    Generated index file: {}
                    '''.format(bai_file)
                )

            logging.info('Getting methylation calls from file')
            calls_per_probe, calls_per_read = bam_to_calls(
                bam_file = file_path,
                probes_df = deepcopy(probes_df),
                margin = margin,
                neg_threshold = neg_threshold,
                pos_threshold = pos_threshold,
            )

            calls_per_probe_file = os.path.join(output_path, file_name + '_probes_methyl_calls.txt')
            calls_per_probe.to_csv(
                calls_per_probe_file,
                header = True, index = False, sep = '\t'
            )

            if save_methyl_read_calls:
                calls_per_read.to_csv(
                    os.path.join(output_path, file_name + '_read_methyl_calls.txt'), 
                    header = True, index = False, sep = '\t'
                )

            merged_output_file = os.path.join(output_path, 'merged_probes_methyl_calls_{}.txt'.format(len(bam_files)))
            if len(bam_files) == 0:
                merge_probes_methyl_calls(
                    [calls_per_probe_file], 
                    merged_output_file
                )
            else:
                merge_probes_methyl_calls(
                    [
                        calls_per_probe_file,
                        os.path.join(output_path, 'merged_probes_methyl_calls_{}.txt'.format(len(bam_files)-1))
                    ], 
                    merged_output_file
                )

            bed_output_file = os.path.join(
                output_path,
                'merged_probes_methyl_calls_{}.bed'.format(len(bam_files))
            )
            probes_methyl_calls_to_bed(
                merged_output_file,
                bed_output_file
            )

            bam_files[file_name] = file_path

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
                    bed_files = [bed_output_file],
                    model_file = model,
                    output_dir = output_path,
                    plot_results = plot_results,
                )

        logging.info('No new bam files found, sleeping for {} seconds'.format(cooldown))
        time.sleep(cooldown)

