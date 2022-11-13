import os
from typing import Optional, List
import logging

from sturgeon.bam import bam_to_bed
import pysam

def bamtobed(
    input_path: List[str],
    output_path: str,
    probes_file: str,
    margin: Optional[int] = 25,
    neg_threshold: Optional[float] = 0.3,
    pos_threshold: Optional[float] = 0.7,
):
    
    logging.info("Sturgeon start up")
    logging.info("Bam to bed program")

    bam_files = list()
    if os.path.isfile(input_path):
        bam_files.append(input_path)
    elif os.path.isdir(input_path):
        for f in os.listdir(input_path):
            if not f.endswith('.bam'):
                continue
            bam_files.append(os.path.join(input_path, f)) 
    else:
        err_msg = '''
        --input-path must be a directory or file, given: {}
        '''.format(input_path)
        raise ValueError(err_msg)

    if not os.path.exists(probes_file):
        err_msg = '''
        --probes-file not found, given: {}
        '''.format(probes_file)
        raise ValueError(err_msg)

    if neg_threshold < 0 or neg_threshold > 1:
        err_msg = '''
        --neg-threshold must be between 0 and 1, given: {}
        '''.format(neg_threshold)
        raise ValueError(err_msg)

    if pos_threshold < 0 or pos_threshold > 1:
        err_msg = '''
        --pos-threshold must be between 0 and 1, given: {}
        '''.format(pos_threshold)
        raise ValueError(err_msg)

    if pos_threshold <= neg_threshold:
        err_msg = '''
        --pos-threshold cannot be smaller or equal to -neg-threshold, given: 
        {} and {}
        '''.format(pos_threshold, neg_threshold)
        raise ValueError(err_msg)

    if margin < 0:
        err_msg = '''
        --margin must be zero or a positive integer, given: {}
        '''.format(margin)
        raise ValueError(err_msg)


    logging.info("Found a total of {} bam files".format(len(bam_files)))
    logging.info("Output will be saved in: {}".format(output_path))

    for bam_file in bam_files:
        bai_file = bam_file + '.bai'
        if not os.path.exists(bai_file):
            logging.info(
                '''
                Index file not found for bam file: {}
                '''.format(bam_file)
            )
            logging.info(
                '''
                Generating index file: {}
                '''.format(bai_file)
            )
            pysam.index(bam_file)

            logging.info(
                '''
                Generated index file: {}
                '''.format(bai_file)
            )

    if not os.path.exists(output_path):
        os.makedirs(output_path)

    bam_to_bed(
        input_path = bam_files,
        output_path = output_path,
        probes_file = probes_file,
        margin = margin,
        neg_threshold = neg_threshold,
        pos_threshold = pos_threshold,
    )







    





    


    