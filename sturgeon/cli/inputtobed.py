import os
from typing import Optional, List
import logging

from sturgeon.callmapping import bam_path_to_bed, mega_path_to_bed
from sturgeon.utils import validate_megalodon_file
import pysam

def filetobed(
    input_path: List[str],
    output_path: str,
    source: str,
    probes_file: str,
    margin: Optional[int] = 25,
    neg_threshold: Optional[float] = 0.3,
    pos_threshold: Optional[float] = 0.7,
):

    logging.info("Sturgeon start up")
    logging.info("File to bed program")

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

    if not os.path.exists(output_path):
        os.makedirs(output_path)

    if source == 'guppy':

        bamtobed(
            input_path = input_path,
            output_path = output_path,
            probes_file = probes_file,
            margin = margin,
            neg_threshold = neg_threshold,
            pos_threshold = pos_threshold,
        )
    
    elif source == 'megalodon':

        megatobed(
            input_path = input_path,
            output_path = output_path,
            probes_file = probes_file,
            margin = margin,
            neg_threshold = neg_threshold,
            pos_threshold = pos_threshold,
        )


def bamtobed(
    input_path: List[str],
    output_path: str,
    probes_file: str,
    margin: Optional[int] = 25,
    neg_threshold: Optional[float] = 0.3,
    pos_threshold: Optional[float] = 0.7,
):
    
    
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

    bam_path_to_bed(
        input_path = bam_files,
        output_path = output_path,
        probes_file = probes_file,
        margin = margin,
        neg_threshold = neg_threshold,
        pos_threshold = pos_threshold,
    )




def megatobed(
    input_path: List[str],
    output_path: str,
    probes_file: str,
    margin: Optional[int] = 25,
    neg_threshold: Optional[float] = 0.3,
    pos_threshold: Optional[float] = 0.7,
):

    logging.info("Megalodon to bed program")

    txt_files = list()
    if os.path.isfile(input_path):
        txt_files.append(input_path)
    elif os.path.isdir(input_path):
        for f in os.listdir(input_path):
            if not f.endswith('.txt'):
                continue
            txt_files.append(os.path.join(input_path, f)) 
    else:
        err_msg = '''
        --input-path must be a directory or file, given: {}
        '''.format(input_path)
        logging.error(err_msg)
        raise ValueError(err_msg)

    mega_files = list()
    for m in txt_files:
        success, msg = validate_megalodon_file(m)
        if success:
            mega_files.append(m)
        else:
            logging.error(
                '''
                File {}, did not pass validation either not a megalodon file or
                an invalid megalodon file. Reason: {}.
                '''.format(m, msg)
            )

    logging.info("Found a total of {} megalodon files".format(len(mega_files)))
    logging.info("Output will be saved in: {}".format(output_path))

    mega_path_to_bed(
        input_path = mega_files,
        output_path = output_path,
        probes_file = probes_file,
        margin = margin,
        neg_threshold = neg_threshold,
        pos_threshold = pos_threshold,
    )


    





    


    