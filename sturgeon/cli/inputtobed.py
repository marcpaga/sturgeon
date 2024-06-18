import os
from typing import Optional, List
import logging
import warnings

from sturgeon.callmapping import (
    bam_path_to_bed, 
    mega_path_to_bed, 
    modkit_path_to_bed,
    modkit_pileup_file_to_bed
)
from sturgeon.utils import validate_megalodon_file, validate_modkit_file

try:
    import pysam
except ImportError:
    warnings.warn('Error loading modbampy, bam functionalities will not work')

def filetobed(
    input_path: List[str],
    output_path: str,
    source: str,
    probes_file: str,
    reference_genome: Optional[str],
    margin: Optional[int] = 25,
    neg_threshold: Optional[float] = 0.3,
    pos_threshold: Optional[float] = 0.7,
    fivemc_code: str = 'm',
):

    logging.info("Sturgeon start up")
    logging.info("File to bed program")

    if reference_genome is not None:
        probes_file = os.path.join(
            os.path.dirname(__file__), 
            '../include/static', 
            'probes_{}.bed'.format(reference_genome)
        )

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

        warnings.warn(
        '''\nUsing this source is NOT recommended.\n
        Please strongly consider using modkit (https://github.com/nanoporetech/modkit/), to extract methylation calls from bam files.\n 
        Otherwise we rely on modbampy, which is currently deprecated. Updates to guppy might break modbampy compatibility and could give WRONG results.                   
        ''')

        bamtobed(
            input_path = input_path,
            output_path = output_path,
            probes_file = probes_file,
            margin = margin,
            neg_threshold = neg_threshold,
            pos_threshold = pos_threshold,
        )
    
    elif source == 'megalodon':

        warnings.warn(
        '''\nUsing this source is NOT recommended.\n
        Please strongly consider using modkit (https://github.com/nanoporetech/modkit/), to extract methylation calls.\n 
        Megalodon is deprecated on ONT's end.                   
        ''')

        megatobed(
            input_path = input_path,
            output_path = output_path,
            probes_file = probes_file,
            margin = margin,
            neg_threshold = neg_threshold,
            pos_threshold = pos_threshold,
        )

    elif source == 'modkit':
        
        modkittobed_extract(
            input_path = input_path,
            output_path = output_path,
            probes_file = probes_file,
            margin = margin,
            neg_threshold = neg_threshold,
            pos_threshold = pos_threshold,
            fivemc_code = fivemc_code,
        )

    elif source == 'modkit_pileup':

        if os.path.exists(output_path):
            os.removedirs(output_path)

        modkittobed_pileup(
            input_path = input_path,
            output_path = output_path,
            probes_file = probes_file,
            margin = margin,
            neg_threshold = neg_threshold,
            pos_threshold = pos_threshold,
            fivemc_code = fivemc_code,
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

def modkittobed_extract(
    input_path: List[str],
    output_path: str,
    probes_file: str,
    margin: Optional[int] = 25,
    neg_threshold: Optional[float] = 0.3,
    pos_threshold: Optional[float] = 0.7,
    fivemc_code: str = 'm',
):
    
    logging.info("Modkit extract to bed program")

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

    modkit_files = list()
    for m in txt_files:
        success, msg = validate_modkit_file(m)
        if success:
            modkit_files.append(m)
        else:
            logging.error(
                '''
                File {}, did not pass validation either not a modkit file or
                an invalid modkit file. Reason: {}.
                '''.format(m, msg)
            )

    logging.info("Found a total of {} modkit files".format(len(modkit_files)))
    logging.info("Output will be saved in: {}".format(output_path))

    modkit_path_to_bed(
        input_path = modkit_files,
        output_path = output_path,
        probes_file = probes_file,
        margin = margin,
        neg_threshold = neg_threshold,
        pos_threshold = pos_threshold,
        fivemc_code = fivemc_code,
    )


def modkittobed_pileup(
    input_path: str,
    output_path: str,
    probes_file: str,
    margin: Optional[int] = 25,
    neg_threshold: Optional[float] = 0.3,
    pos_threshold: Optional[float] = 0.7,
    fivemc_code: str = 'm',
):
    
    logging.info("Modkit pileup to bed program")
    if not os.path.isfile(input_path):
        err_msg = '''
        --input-path must be file, given: {}
        '''.format(input_path)
        logging.error(err_msg)
        raise ValueError(err_msg)

    logging.info("Input modkit file: {}".format(input_path))
    logging.info("Output will be saved in: {}".format(output_path))

    modkit_pileup_file_to_bed(
        input_file = input_path,
        output_file = output_path,
        probes_file = probes_file,
        margin = margin,
        neg_threshold = neg_threshold,
        pos_threshold = pos_threshold,
        fivemc_code = fivemc_code,
    )



    


    