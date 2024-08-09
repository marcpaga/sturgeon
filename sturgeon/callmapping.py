import os
from pathlib import Path
from copy import deepcopy
from typing import Optional, List
import logging
from contextlib import contextmanager
import warnings

import pandas as pd
import numpy as np
try:
    from modbampy import ModBam
except ImportError:
    warnings.warn('Error loading modbampy, bam functionalities will not work')

from sturgeon.utils import read_probes_file

@contextmanager
def SuppressPandasWarning():
    with pd.option_context("mode.chained_assignment", None):
        yield


def get_methyl_calls_per_read(
    bam_file: str,
    chromosome: str,
) -> pd.DataFrame:
    """Get the methylation status of probes from a guppy methylation bam file
    of a chromosome.

    Args:
        bam_file (str): path to the bam file
        chromosome (str): chromosome to be evaluated

    Returns:
        pd.DataFrame with methylation calls per position
    """

    results = {
        'read_id': list(),
        'chr': list(),
        'reference_pos': list(),
        'strand': list(),
        'score': list()
    }

    if not isinstance(chromosome, str):
        chromosome = str(chromosome)

    if not chromosome.startswith('chr'):
        chromosome = 'chr'+str(chromosome)

    with ModBam(bam_file) as bam:

        st = 0
        nd = 250000000

        for read in bam.reads(chromosome, st, nd):

            for pos_mod in read.mod_sites:
                read_id, reference_pos, _, strand, _, _, _, score = pos_mod

                if reference_pos == -1:
                    continue

                results['read_id'].append(read_id)
                results['chr'].append(chromosome[3:])
                results['reference_pos'].append(reference_pos)
                results['strand'].append(strand)
                results['score'].append((1+score)/256)

    results = pd.DataFrame(results)
    results = results.sort_values(['reference_pos'])

    return results

def map_methyl_calls_to_probes_chr(
    probes_df: pd.DataFrame,
    methyl_calls_per_read: pd.DataFrame,
    margin: int,
    neg_threshold: float,
    pos_threshold: float,
) -> pd.DataFrame:
    """Maps calls per read to probe locations in a chromosome
    """

    probes_df = probes_df[probes_df['start'] > -1]
    probes_df = probes_df.sort_values(['start'])
    probes_df.reset_index(inplace = True, drop = True)
    methyl_calls_per_read = methyl_calls_per_read.sort_values(['reference_pos'])
    methyl_calls_per_read.reset_index(inplace = True, drop = True)

    starts = np.array(probes_df['start']) - margin
    ends = starts + margin * 2 + 1

    queries = np.array(methyl_calls_per_read['reference_pos'])
    scores = np.array(methyl_calls_per_read['score'])

    s = np.searchsorted(queries, starts, 'left')
    n = np.searchsorted(queries, ends, 'right')

    f = np.where(s != n)[0]

    s = s[f]
    n = n[f]

    with SuppressPandasWarning():
        for ss, nn, ff in zip(s, n, f):

            current_scores = scores[ss:nn]
            bin_scores = np.zeros(current_scores.shape)

            bin_scores[current_scores > pos_threshold] = 1
            bin_scores[current_scores < neg_threshold] = -1

            bin_scores = bin_scores[bin_scores != 0]
            if len(bin_scores) == 0:
                continue

            final_score = int(np.median(bin_scores))

            if final_score == 1:
                probes_df.loc[ff, 'methylation_calls'] += 1
            elif final_score == -1:
                probes_df.loc[ff, 'unmethylation_calls'] += 1
            else:
                continue

        probes_df['total_calls'] = probes_df['methylation_calls'] + probes_df['unmethylation_calls']

    return probes_df

def merge_probes_methyl_calls(
    file_list: List[str],
    output_file: str,
) -> pd.DataFrame:
    """Merge a list of methylation calls into a single one by adding all
    calls accross files

    Args:
        file_list (list): with the files to be merged
        output_file (str): file path to save the merged result
    """

    for file_num, file in enumerate(file_list):

        logging.info('Merging file: {}'.format(file))
        df = pd.read_csv(
            file,
            header = 0,
            index_col = None,
            sep = '\t',
        )
        df = df.sort_values(['chr', 'start'])

        if file_num == 0:
            output_df = deepcopy(df)
            continue

        for column in ['methylation_calls', 'unmethylation_calls', 'total_calls']:
            output_df[column] += df[column]

    logging.info('Saving merged file: {}'.format(output_file))
    output_df.to_csv(
        output_file, header = True, index = False, sep = '\t'
    )

    return output_df

def probes_methyl_calls_to_bed(
    input_file: str,
    output_file: str,
) -> pd.DataFrame:

    bed_df = {
        "chrom": list(),
        "chromStart": list(),
        "chromEnd": list(),
        "methylation_call": list(),
        "probe_id": list(),
    }

    probes_df = pd.read_csv(
        input_file,
        header = 0,
        index_col = None,
        sep = '\t',
    )

    probes_df = probes_df[probes_df['total_calls'] > 0]

    for _, row in probes_df.iterrows():
        if row.methylation_calls == row.unmethylation_calls:
            continue

        if row.methylation_calls > row.unmethylation_calls:
            m = 1
        else:
            m = 0

        bed_df['chrom'].append(row.chr)
        bed_df['chromStart'].append(row.start)
        bed_df['chromEnd'].append(row.end)
        bed_df['methylation_call'].append(m)
        bed_df['probe_id'].append(row.ID_REF)

    bed_df = pd.DataFrame(bed_df)

    logging.info('Total measured array CpG sites: {}'.format(bed_df.shape[0]))
    logging.info('Saving bed file to: {}'.format(output_file))
    bed_df.to_csv(
        output_file,
        header = True,
        index = False,
        sep = '\t'
    )

    return bed_df

def bam_to_calls(
    bam_file: str,
    probes_df: pd.DataFrame,
    margin: int,
    neg_threshold: float,
    pos_threshold: float,
):

    chromosomes = np.unique(probes_df['chr'])

    probes_df['methylation_calls'] = 0
    probes_df['unmethylation_calls'] = 0
    probes_df['total_calls'] = 0

    calls_per_probe = list()
    calls_per_read = list()
    for chrom in chromosomes:

        calls_per_read_df = get_methyl_calls_per_read(
            bam_file,
            chrom.item(),
        )

        calls_per_probe_df = map_methyl_calls_to_probes_chr(
            probes_df = probes_df[probes_df['chr'] == chrom.item()],
            methyl_calls_per_read = calls_per_read_df,
            margin = margin,
            neg_threshold = neg_threshold,
            pos_threshold = pos_threshold,
        )

        chrom_num = np.unique(calls_per_probe_df['chr']).item()
        calls = calls_per_probe_df['total_calls'].sum()
        logging.debug(
            '''
            Found a total of {} methylation calls on chromosome {}
            '''.format(calls, chrom_num)
        )
        calls_per_probe.append(calls_per_probe_df)
        calls_per_read.append(calls_per_read_df)

    calls_per_probe = pd.concat(calls_per_probe)
    calls_per_read = pd.concat(calls_per_read)

    return calls_per_probe, calls_per_read

def bam_path_to_bed(
    input_path: List[str],
    output_path: str,
    probes_file: str,
    margin: Optional[int] = 25,
    neg_threshold: Optional[float] = 0.3,
    pos_threshold: Optional[float] = 0.7,
):

    probes_df = read_probes_file(probes_file)

    output_files = list()

    for bam_file in input_path:

        logging.info(
            '''
            Getting methylation calls from: {}
            '''.format(bam_file)
        )

        bam_name = Path(bam_file).stem

        output_file = os.path.join(
            output_path,
            bam_name + '_probes_methyl_calls.txt'
        )
        output_files.append(output_file)

        if os.path.exists(output_file):
            logging.info(
                '''
                Skipping, probe_methyl_calls exists: {}
                '''.format(output_file)
            )
            continue

        probes_methyl_df = deepcopy(probes_df)
        logging.info('Processing bam file: {}'.format(bam_file))

        calls_per_probe, calls_per_read = bam_to_calls(
            bam_file = bam_file,
            probes_df = probes_methyl_df,
            margin = margin,
            neg_threshold = neg_threshold,
            pos_threshold = pos_threshold,
        )

        calls_per_probe.to_csv(
            output_file, header = True, index = False, sep = '\t'
        )
        ofile = os.path.join(
            output_path,
            bam_name + '_read_methyl_calls.txt'
        )
        calls_per_read.to_csv(
            ofile, header = True, index = False, sep = '\t'
        )

    merged_output_file = os.path.join(
        output_path,
        'merged_probes_methyl_calls.txt'
    )
    merge_probes_methyl_calls(
        output_files,
        merged_output_file
    )

    bed_output_file = os.path.join(
        output_path,
        'merged_probes_methyl_calls.bed'
    )
    probes_methyl_calls_to_bed(
        merged_output_file,
        bed_output_file
    )


def mega_file_to_bed(
    input_file: str,
    probes_file: str,
    margin: Optional[int] = 25,
    neg_threshold: Optional[float] = 0.3,
    pos_threshold: Optional[float] = 0.7,
) -> pd.DataFrame:

    probes_df = read_probes_file(probes_file)

    logging.info(
        '''
        Getting methylation calls from: {}
        '''.format(input_file)
    )

    mega_calls = pd.read_csv(
        input_file,
        header = 0,
        index_col = None,
        sep = '\t',
    )
    mega_calls = mega_calls.rename(columns = {
        'chrm':'chr',
        'pos': 'reference_pos',
    })
    mega_calls['score'] = np.exp(mega_calls['mod_log_prob'])
    mega_calls = mega_calls.drop(columns = [
        'mod_log_prob',
        'can_log_prob',
        'mod_base',
    ])

    probes_methyl_df = deepcopy(probes_df)
    logging.info('Processing megalodon file: {}'.format(input_file))

    chromosomes = np.unique(probes_df['chr'])

    probes_methyl_df['methylation_calls'] = 0
    probes_methyl_df['unmethylation_calls'] = 0
    probes_methyl_df['total_calls'] = 0

    calls_per_probe = list()
    for chrom in chromosomes:

        calls_per_probe_chr =  map_methyl_calls_to_probes_chr(
            probes_df =  probes_methyl_df[probes_methyl_df['chr'] == chrom.item()],
            methyl_calls_per_read = mega_calls[mega_calls['chr'] == 'chr'+str(chrom.item())],
            margin = margin,
            neg_threshold = neg_threshold,
            pos_threshold = pos_threshold,
        )
        calls_per_probe.append(calls_per_probe_chr)

        calls = calls_per_probe_chr['total_calls'].sum()
        logging.debug(
            '''
            Found a total of {} methylation array sites on chromosome {}
            '''.format(calls, chrom)
        )

    calls_per_probe = pd.concat(calls_per_probe)
    return calls_per_probe


def mega_path_to_bed(
    input_path: List[str],
    output_path: str,
    probes_file: str,
    margin: Optional[int] = 25,
    neg_threshold: Optional[float] = 0.3,
    pos_threshold: Optional[float] = 0.7,
):


    output_files = list()
    for mega_file in input_path:

        logging.info(
            '''
            Found methylation megalodon file: {}
            '''.format(mega_file)
        )

        mega_name = Path(mega_file).stem

        output_file = os.path.join(
            output_path,
            mega_name + '_probes_methyl_calls.txt'
        )
        output_files.append(output_file)

        if os.path.exists(output_file):
            logging.info(
                '''
                Skipping, probe_methyl_calls exists: {}
                '''.format(output_file)
            )
            continue

        calls_per_probe = mega_file_to_bed(
            input_file = mega_file,
            probes_file = probes_file,
            margin = margin,
            neg_threshold = neg_threshold,
            pos_threshold = pos_threshold,
        )

        calls_per_probe.to_csv(
            output_file, header = True, index = False, sep = '\t'
        )

    merged_output_file = os.path.join(
        output_path,
        'merged_probes_methyl_calls.txt'
    )
    merge_probes_methyl_calls(
        output_files,
        merged_output_file
    )

    bed_output_file = os.path.join(
        output_path,
        'merged_probes_methyl_calls.bed'
    )
    probes_methyl_calls_to_bed(
        merged_output_file,
        bed_output_file
    )


def modkit_file_to_bed(
    input_file: str,
    probes_file: str,
    margin: Optional[int] = 25,
    neg_threshold: Optional[float] = 0.3,
    pos_threshold: Optional[float] = 0.7,
    fivemc_code:str = 'm',
) -> pd.DataFrame:

    mandatory_columns = [
        "read_id",
        "chrom",
        "ref_position",
        "mod_qual",
        "mod_code",
        "canonical_base",
        "modified_primary_base",
    ]
    modkit_df = pd.read_csv(
        input_file,
        sep = '\t',
        header = 0,
        index_col = None,
        usecols= mandatory_columns
    )

    modkit_df = modkit_df[modkit_df['mod_code'] == fivemc_code]
    modkit_df = modkit_df[modkit_df['canonical_base'] == 'C']
    modkit_df = modkit_df[modkit_df['modified_primary_base'] == 'C']

    modkit_df = modkit_df.rename(columns={
        'chrom': 'chr',
        'ref_position': 'reference_pos',
        'mod_qual': 'score'
    })
    modkit_df = modkit_df.drop(columns=[
        'mod_code',
        'canonical_base',
        'modified_primary_base'
    ])
    modkit_df = modkit_df[modkit_df['reference_pos'] != -1]
    modkit_df = modkit_df[modkit_df['chr'] != '.']

    probes_df = read_probes_file(probes_file)

    probes_methyl_df = deepcopy(probes_df)
    logging.info('Processing modkit file: {}'.format(input_file))

    chromosomes = np.unique(probes_df['chr'])

    probes_methyl_df['methylation_calls'] = 0
    probes_methyl_df['unmethylation_calls'] = 0
    probes_methyl_df['total_calls'] = 0

    calls_per_probe = list()
    for chrom in chromosomes:

        calls_per_probe_chr =  map_methyl_calls_to_probes_chr(
            probes_df =  probes_methyl_df[probes_methyl_df['chr'] == chrom.item()],
            methyl_calls_per_read = modkit_df[modkit_df['chr'] == 'chr'+str(chrom.item())],
            margin = margin,
            neg_threshold = neg_threshold,
            pos_threshold = pos_threshold,
        )
        calls_per_probe.append(calls_per_probe_chr)

        calls = calls_per_probe_chr['total_calls'].sum()
        logging.debug(
            '''
            Found a total of {} methylation array sites on chromosome {}
            '''.format(calls, chrom)
        )

    calls_per_probe = pd.concat(calls_per_probe)
    return calls_per_probe


def modkit_path_to_bed(
    input_path: List[str],
    output_path: str,
    probes_file: str,
    margin: Optional[int] = 25,
    neg_threshold: Optional[float] = 0.3,
    pos_threshold: Optional[float] = 0.7,
    fivemc_code:str = 'c',
):

    output_files = list()
    for modkit_file in input_path:

        logging.info(
            '''
            Found methylation modkit file: {}
            '''.format(modkit_file)
        )

        modkit_name = Path(modkit_file).stem

        output_file = os.path.join(
            output_path,
            modkit_name + '_probes_methyl_calls.txt'
        )
        output_files.append(output_file)

        if os.path.exists(output_file):
            logging.info(
                '''
                Skipping, probe_methyl_calls exists: {}
                '''.format(output_file)
            )
            continue

        calls_per_probe = modkit_file_to_bed(
            input_file = modkit_file,
            probes_file = probes_file,
            margin = margin,
            neg_threshold = neg_threshold,
            pos_threshold = pos_threshold,
            fivemc_code = fivemc_code,
        )

        if calls_per_probe is None:
            logging.info(
                '''
                Modkit file processing gave an error, skipping: {}
                '''.format(modkit_file)
            )
            continue

        calls_per_probe.to_csv(
            output_file, header = True, index = False, sep = '\t'
        )

    merged_output_file = os.path.join(
        output_path,
        'merged_probes_methyl_calls.txt'
    )
    merge_probes_methyl_calls(
        output_files,
        merged_output_file
    )

    bed_output_file = os.path.join(
        output_path,
        'merged_probes_methyl_calls.bed'
    )
    probes_methyl_calls_to_bed(
        merged_output_file,
        bed_output_file
    )


def modkit_pileup_file_to_bed(
    input_file: str,
    output_file: str,
    probes_file: str,
    margin: Optional[int] = 25,
    neg_threshold: Optional[float] = 0.3,
    pos_threshold: Optional[float] = 0.7,
    fivemc_code:str = 'm',
) -> pd.DataFrame:

    column_names = [
        "chrom",
        "chromStart",
        "chromEnd",
        "mod_code",
        "score_bed",
        "strand",
        "thickStart",
        "thickEnd",
        "color",
        "valid_cov",
        "percent_modified",
        "n_mod",
        "n_canonical",
        "n_othermod",
        "n_delete",
        "n_fail",
        "n_diff",
        'n_nocall'
    ]
    modkit_df = pd.read_csv(
        input_file,
        delim_whitespace = True,
        header = None,
        index_col = None,
    )

    try:
        assert modkit_df.shape[1] == len(column_names)
    except AssertionError:
        err_msg = """
        Invalid modkit pileup file, number of columns does not match, expected {}, got {}.
        """.format(len(column_names), modkit_df.shape[1])

        raise AssertionError(err_msg)

    modkit_df.columns = column_names
    modkit_df = modkit_df[modkit_df['mod_code'] == fivemc_code]

    modkit_df = modkit_df.rename(columns={
        'chrom': 'chr',
        'chromStart': 'reference_pos',
        'percent_modified': 'score'
    })
    modkit_df['score'] /= 100

    modkit_df = modkit_df.drop(columns=[
        'mod_code',
        'thickStart',
        'thickEnd',
        'color',
        'valid_cov',
        'n_mod',
        "n_canonical",
        "n_othermod",
        "n_delete",
        "n_fail",
        "n_diff",
        'n_nocall'
    ])
    modkit_df = modkit_df[modkit_df['reference_pos'] != -1]
    modkit_df = modkit_df[modkit_df['chr'] != '.']

    probes_df = read_probes_file(probes_file)

    probes_methyl_df = deepcopy(probes_df)
    logging.info('Processing modkit file: {}'.format(input_file))

    chromosomes = np.unique(probes_df['chr'])

    probes_methyl_df['methylation_calls'] = 0
    probes_methyl_df['unmethylation_calls'] = 0
    probes_methyl_df['total_calls'] = 0

    calls_per_probe = list()
    for chrom in chromosomes:

        calls_per_probe_chr =  map_methyl_calls_to_probes_chr(
            probes_df =  probes_methyl_df[probes_methyl_df['chr'] == chrom.item()],
            methyl_calls_per_read = modkit_df[modkit_df['chr'] == 'chr'+str(chrom.item())],
            margin = margin,
            neg_threshold = neg_threshold,
            pos_threshold = pos_threshold,
        )
        calls_per_probe.append(calls_per_probe_chr)

        calls = calls_per_probe_chr['total_calls'].sum()
        logging.debug(
            '''
            Found a total of {} methylation array sites on chromosome {}
            '''.format(calls, chrom)
        )

    calls_per_probe = pd.concat(calls_per_probe)

    calls_per_probe.to_csv(
        output_file+'.tmp', header = True, index = False, sep = '\t'
    )

    calls_per_probe = calls_per_probe.rename(columns={
        'chr': 'chrom',
        'start': 'chromStart',
        'end': 'chromEnd',
        'ID_REF': 'probe_id',
        'methylation_calls': 'methylation_call'
    })

    calls_per_probe = calls_per_probe[calls_per_probe['total_calls'] > 0]
    calls_per_probe = calls_per_probe[
        ['chrom', 'chromStart', 'chromEnd', 'methylation_call', 'probe_id']
    ]

    calls_per_probe.to_csv(
        output_file, header = True, index = False, sep = '\t'
    )
