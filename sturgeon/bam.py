import os
from pathlib import Path
from copy import deepcopy
from typing import Optional, List
import logging

import pandas as pd
import numpy as np
from modbampy import ModBam

def get_methyl_calls_per_read( 
    bam_file: str, 
    chromosome: str,
) -> pd.DataFrame:
    """Get the methylation status of probes from a guppy methylation bam file

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

def map_methyl_calls_to_probes(
    chr_probes_df: pd.DataFrame,
    bam_file: str, 
    chromosome: str,
    margin: int, 
    neg_threshold: float,
    pos_threshold: float,
) -> pd.DataFrame:
    """Get the methylation status of probes from a guppy methylation bam file

    Args:
        chr_probes_df (pd.DataFrame): data.frame with the probe locations for
            one chromosome. Expected columns:
                - chr
                - start
                - end
                - methylation_calls
                - unmethylation_calls
                - total_calls
        bam_file (str): path to the bam file
        chromosome (str): chromosome to be evaluated
        margin (int): location margin around the probe location, methylation
            score will be averaged.
        neg_threshold (float): positions with scores lower than this threshold
            will be considered non-methylated
        pos_threshold (float): positions with scores higher than this threshold
            will be considered methylated

    Returns:
        Same pd.DataFrame as the input with the additional calls to each probe
    """

    methyl_calls_per_read = get_methyl_calls_per_read(
        bam_file,
        chromosome,
    )

    chr_probes_df.reset_index(inplace = True, drop = True)

    starts = np.array(chr_probes_df['start']) - margin
    ends = starts + margin * 2 + 1

    queries = np.array(methyl_calls_per_read['reference_pos'])
    scores = np.array(methyl_calls_per_read['score'])

    s = np.searchsorted(queries, starts, 'left')
    n = np.searchsorted(queries, ends, 'right')

    f = np.where(s != n)[0]

    s = s[f]
    n = n[f]

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
            chr_probes_df.loc[ff, 'methylation_calls'] += 1
        elif final_score == -1:
            chr_probes_df.loc[ff, 'unmethylation_calls'] += 1
        else:
            continue

    chr_probes_df.loc[:, 'total_calls'] = chr_probes_df.loc[:, 'methylation_calls'] + chr_probes_df.loc[:, 'unmethylation_calls']
    return chr_probes_df, methyl_calls_per_read

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

    calls_per_probe = list()
    calls_per_read = list()
    for chrom in chromosomes:
        calls_per_probe_df, calls_per_read_df = map_methyl_calls_to_probes(
            probes_df[probes_df['chr'] == chrom.item()], 
            bam_file, 
            chrom.item(),
            margin, 
            neg_threshold,
            pos_threshold,
        )

        chrom_num = np.unique(calls_per_probe_df['chr']).item()
        calls = calls_per_probe_df['total_calls'].sum()
        logging.info(
            '''
            Found a total of {} methylation calls on chromosome {}
            '''.format(calls, chrom_num)
        )
        calls_per_probe.append(calls_per_probe_df)
        calls_per_read.append(calls_per_read_df)
            
    calls_per_probe = pd.concat(calls_per_probe)
    calls_per_read = pd.concat(calls_per_read)

    return calls_per_probe, calls_per_read

def bam_to_bed(
    input_path: List[str],
    output_path: str,
    probes_file: str,
    margin: Optional[int] = 25,
    neg_threshold: Optional[float] = 0.3,
    pos_threshold: Optional[float] = 0.7,
):

    probes_df = pd.read_csv(
        probes_file, 
        header = 0, 
        index_col = None, 
        sep = ' ',
    )

    probes_df['methylation_calls'] = 0
    probes_df['unmethylation_calls'] = 0
    probes_df['total_calls'] = 0

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

                
        
