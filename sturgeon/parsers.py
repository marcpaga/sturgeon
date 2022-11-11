import os
import argparse

from sturgeon.utils import get_available_models

def register_predict(parser):

    subparser = parser.add_parser(
        'predict',
        description = 'Predict methylation samples in bed format',
        help = 'Predict methylation samples in bed format',
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
    )

    subparser.add_argument(
        '-i', '--input-path',
        type = str,
        required = True,
        help='Path to file (bed) or directory where the bed files are'
    )
    subparser.add_argument(
        '-o', '--output-path',
        type = str,
        required = True,
        help='Path where to save the results'
    )
    subparser.add_argument(
        '-m', '--model-files',
        type=str,
        required = True,
        nargs='+',
        help= '''
        Model file (zip) to be used to predict. More than one can be specified.
        These can be a path to the zip file, or one of the following built in
        models: {} '''.format(get_available_models(print_str = True)),
    )
    subparser.add_argument(
        '-p', '--plot-results',
        action='store_true',
        help='Also plot the results of the predictions'
    )

    subparser.set_defaults(func=run_predict)

def run_predict(args):

    from sturgeon.cli import predict

    predict.predict(
        input_path = args.input_path,
        output_path = args.output_path,
        model_files = args.model_files,
        plot_results = args.plot_results
    )

def register_livebam(parser):

    subparser = parser.add_parser(
        'livebam',
        description = '''
        Run a process that watches over a directory and predicts as the samples 
        are written into it. It assumes that each file written in the folder 
        comes from the same sample and are of a later timepoint.
        This program never ends, as it keeps watching the input folder for new
        files, therefore this has to be terminated manually.''',
        help = 'Predict methylation samples from format',
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
    )

    subparser.add_argument(
        '-i', '--input-path',
        type = str,
        required = True,
        help='Path where the bam files are written by the basecaller'
    )
    subparser.add_argument(
        '-o', '--output-path',
        type = str,
        required = True,
        help='Path where to save the results'
    )
    subparser.add_argument(
        '-m', '--model-files',
        type=str,
        required = True,
        nargs='+',
        help= '''
        Model file (zip) to be used to predict. More than one can be specified.
        These can be a path to the zip file, or one of the following built in
        models: {} '''.format(get_available_models(print_str = True)),
    )
    subparser.add_argument(
        '--probes-file',
        type = str,
        default = os.path.join(os.path.dirname(__file__), 'include/static', 'probes.bed'),
        help = 'Bed file with probe names and genomic locations'
    )
    subparser.add_argument(
        '--margin',
        type = int,
        default = 25,
        help='Neighbor methylation calls to consider when evaluating a probe location'
    )
    subparser.add_argument(
        '--neg-threshold',
        type = float,
        default = 0.3,
        help='Positions with scores below this threshold will be considered non-methylated'
    )
    subparser.add_argument(
        '--pos-threshold',
        type = float,
        default = 0.7,
        help='Positions with scores above this threshold will be considered methylated'
    )
    subparser.add_argument(
        '--save-methylreadcalls',
        action = 'store_true',
        help= 'Save as a text file methylation calls per read per bam file'
    )
    subparser.add_argument(
        '-p', '--plot-results',
        action='store_true',
        help='Also plot the results of the predictions'
    )
    subparser.add_argument(
        '--cooldown',
        type = int,
        default = 10,
        help = 'Seconds in between checking for a new bam file'
    )

    subparser.set_defaults(func=run_livebam)

def run_livebam(args):

    from sturgeon.cli import livebam

    livebam.livebam(
        input_path = args.input_path,
        output_path = args.output_path,
        model_files = args.model_files,
        probes_file = args.probes_file,
        margin = args.margin,
        neg_threshold = args.neg_threshold,
        pos_threshold = args.pos_threshold,
        save_methyl_read_calls = args.save_methylreadcalls,
        plot_results = args.plot_results,
        cooldown = args.cooldown,
    )

def register_bamtobed(parser):

    subparser = parser.add_parser(
        'bamtobed',
        description = 'Convert a modbam file, output of guppy to bed format',
        help = 'Convert a modbam file, output of guppy to bed format',
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
    )

    subparser.add_argument(
        '-i', '--input-path',
        type = str,
        required = True,
        help='Path to file (bam) or directory where the bam files are'
    )
    subparser.add_argument(
        '-o', '--output-path',
        type = str,
        required = True,
        help='Path where to save the output'
    )

    subparser.add_argument(
        '--probes-file',
        type = str,
        default = os.path.join(os.path.dirname(__file__), 'include/static', 'probes.bed'),
        help = 'Bed file with probe names and genomic locations'
    )
    subparser.add_argument(
        '--margin',
        type = int,
        default = 25,
        help='Neighbor methylation calls to consider when evaluating a probe location'
    )
    subparser.add_argument(
        '--neg-threshold',
        type = float,
        default = 0.3,
        help='Positions with scores below this threshold will be considered non-methylated'
    )
    subparser.add_argument(
        '--pos-threshold',
        type = float,
        default = 0.7,
        help='Positions with scores above this threshold will be considered methylated'
    )
    subparser.add_argument(
        '--save-methylreadcalls',
        action = 'store_true',
        help= 'Save as a text file methylation calls per read per bam file'
    )

    subparser.set_defaults(func=run_bamtobed)

def run_bamtobed(args):

    from sturgeon.cli import bamtobed

    bamtobed.bamtobed(
        input_path = args.input_path,
        output_path = args.output_path,
        probes_file = args.probes_file,
        margin = args.margin,
        neg_threshold = args.neg_threshold,
        pos_threshold = args.pos_threshold,
        save_methyl_read_calls = args.save_methylreadcalls
    )

def register_models(parser):

    subparser = parser.add_parser(
        'models',
        description = 'Get info on available models',
        help = 'Check and add models',
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
    )
    subparser.add_argument(
        "-a", "--action",
        type = str,
        choices = ['list', 'add', 'delete'],
        default = 'list',
        help = 'What to do, list the models, add models or delete models'
    )

    subparser.add_argument(
        '--model-files',
        type=str,
        required = False,
        nargs='+',
        help='Model files to be added or deleted'
    )

    subparser.set_defaults(func=run_models)


def run_models(args):

    from sturgeon.cli import models

    models.actions_models(
        action = args.action,
        model_files = args.model_files
    )
