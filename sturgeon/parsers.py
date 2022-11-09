import os
import argparse

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
        help='Model file (zip) to be used to predict. More than one can be specified'
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

def register_watch(parser):

    subparser = parser.add_parser(
        'watch',
        description = '''
        Run a process that watches over a directory and predicts as the samples 
        are written into it. It assumes that each file written in the folder 
        comes from the same sample and are of a later timepoint''',
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
        help='Model file (zip) to be used to predict. More than one can be specified'
    )
    subparser.add_argument(
        '-p', '--plot-results',
        action='store_true',
        help='Also plot the results of the predictions'
    )

    subparser.set_defaults(func=run_watch)

def run_watch(args):

    from sturgeon.cli import watch

    watch(
        input_path = args.input_path,
        output_path = args.output_path,
        model_files = args.model_files,
        plot_results = args.plot_results
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
        '--processes',
        type = int,
        default = 1,
        help='Number of parallel processes to run'
    )

    subparser.set_defaults(func=run_bamtobed)

def run_bamtobed(args):

    from sturgeon.cli import bamtobed

    bamtobed.bamtobed(
        input_path = args.input_path,
        output_path = args.output_path,
        probes_file = os.path.join(os.path.dirname(__file__), 'static', 'probes.bed'),
        margin = args.margin,
        neg_threshold = args.neg_threshold,
        pos_threshold = args.pos_threshold,
        processes = args.processes
    )