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
        '-m', '--model-file',
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

    subparser.set_defaults(fun=run_predict)

def run_predict(args):

    from sturgeon.cli import predict

    predict(
        input_path = args.input_path,
        output_path = args.output_path,
        model_file = args.model_file,
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
        '-m', '--model-file',
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

    subparser.set_defaults(fun=run_watch)

def run_watch(args):

    from sturgeon.cli import watch

    watch(
        input_path = args.input_path,
        output_path = args.output_path,
        model_file = args.model_file,
        plot_results = args.plot_results
    )