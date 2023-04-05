"""
Main program to access all subprograms
"""
import sys
import os
import time
import argparse

from sturgeon import __version__
from sturgeon.logger import setup_logging

def run():

    parser = argparse.ArgumentParser(
        prog="Sturgeon",
        description='''
        Sturgeon - CNS classifer - 
        \n
        \n
        This software is provided as indicated in the LICENSE file, by using
        this software you agree with the terms in the LICENSE file. - 
        \n
        \n
        © 2022 de Ridder lab
        ''',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version="Sturgeon version: {}".format(__version__),
        help="Show Sturgeon version and exit.",
    )
    parser.add_argument(
        "--no-logfile",
        action='store_true',
        help="Will not write the log to a file",
    )
    parser.set_defaults(func=lambda _: parser.print_help())

    from sturgeon.parsers import (
        register_predict,
        register_live,
        register_inputtobed,
        register_models,
    )
    
    subparsers = parser.add_subparsers(title="sub-commands")
    register_predict(subparsers)
    register_live(subparsers)
    register_inputtobed(subparsers)
    register_models(subparsers)

    args = parser.parse_args()
    
    cmd_func = args.func

    # skip making a logfile
    if args.no_logfile:
        log_file = None
    else:
        log_dir = os.path.join(os.getcwd(),'logs')
        log_file = os.path.join(log_dir, time.strftime("%Y%m%d-%H%M%S") + '.log')
        if not os.path.exists(log_dir):
            os.mkdir(log_dir)

    log_setup_success = setup_logging(
	    logfile_file = log_file,
    )
    if not log_setup_success:
        print('Failed to setup the log')
        sys.exit(1)
    
    cmd_func(args)
    

if __name__ == "__main__":
    run()