"""
Main program to access all subprograms
"""
import logging
import sys
import os
import time
import argparse

from sturgeon import __version__
from sturgeon.logger import setup_logging

def run():

    parser = argparse.ArgumentParser(
        prog="Sturgeon",
        description="Sturgeon\nCNS classifer",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version="Sturgeon version: {}".format(__version__),
        help="Show Sturgeon version and exit.",
    )
    parser.set_defaults(func=lambda _: parser.print_help())

    from sturgeon.parsers import (
        register_predict,
        register_livebam,
        register_bamtobed,
        register_models,
    )
    
    subparsers = parser.add_subparsers(title="sub-commands")
    register_predict(subparsers)
    register_livebam(subparsers)
    register_bamtobed(subparsers)
    register_models(subparsers)

    args = parser.parse_args()
    
    cmd_func = args.func

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