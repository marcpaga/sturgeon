"""
Main program to access all subprograms
"""

import argparse

from sturgeon import __version__

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

    # the path to the configuration file

    from sturgeon.parsers import (
        register_predict,
        register_watch,
    )
    
    subparsers = parser.add_subparsers(title="sub-commands")
    register_predict(subparsers)
    register_watch(subparsers)

    args = parser.parse_args()
    cmd_func = args.func
    cmd_func(args)


if __name__ == "__main__":
    run()