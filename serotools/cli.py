#!/usr/bin/env python3

"""Console script for SeroTools."""

import argparse
import logging
import sys

from serotools import serotools as sero
from serotools.__init__ import __version__

# Ignore flake8 errors in this module
# flake 8: noqa

def parse_arguments(system_args):
    """Parse command line arguments.

    Parameters
    ----------
    system_args : list
        List of command line arguments, usually sys.argv[1:].

    Returns
    -------
    Namespace
        Command line arguments are stored as attributes of a Namespace.
    """

    description = """Tools for accessing the White-Kauffmann-Le Minor (WKLM) Salmonella 
                     serotyping scheme and for comparing serovars for congruency."""
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-v", "--version", action="version", version="%(prog)s version " + __version__)
    subparsers = parser.add_subparsers(dest="subparser_name", help=None, metavar="subcommand")
    subparsers.required = True

    formatter_class = argparse.ArgumentDefaultsHelpFormatter

    help_str = """Query the WKLM database with a serovar name/formula in order to
                     retrieve a corresponding name/formula."""
    description = help_str
    subparser = subparsers.add_parser("query", formatter_class=formatter_class, description=description, help=help_str)
    subparser.add_argument("-i", "--input",   dest="in_file", type=str, metavar="INPUT_FILE",  help="Specify an input file with one query (serovar or antigenic formula) per line.")
    subparser.add_argument("-s", "--serovar", dest="serovar", type=str, metavar="SEROVAR", help="Specify a query (serovar name or antigenic formula).")
    subparser.set_defaults(func=query_command)

    help_str = "Compare two serovars for congruency."
    description = help_str + """Possible results are exact match, 
                     congruent, minimally congruent, or incongruent."""
    subparser = subparsers.add_parser("compare", formatter_class=formatter_class, description=description, help=help_str)
    subparser.add_argument(      "--header",dest="header",     action="store_true",           help="The input file includes a header line.")
    subparser.add_argument("-i", "--input", dest="input_file", type=str, metavar="INPUT_FILE",  help="Specify a tab-delimited input file with two columns of serovars for comparison.")
    subparser.add_argument("-1", "--subj",  dest="subj",       type=str, metavar="SEROVAR1",  help="Specify the first serovar for comparison.")
    subparser.add_argument("-2", "--query", dest="query",      type=str, metavar="SEROVAR2", help="Specify the second serovar for comparison.")
    subparser.set_defaults(func=compare_command)

    help_str = """Predict which serovars are mostly likely to be represented by 
                  (are minimally congruent with) a serovar name or antigenic formula."""
    description = help_str
    subparser = subparsers.add_parser("predict", formatter_class=formatter_class, description=description, help=help_str)
    subparser.add_argument("-i", "--input",   dest="in_file", type=str, metavar="INPUT_FILE",  help="Specify an input file with one query (serovar or antigenic formula) per line.")
    subparser.add_argument("-s", "--serovar", dest="serovar", type=str, metavar="SEROVAR", help="Specify a query (serovar name or antigenic formula).")
    subparser.set_defaults(func=predict_command)

    args = parser.parse_args(system_args)
    return args


def compare_command(args):
    """Compare two serovars for congruency. Possible results are exact match, 
       congruent, minimally congruent, or incongruent.""
    Parameters
    ----------
    args : Namespace
        Command line arguments stored as attributes of a Namespace, usually
        parsed from sys.argv, but can be set programmatically for unit testing
        or other purposes.
    """
    sero.compare(args.in_file, args.subj, args.query, args.header)


def predict_command(args):
    """Predicts which serovars are mostly likely to be represented by 
    (are minimally congruent with) a serovar name or antigenic formula.
    ----------
    args : Namespace
        Command line arguments stored as attributes of a Namespace, usually
        parsed from sys.argv, but can be set programmatically for unit testing
        or other purposes.
    """
    sero.predict(args.in_file, args.serovar)


def query_command(args):
    """Query the WKLM database with a serovar name/formula in order to
       retrieve corresponding name/formula.
    ----------
    args : Namespace
        Command line arguments stored as attributes of a Namespace, usually
        parsed from sys.argv, but can be set programmatically for unit testing
        or other purposes.
    """
    sero.query(args.in_file, args.serovar)


def run_command_from_args(args):
    """Run a subcommand with previously parsed arguments in an argparse namespace.
    This function is intended to be used for unit testing.
    Parameters
    ----------
    args : Namespace
        Command line arguments are stored as attributes of a Namespace.
        The args are obtained by calling parse_argument_list().
    Returns
    -------
    Returns 0 on success if it completes with no exceptions.
    """
    return args.func(args)  # this executes the function previously associated with the subparser with set_defaults


def run_from_line(line):
    """Run a command with a command line.
    This function is intended to be used for unit testing.
    Parameters
    ----------
    line : str
        Command line.
    Returns
    -------
    Returns 0 on success if it completes with no exceptions.
    """
    argv = line.split()
    args = parse_arguments(argv)
    return args.func(args)  # this executes the function previously associated with the subparser with set_defaults


def main():
    """This is the main function which is turned into an executable
    console script by the setuptools entry_points.  See setup.py.
    To run this function as a script, first install the package:
        $ python setup.py develop
        or
        $ pip install --user vcftoolz
    Parameters
    ----------
    This function must not take any parameters
    Returns
    -------
    The return value is passed to sys.exit().
    """
    enable_log_timestamps = False
    if enable_log_timestamps:
        logging.basicConfig(format="%(asctime)s %(message)s", datefmt="%Y-%m-%d %H:%M:%S", level=logging.INFO)
    else:
        logging.basicConfig(format="%(message)s", level=logging.INFO)
    args = parse_arguments(sys.argv[1:])
    return args.func(args)  # this executes the function previously associated with the subparser with set_defaults


# This snippet lets you run the cli without installing the entrypoint.
if __name__ == "__main__":
    sys.exit(main())
