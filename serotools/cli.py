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

    description = """Tools for accessing the White-Kauffmann-Le Minor (WKL) Salmonella 
                     serotyping scheme and for comparing serovars for congruency."""
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-v", "--version", action="version", version="%(prog)s version " + __version__)
    subparsers = parser.add_subparsers(dest="subparser_name", help=None, metavar="subcommand")
    subparsers.required = True

    formatter_class = argparse.ArgumentDefaultsHelpFormatter

    help_str = """Query the WKL database with one or more serovar names or antigenic formulas."""
    description = help_str
    subparser = subparsers.add_parser("query", formatter_class=formatter_class, description=description, help=help_str)
    subparser.add_argument("-i", "--input",   dest="in_file", type=str,            help="Specify an input file with one query (serovar or antigenic formula) per line.")
    subparser.add_argument("-s", "--serovar", dest="serovar", type=str,            help="Specify a query (serovar name or antigenic formula).")
    subparser.add_argument("-e", "--exact",   dest="exact",   action="store_true", help="Find exact matches only.")
    subparser.set_defaults(func=query_command)

    help_str = "Compare one of more pairs of serovars for congruency."
    description = help_str
    subparser = subparsers.add_parser("compare", formatter_class=formatter_class, description=description, help=help_str)
    subparser.add_argument(      "--header",dest="header",  action="store_true",  help="The input file includes a header line.")
    subparser.add_argument("-i", "--input", dest="in_file", type=str,             help="Specify a tab-delimited input file with two columns of serovars for comparison.")
    subparser.add_argument("-1", "--subj",  dest="subj",    type=str,             help="Specify the first serovar for comparison.")
    subparser.add_argument("-2", "--query", dest="query",   type=str,             help="Specify the second serovar for comparison.")
    subparser.set_defaults(func=compare_command)

    help_str = """Determine the most abundant serovar(s) for one or more clusters of isolates."""
    description = help_str
    subparser = subparsers.add_parser("cluster", formatter_class=formatter_class, description=description, help=help_str)
    subparser.add_argument("-i", "--input",     dest="in_file", type=str, help="Specify a tab-delimited input file in which each line contains two fields: a cluster id and a serovar designation, respectively.")
    subparser.add_argument("-s", "--sortby",    dest="sort_by", type=str, help="One or more comma-delim options for ordered sort results. Options = m (min_con), c (congruent), e (exact), i (init). Default = c,e,i.")                               
    subparser.add_argument("-v", "--verbosity", dest="v",       type=int, help="Verbosity of output. 1 - serovar info. 2 - serovar and abundance. 3 - all metrics.")
    subparser.set_defaults(func=cluster_command)

    args = parser.parse_args(system_args)
    return args


def cluster_command(args):
    """Determine the most abundant serovar(s) for one or more clusters of closely related 
       isolates.
    Parameters
    ----------
    args : Namespace
        Command line arguments stored as attributes of a Namespace, usually
        parsed from sys.argv, but can be set programmatically for unit testing
        or other purposes.
    """
    sero.cluster(args.in_file, args.sort_by, args.v)


def compare_command(args):
    """Compare one of more pairs of serovars for congruency.
    Parameters
    ----------
    args : Namespace
        Command line arguments stored as attributes of a Namespace, usually
        parsed from sys.argv, but can be set programmatically for unit testing
        or other purposes.
    """
    sero.compare(args.in_file, args.subj, args.query, args.header)


def query_command(args):
    """Query the WKL database with one or more serovar names or antigenic formulas.
    ----------
    args : Namespace
        Command line arguments stored as attributes of a Namespace, usually
        parsed from sys.argv, but can be set programmatically for unit testing
        or other purposes.
    """
    sero.query(args.in_file, args.serovar, args.exact)


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
