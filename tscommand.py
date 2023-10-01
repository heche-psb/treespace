#!/usr/bin/python3
import click
import logging
import sys
import os
import warnings
from timeit import default_timer as timer
import pkg_resources
from rich.logging import RichHandler
__version__ = pkg_resources.require("treespace")[0].version


# cli entry point
@click.group(context_settings={'help_option_names': ['-h', '--help']})
@click.option('--verbosity', '-v', type=click.Choice(['info', 'debug']),
    default='info', help="Verbosity level, default = info.")
def cli(verbosity):
    """
    treespace - Copyright (C) 2023-2024 Hengchi Chen\n
    Contact: heche@psb.vib-ugent.be
    """
    logging.basicConfig(
        format='%(message)s',
        handlers=[RichHandler()],
        datefmt='%H:%M:%S',
        level=verbosity.upper())
    logging.info("This is treespace v{}".format(__version__))
    pass


@cli.command(context_settings={'help_option_names': ['-h', '--help']})
@click.argument('config', type=click.Path(exists=True))
@click.option('--outdir', '-o', default='treespace_out', show_default=True, help='output directory')
def treerun(**kwargs):
    """
    Building gene trees
    """
    _treerun(**kwargs)

def _treerun(config,outdir):
    from treespace.tree import Config_Hauler
    start = timer()
    Config_Hauler(config,outdir)
    end = timer()
    logging.info("Total run time: {} min".format(round((end-start)/60,2)))

if __name__ == '__main__':
    cli()
