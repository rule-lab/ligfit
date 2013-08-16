# -*- coding: utf-8 -*-

"""LigFit

Usage:
  ligfit.py <input> [<output>]
  ligfit.py --help
  ligfit.py --version

Arguments:
    INPUT   input file path
    OUTPUT  output file path, defaults to ligfit_INPUT

Options:
  -h --help            show this help message and exit
  --version            show version and exit
  -v --verbose         print status messages
  -q --quiet           report only file names

"""
from docopt import docopt


def main():
    pass

if __name__ == '__main__':
    arguments = docopt(__doc__, version='0.0.1')
    main()