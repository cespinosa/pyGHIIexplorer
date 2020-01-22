#! /usr/bin/env python

__author__ = "Carlos Espinosa-Ponce"
__version__ = "0.1.0"
__license__ = "BSD"

import argparse

def main(args):
    fe_file = args.fe_file
    nHa = args.nHa
    print(fe_file)
    print(nHa)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='pyGHIIexplorer', description='TBW')

    parser.add_argument("fe_file", help="Flux elines file", type=str)
    parser.add_argument("nHa", help="Channel of Ha", type=int)
    parser.add_argument("--version", action="version",
                        version="%(prog)s"\
                        "(version {version})".format(version=__version__))
    args = parser.parse_args()
    main(args)
