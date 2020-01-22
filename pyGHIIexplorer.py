#! /usr/bin/env python

__author__ = "Carlos Espinosa-Ponce"
__version__ = "0.1.0"
__license__ = "BSD"

import os
import argparse
from astropy.io import fits

def read_flux_elines_cubes(file_path, header=True):
    if header:
        data, header = fits.getdata(file_path, header=header)
        return header, data
    else:
        data = fits.getdata(file_path, header=header)
        return data

def read_fits_slice(file_path, n_param, header=True):
    if header:
        head, data = read_flux_elines_cubes(file_path, header=header)        
    else:
        data = read_flux_elines_cubes(file_path, header=header)
    map_slice = data[n_param, :, :]
    if header:
        return head, map_slice
    else:
        return map_slice

def main(args):
    fe_file = args.fe_file
    nHa = args.nHa
    print(fe_file)
    print(nHa)
    try:
        Ha_map = read_fits_slice(fe_file, nHa, header=False)
    except FileNotFoundError as err:
        print("File Not Found Error: {0}".format(err))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='pyGHIIexplorer', description='TBW')

    parser.add_argument("fe_file", help="Flux elines file", type=str)
    parser.add_argument("nHa", help="Channel of Ha", type=int)
    parser.add_argument("--version", action="version",
                        version="%(prog)s"\
                        "(version {version})".format(version=__version__))
    args = parser.parse_args()
    main(args)
