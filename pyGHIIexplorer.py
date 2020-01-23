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

def get_flux_min(Ha_map, plot=False):
    map_data = Ha_map
    mask_flux = map_data != 0
    mask_neg_flux = map_data < 0
    flux = map_data[mask_flux]
    # flux_neg = map_data[mask_neg_flux]
    # flux_neg_std = np.std(flux_neg)

    values, bins = np.histogram(flux, 2000, density=True)
    center_bins = [bins[i] + (bins[i+1]-bins[i])/2 \
                   for i in np.arange(0, len(bins)-1)]
    center_bins = np.array(center_bins)
    fit_params, fit_errors = gaussian_fit(center_bins, values)
    (mu, sigma) = (fit_params[1], fit_params[2])
    if plot:
        fig = Figure()
        canvas = FigureCanvas(fig)
        ax = fig.add_subplot()
        n, bins, patches = ax.hist(flux, 2000, density=True, facecolor='green',
                                   alpha=0.75)
        # x_values = np.linspace(bins[0], bins[-1], 2000)
        l = ax.plot(center_bins, n, 'r--', linewidth=2,
                    label=r'$\mathrm{Gaussian\ fit:}\ \mu=%.5f,\ \sigma=%.5f$'\
                    %(mu, sigma))
        ax.scatter(center_bins, values)
        ax.legend(loc='best')
        ax.set_xlabel('Flux')
        ax.set_ylabel('Number of pixels')
        ax.set_xlim([-0.015, 0.03])
        ax.set_ylim([0, 300])
        canvas.print_figure("histogram_sigma_flux.png", format="png")
    return np.abs(sigma)
    
def gaussian(x, amp, xc, sigma):
    return amp*np.exp( -(x-xc)**2 / (2*sigma**2)) / np.sqrt(2*np.pi*sigma**2)
        
def gaussian_fit(x, y, sigma = 0.001):
    # define some initial guess values for the fit routine
    # sigma = 0.001
    amp = x.max()
    xc = 0
    guess_vals = [amp, xc, sigma]

    # perform the fit and calculate fit parameter errors from covariance matrix
    fit_params, cov_mat = curve_fit(gaussian, x, y, p0=guess_vals)
    fit_errors = np.sqrt(np.diag(cov_mat))

    # manually calculate R-squared goodness of fit
    fit_residual = y - gaussian(x, *fit_params)
    fit_Rsquared = 1 - np.var(fit_residual)/np.var(y)
    return fit_params, fit_errors
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
