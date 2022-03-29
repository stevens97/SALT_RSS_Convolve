'''
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SALT RSS Image Convolution

This program convolves 2 SALT RSS spectra images to the lowest resolution
of the two.

The files required are:
(1) The 2D spectrum FITS file of the first source.
(2) The 2D ARC FITS file of the first source.
(3) The 2D spectrum FITS file of the second source.
(4) The 2D ARC FITS file of the second source.

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
'''

'''
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
Import Libraries
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
'''

import os # For bash commands
from scipy.ndimage import gaussian_filter # For gaussian filtering
from astropy.io import fits as fits # For FITS file handling
import numpy as np # For array handling
from scipy.optimize import curve_fit # For curve fitting
from pyraf import iraf # For IRAF commands
from pathlib import Path # To extract filenames

'''
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
Load IRAF Libraries
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
'''

iraf.images()
iraf.images.imutil()
iraf.images.imgeom()
iraf.images.immatch()


'''
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
Read FITS file
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
'''
def read_FITS(file, centre):
    '''

    Opens the given FITS file and extracts the central spectrum of the astronomical source.

    :param file [String]: File name of FITS file to be opened.
    :param centre [Float]: Number of the central pixel of the astronomical source on the FITS image.

    :return: wave [Float array]: Wavelength.
    :return: flux [Float array]: Flux at current wavelength.
    '''

    # Set centre pixel
    centre = int(float(centre))

    # Get basename of file
    basename = Path(file).stem

    # Create temporary file
    tmp = '{}_tmp.fits'.format(basename)
    if os.path.isfile('tmp_ctr.fits'):
        iraf.images.imutil.imdelete(images=tmp, verify='No', mode='ql')

    # Extract central aperture of the FITS file
    iraf.images.imgeom.blkavg(input='{}[*,{}]'.format(file, centre), output=tmp, option='average',
                              b1=1, b2=1, b3=1, b4=1, b5=1, b6=1, b7=1, mode='ql')

    # Open FITS file and extract wavelength & flux
    hdu = fits.open(tmp)
    hdr = hdu[0].header
    flux = hdu[0].data
    flux = np.array(flux, dtype=np.float64)
    start_wave = hdr['CRVAL1']  # Initial wavelenghth
    step = hdr['CDELT1']  # Increment per pixel
    w0, dw, n = start_wave, step, len(flux)
    w = start_wave + step * n
    wave = np.linspace(w0, w, n, endpoint=False)
    if os.path.isfile(tmp):
        iraf.images.imutil.imdelete(images=tmp, verify='No', mode='ql')

    return wave, flux

'''
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
FWHM Calculation with 4333 Angstrom Line
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
'''
def FWHM_calc(ARC, centre):
    '''

    Calculates the instrumental FWHM resolution of the telescope (in this case the Southern African Large Telescope),
    given the ARC FITS file. The FWHM is calculated based off the 4333 Angstrom line. This code will only work if the given
    ARC file includes the 4333 Angstrom emission line.

    :param ARC [String]: File name of ARC FITS file to be opened.
    :param centre [Float]: Number of the central pixel of the astronomical source on the FITS image.

    :return: FWHM [float]: Instrumental resolution of the telescope.
    :return: FWHM [float]: Error of the instrumental resolution of the telescope.

    '''

    # Read FITS file
    wave, flux = read_FITS(ARC, centre)

    # Set wavelength range for fitting the Gaussian.
    # In this case we fit the 4333 Angstrom emission line.
    line_4333 = [(wave >= 4320) & (wave <= 4340)]
    flux = flux[tuple(line_4333)]
    wave = wave[tuple(line_4333)]

    # Determine which data points have a signal above 5 sigma.
    n = len(flux)
    noise = 0.6052697 * np.median(np.abs(2.0 * flux[2:n - 2] - flux[0:n - 4] - flux[4:n]))
    flux_base = flux[flux < 5 * noise]
    wave_base = wave[flux < 5 * noise]

    # Fit the baseline (all data points with a signal less than 5 sigma are considered to be the baseline).
    def straight_line(x, m, c):
        return m * x + c
    coeff, var_matrix = curve_fit(straight_line, wave_base, flux_base)
    m = coeff[0]
    c = coeff[1]
    baseline = m * wave + c

    # Subtract the baseline
    smooth_flux = flux - baseline

    # Fit the Guassian
    def Gauss(x, a, mu, sigma):
        return a * np.exp(-(x - mu) ** 2 / (2 * sigma ** 2))
    amp = np.max(flux)
    mean = np.mean(wave)
    sigma = np.std(wave)
    coeff, var_matrix = curve_fit(Gauss, wave, smooth_flux, p0=[amp, mean, sigma])
    sigma = coeff[2]
    sigma_err = var_matrix[2][2]
    FWHM = 2.355 * sigma
    FWHM_err = 2.355 * sigma_err

    return FWHM, FWHM_err


'''
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
Convolve Images
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
'''

def convolve(work_dir, src1, ARC1, ctr1, src2, ARC2, ctr2):
    '''

    Convolve two image sources together, such that the two images match in resolution.
    The higher resolution image is set to match the lower resolution image.

    :param src1 [String]: Complete path and filename of the first image source (FITS file).
    :param ARC1 [String]: Complete path and filename of the first ARC image (FITS file).
    :param ctr1 [Float]: Number of the central pixel of the first astronomical source.
    :param src2 [String]: Complete path and filename of the second image source (FITS file).
    :param ARC2 [String]: Complete path and filename of the first ARC image (FITS file).
    :param ctr1 [Float]: Number of the central pixel of the second astronomical source.

    :return: None
    '''

    # Open first FITS image
    hdu_1 = fits.open('{}'.format(src1))

    # Open second FITS image
    hdu_2 = fits.open('{}'.format(src2))

    # Calculate instrumental FWHM resolutions of both images
    FWHM_1, FWHM_1_err = FWHM_calc(ARC1, ctr1)
    FWHM_2, FWHM_2_err = FWHM_calc(ARC2, ctr2)

    print('FWHM_1 = {} +\- {}'.format(FWHM_1, FWHM_1_err))
    print('FWHM_2 = {} +\- {}'.format(FWHM_2, FWHM_2_err))

    # Calculate the difference
    delta_sigma = np.abs(FWHM_1 - FWHM_2)
    print('FWHM Difference = {}'.format(delta_sigma))

    # Convolve the images to match the lowest resolution
    if FWHM_1 < FWHM_2:
        for i in range(len(hdu_1[0].data)):
            hdu_1[0].data[i] = gaussian_filter(hdu_1[0].data[i], sigma=delta_sigma)
    elif FWHM_2 < FWHM_1:
        for i in range(len(hdu_2[0].data)):
            hdu_2[0].data[i] = gaussian_filter(hdu_2[0].data[i], sigma=delta_sigma)

    # Get basenames of source files
    name1 = Path(src1).stem
    name2 = Path(src2).stem
    new1 = '{}{}_new.fits'.format(work_dir, name1)
    new2 = '{}{}_new.fits'.format(work_dir, name2)

    # Write new files
    hdu_1.writeto(new1)
    hdu_2.writeto(new2)

    return None


'''
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
Run Program
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
'''

# Set work_dir
work_dir = raw_input('\nPlease enter the filepath where you would like to save the images to:  ')

# Get 1st Galaxy File
src1 = raw_input('\nEnter file path & name of the first GALAXY FITS file (with full path and extension):  ')

# Get 1st ARC File
ARC1 = raw_input('\nEnter file path & name of the first ARC FITS file (with full path and extension):  ')

# Get 1st File centre pixel
ctr1 = raw_input('\nEnter the centre pixel of the first galaxy:  ')

# Get 2nd Galaxy File
src2 = raw_input('\nEnter file path & name of the second GALAXY FITS file (with full path and extension):  ')

# Get 2nd ARC File
ARC2 = raw_input('\nEnter file path & name of the second ARC FITS file (with full path and extension):  ')

# Get 2nd File centre pixel
ctr2 = raw_input('\nEnter the centre pixel of the second galaxy:  ')

# Convolve Files
convolve(work_dir, src1, ARC1, ctr1, src2, ARC2, ctr2)
