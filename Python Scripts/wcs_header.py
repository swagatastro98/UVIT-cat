import numpy as np
import pandas as pd
from glob import glob
from astropy.io import fits
from astropy.wcs import WCS

def wcs_header():
    
    path_to_fuv = ''

    path_to_nuv = '/Users/swagat98/Documents/Combine_cat/NUV-cat/' # path for the FUV files
    path_to_fuv = '/Users/swagat98/Documents/Combine_cat/FUV-cat_original/' # path for the NUV files

    fuv = glob(f'{path_to_fuv}*.fits')
    nuv = glob(f'{path_to_nuv}*.fits')

    count = 0
    for x in fuv:
        
        imgfuv = fits.open(x)

        header = imgfuv[0].header
        bin_hdr = imgfuv[1].header

        obj = bin_hdr['OBJECT']
        filter = bin_hdr['FILTER']
        ra_point = bin_hdr['RA_PNT']
        dec_point = bin_hdr['DEC_PNT']
        # exp_time = bin_hdr['EXP_TIME']

        binary_extension = {'OBJECT':obj,'filter':filter,'RA_PNT':ra_point,'DEC_PNT':dec_point}

        for key,value in binary_extension.items():
            header[key] = value
        # create new fits

        primary = fits.PrimaryHDU(header=header)

        final = fits.HDUList(primary)

        count += 1
    # print(bin_hdr.keys)
        final.writeto(f'/Users/swagat98/Documents/UVIT-cat/Corrected_FITS_header/{obj}_{count}_fuv_wcs_hdr.fits',overwrite=True)
    print('FUV WCS Files created')

    count1 = 0
        
    for y in nuv:
        imgfuv = fits.open(y)

        header = imgfuv[0].header
        bin_hdr = imgfuv[1].header

        obj = bin_hdr['OBJECT']
        filter = bin_hdr['FILTER']
        ra_point = bin_hdr['RA_PNT']
        dec_point = bin_hdr['DEC_PNT']
        # exp_time = bin_hdr['EXP_TIME']

        binary_extension = {'OBJECT':obj,'filter':filter,'RA_PNT':ra_point,'DEC_PNT':dec_point}

        for key,value in binary_extension.items():
            header[key] = value
        # create new fits

        primary = fits.PrimaryHDU(header=header)
    #
        # final = fits.HDUList(primary)
        count += 1 
    # print(bin_hdr.keys)
        final.writeto(f'/Users/swagat98/Documents/UVIT-cat/Corrected_FITS_header/{obj}_{count1}_nuv_wcs_hdr.fits',overwrite=True)
    print('NUV WCS Files created')







wcs_header()
