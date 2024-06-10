
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table
import pandas as pd
from glob import glob
import os
import subprocess
def add_keyword_comment():

    


    img = fits.open( '/Users/swagat98/Documents/UVIT-cat/UVIT-cat_merged.fits',mode='update')
    img.info()

    data = img[1].data
    header = img[1].header
    primary_header = img[0].header
    primary_data = img[0].data

    binary = fits.BinTableHDU(data,header=header)

    binary.header.comments['TTYPE1'] = 'Source ID of the detected sources'
    binary.header.comments['TTYPE2'] = 'Right Ascension of source in J2000'
    binary.header.comments['TTYPE3'] = 'Declination of source in J2000'
    binary.header.comments['TTYPE4'] = 'FUV F1 Filter Magnitude Auto'
    binary.header.comments['TTYPE5'] = 'FUV F1 Filter Magnitude Auto Error'
    binary.header.comments['TTYPE6'] = 'FUV F1 Filter Isophotal Magnitude' 
    binary.header.comments['TTYPE7'] = 'FUV F1 Filter Isophotal Magnitude Error'
    binary.header.comments['TTYPE8'] = 'Extraction Flags'
    binary.header.comments['TTYPE9'] = 'FUV F1 Filter Flux Auto'
    binary.header.comments['TTYPE10'] = 'FUV F1 Filter Flux Auto Error'
    binary.header.comments['TTYPE11'] = 'Isophotal Image Major Axis (pixels)'
    binary.header.comments['TTYPE12'] = 'Isophotal Image Minor Axis (pixels)'
    binary.header.comments['TTYPE13'] = 'Isophotal Image Major Axis (World units)'
    binary.header.comments['TTYPE14'] = 'Isophotal Image Minor Axis (World units)'
    binary.header.comments['TTYPE15'] = 'Object position along X axis (pixels)'
    binary.header.comments['TTYPE16'] = 'Object position along Y axis (pixels)'
    binary.header.comments['TTYPE17'] = 'Object position along X axis (world units)'
    binary.header.comments['TTYPE18'] = 'Object position along Y axis (world units)'
    binary.header.comments['TTYPE19'] = 'Isophotal Image Major Axis error (pixels)'
    binary.header.comments['TTYPE20'] = 'Isophotal Image Minor Axis error (pixels)'
    binary.header.comments['TTYPE21'] = 'Isophotal Image Major Axis error (world units)'
    binary.header.comments['TTYPE22'] = 'Isophotal Image Minor Axis error (world units)'
    binary.header.comments['TTYPE23'] = 'Covariance of position between x and y'
    binary.header.comments['TTYPE24'] = 'FWHM of object assuming Gaussian core (degrees)'
    binary.header.comments['TTYPE25'] = 'FWHM of object assuming Gaussian core (pixels)'
    binary.header.comments['TTYPE26'] = 'Star/Galaxy Classifier'
    binary.header.comments['TTYPE27'] = 'Ellipticity'
    binary.header.comments['TTYPE28'] = 'Position angle in J2000'
    binary.header.comments['TTYPE29'] = 'Kron radius in units of A and B'
    binary.header.comments['TTYPE30'] = 'Galactic Right Ascension'
    binary.header.comments['TTYPE31'] = 'Galactic Declination'
    binary.header.comments['TTYPE32'] = 'Distance of object from centre of FOV'
    binary.header.comments['TTYPE33'] = 'FUV F2 Filter Magnitude Auto'
    binary.header.comments['TTYPE34'] = 'FUV F2 Filter Magnitude Auto Error'
    binary.header.comments['TTYPE35'] = 'FUV F2 Filter Isophotal Magnitude'
    binary.header.comments['TTYPE36'] = 'FUV F2 Filter Isophotal Magnitude Error'
    binary.header.comments['TTYPE37'] = 'FUV F2 Filter Flux Auto'
    binary.header.comments['TTYPE38'] = 'FUV F2 Filter Flux Auto Error'
    binary.header.comments['TTYPE39'] = 'FUV F3 Filter Magnitude Auto'
    binary.header.comments['TTYPE40'] = 'FUV F3 Filter Magnitude Auto Error'
    binary.header.comments['TTYPE41'] = 'FUV F3 Filter Isophotal Magnitude'
    binary.header.comments['TTYPE42'] = 'FUV F3 Filter Isophotal Magnitude Error'
    binary.header.comments['TTYPE43'] = 'FUV F3 Filter Flux Auto'
    binary.header.comments['TTYPE44'] = 'FUV F3 FILTER Flux Auto Error'
    binary.header.comments['TTYPE45'] = 'FUV F5 Filter Magnitude Auto'
    binary.header.comments['TTYPE46'] = 'FUV F5 Filter Magnitude Auto Error'
    binary.header.comments['TTYPE47'] = 'FUV F5 Filter Isophotal Magnitude'
    binary.header.comments['TTYPE48'] = 'FUV F5 Filter Isophotal Magnitude Error'
    binary.header.comments['TTYPE49'] = 'FUV F5 Filter Flux Auto'
    binary.header.comments['TTYPE50'] = 'FUV F5 Filter Flux Auto Error'
    binary.header.comments['TTYPE51'] = 'FUV F7 Filter Magnitude Auto'
    binary.header.comments['TTYPE52'] = 'FUV F7 Filter Magnitude Auto Error'
    binary.header.comments['TTYPE53'] = 'FUV F7 Filter Isophotal Magnitude'
    binary.header.comments['TTYPE54'] = 'FUV F7 Filter Isophotal Magnitude Error'
    binary.header.comments['TTYPE55'] = 'FUV F7 Filter Flux Auto'
    binary.header.comments['TTYPE56'] = 'FUV F7 Filter Flux Auto Error'
    binary.header.comments['TTYPE57'] = 'NUV F1 Filter Magnitude Auto'
    binary.header.comments['TTYPE58'] = 'NUV F1 Filter Magnitude Auto Error'
    binary.header.comments['TTYPE59'] = 'NUV F1 Filter Isophotal Magnitude'
    binary.header.comments['TTYPE60'] = 'NUV F1 Filter Isophotal Magnitude Error'
    binary.header.comments['TTYPE61'] = 'NUV F1 Filter Flux Auto'
    binary.header.comments['TTYPE62'] = 'NUV F1 Filter Flux Auto Error'
    binary.header.comments['TTYPE63'] = 'NUV F2 Filter Magnitude Auto'
    binary.header.comments['TTYPE64'] = 'NUV F2 Filter Magnitude Auto Error'
    binary.header.comments['TTYPE65'] = 'NUV F2 Filter Isophotal Magnitude'
    binary.header.comments['TTYPE66'] = 'NUV F2 Filter Isophotal Magnitude Error'
    binary.header.comments['TTYPE67'] = 'NUV F2 Filter Flux Auto'
    binary.header.comments['TTYPE68'] = 'NUV F2 Filter Flux Auto Error'
    binary.header.comments['TTYPE69'] = 'NUV F3 Filter Magnitude Auto'
    binary.header.comments['TTYPE70'] = 'NUV F3 Filter Magnitude Auto Error'
    binary.header.comments['TTYPE71'] = 'NUV F3 Filter Isophotal Magnitude'
    binary.header.comments['TTYPE72'] = 'NUV F3 Filter Isophotal Magnitude Error'
    binary.header.comments['TTYPE73'] = 'NUV F3 Filter Flux Auto'
    binary.header.comments['TTYPE74'] = 'NUV F3 Filter Flux Auto Error'
    binary.header.comments['TTYPE75'] = 'NUV F5 Filter Magnitude Auto'
    binary.header.comments['TTYPE76'] = 'NUV F5 Filter Magnitude Auto Error'
    binary.header.comments['TTYPE77'] = 'NUV F5 Filter Isophotal Magnitude'
    binary.header.comments['TTYPE78'] = 'NUV F5 Filter Isophotal Magnitude Error'
    binary.header.comments['TTYPE79'] = 'NUV F5 Filter Flux Auto'
    binary.header.comments['TTYPE80'] = 'NUV F5 Filter Flux Auto Error'
    binary.header.comments['TTYPE81'] = 'NUV F6 Filter Magnitude Auto'
    binary.header.comments['TTYPE82'] = 'NUV F6 Filter Magnitude Auto Error'
    binary.header.comments['TTYPE83'] = 'NUV F6 Filter Isophotal Magnitude'
    binary.header.comments['TTYPE84'] = 'NUV F6 Filter Isophotal Magnitude Error'
    binary.header.comments['TTYPE85'] = 'NUV F6 Filter Flux Auto'
    binary.header.comments['TTYPE86'] = 'NUV F6 Filter Flux Auto Error'


# we have to make edit to the BinTableHDU, and then the variable.header.comments to make the changes

    img.flush()    # this line is essential to update the existing keywords of the header values but also we cannot use the write to as it replaces the comments
    img.close()
    
    # binary = fits.BinTableHDU(data=primary_data,header=primary_header)
    primary = fits.PrimaryHDU(data=primary_data,header=primary_header)

    final = fits.HDUList([primary,binary])
    final.writeto('/Users/swagat98/Documents/UVIT-cat/UVIT-cat_merged.fits',overwrite=True)

add_keyword_comment()
