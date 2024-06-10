import numpy as np
import pandas as pd
from glob import glob

def wcs_header():
    
    path_to_fuv = ''

    path_to_nuv = '/Users/swagat98/Documents/Combine_cat/NUV-cat/' # path for the FUV files
    path_to_fuv = '/Users/swagat98/Documents/Combine_cat/FUV-cat_original/' # path for the NUV files

    fuv = glob(f'{path_to_fuv}*.fits')
    nuv = glob(f'{path_to_nuv}*.fits')

    for x in fuv:
        
        imgfuv = fits.open(x)

        header = imgfuv[1].header







wcs_header()
