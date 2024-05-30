from astropy.io import fits
import numpy as np
import pandas as pd
from astropy.table import Table

def cat_check(path):
    
    img = fits.open(path)
    
    data = img[1].data
    header = img[1].header
    
    tabledf = Table(data)
    df = tabledf.to_pandas()
    
    print(df.keys())
    # return df

    length_f1 = '' 
    length_f2 = ''
    length_f3 = ''
    length_f5 = ''
    length_f7 = ''

    length_n1 = ''
    length_n2 = ''
    length_n3 = ''
    length_n5 = ''
    length_n6 = ''

    # print(df[df['F1_MAG_AUTO']>0].shape )
    length_f1 += str( df[df['F1_MAG_AUTO']>0].shape ) 
    length_f2 += str( df[df['F2_MAG_AUTO']>0].shape ) 
    length_f3 += str( df[df['F3_MAG_AUTO']>0].shape)
    length_f5 += str( df[df['F5_MAG_AUTO']>0].shape)
    length_f7 += str( df[df['F7_MAG_AUTO']>0].shape)


    length_n1 += str( df[df['N1_MAG_AUTO']>0].shape)
    length_n2 += str( df[df['N2_MAG_AUTO']>0].shape)
    length_n3 += str( df[df['N3_MAG_AUTO']>0].shape)
    length_n5 += str( df[df['N5_MAG_AUTO']>0].shape)
    length_n6 += str( df[df['N6_MAG_AUTO']>0].shape)


    print(length_f1)
    print(length_f2)
    print(length_f3)
    print(length_f5)
    print(length_f7)

    print(length_n1)
    print(length_n2)
    print(length_n3)
    print(length_n5)
    print(length_n6)

    # print()
    # print(df[ df['F1_MAG_AUTO']<0 ] )
    print(df['F1_MAG_AUTO'].isna().sum()) # before removing the nan values 
    newdf = df.fillna(-9999)
    # print(newdf['F1_MAG_AUTO'].isna().sum())
    print(newdf['F1_MAG_AUTO'].isna().sum()) # after removing the nan values


    # print(newdf)


    print(header)    



cat_check('/Users/swagat98/Downloads/UVIT-cat_main.fits')
