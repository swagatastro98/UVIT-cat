
from astropy.io import fits
import numpy as np
import pandas as pd
from astropy import units as u
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy.coordinates import Angle

def test_cat(path):

    img = fits.open(path)

    pri_header = img[0].header
    pri_data = img[0].data
    header = img[1].header
    data = img[1].data

    newdf = Table(data)
    df = newdf.to_pandas()

    only_f1_df = df[df['F7_M_A']>0]

    print( df.head(500) )

    new500df = only_f1_df[:500]

    binary_df = Table.from_pandas(new500df)

    # all_f1_df = Table.from_pandas(only_f1_df) 

    primary = fits.PrimaryHDU(data=pri_data,header=pri_header)
    binary = fits.BinTableHDU(data=binary_df,header=header)

    final = fits.HDUList([primary,binary])
    final.writeto('/Users/swagat98/Documents/testdf_f7.fits',overwrite=True)




# test_cat('/Users/swagat98/Documents/UVIT-cat/UVIT-cat.fits')

def test(path):

    img = fits.open(path)

    pri_hdr = img[0].header
    header = img[1].header
    data = img[1].data

    # print(header.keys)

    tabledf = Table(data)
    df = tabledf.to_pandas()
    print(df.shape)

    # print(pri_hdr.keys)


# test('/Users/swagat98/Documents/UVIT-cat/UVIT-cat_merged.fits')


def add_columns(path):

    img = fits.open(path)

    data = img[1].data
    header = img[1].header
    
    newdf = Table(data)
    df = newdf.to_pandas()

    # print(df.columns)

    df['N5_M_A'] = np.nan
    df['N5_MER_A'] = np.nan
    df['N5_M_I'] = np.nan
    df['N5_MER_I'] = np.nan
    df['N5_F_A'] = np.nan
    df['N5_FER_A'] = np.nan

    df['F5_M_A'] = np.nan
    df['F5_MER_A'] = np.nan
    df['F5_M_I'] = np.nan
    df['F5_MER_I'] = np.nan
    df['F5_F_A'] = np.nan
    df['F5_FER_A'] = np.nan

    # df['']
    print(df.shape)

    print(df)

    # print(df['N5_M_A'].shape)
# add_columns('/Users/swagat98/UVIT-cat.fits')



def find_filterwise_sources(path):

    
    img = fits.open(path)

    pri_data = img[0].data
    pri_hdr  = img[0].header
    bin_data = img[1].data
    bin_hdr  = img[1].header

    # print(bin_hdr.keys)

    newdf = Table(bin_data)
    df    = newdf.to_pandas()

    print(f'Total sources in the dataframe: { df.shape }')

    f1_df = df[df['F1_M_A']>0]
    f2_df = df[df['F2_M_A']>0]
    f3_df = df[df['F3_M_A']>0]
    f5_df = df[df['F5_M_A']>0]
    f7_df = df[df['F7_M_A']>0]

    print(f'F1 length: {f1_df.shape}')
    print(f'F2 length: {f2_df.shape}')
    print(f'F3 length: {f3_df.shape}')
    print(f'F5 length: {f5_df.shape}')
    print(f'F7 length: {f7_df.shape}')

    n1_df = df[df['N1_M_A']>0]
    n2_df = df[df['N2_M_A']>0]
    n3_df = df[df['N3_M_A']>0]
    n5_df = df[df['N5_M_A']>0]
    n6_df = df[df['N6_M_A']>0]

    print(f'N1 length: {n1_df.shape}')
    print(f'N2 length: {n2_df.shape}')
    print(f'N3 length: {n3_df.shape}')
    print(f'N5 length: {n5_df.shape}')
    print(f'N6 length: {n6_df.shape}')






find_filterwise_sources('/Users/swagat98/Documents/UVIT-cat/UVIT-cat_merged.fits')
