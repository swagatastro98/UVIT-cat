from astropy.io import fits
import numpy as np
import pandas as pd
from astropy.table import Table
from glob import glob


def merge_all(path,output_file):
    
    files = glob(f'{path}*.fits')
    
    dflist = []


    img = fits.open(f'/Users/swagat98/Documents/UVIT-cat/UVIT-cat.fits')
    # img.info()
    # 
    #
    dat_main = img[1].data
    binary_header_main = img[1].header
    primary_header_main = img[0].header



    for x in files:
    
        img = fits.open(x)
    
        img.info()

        data = img[1].data
        header = img[1].header

        # obj = header['OBJECT']
        
        newdf = Table(data)
        df = newdf.to_pandas()
        
        dflist.append(df)
        
        print('For file ')
        print('The columns are: ')
        print(df.columns)

    #checks for nan value in each file before combining: #remove comment to run
        # for m in df.columns:
        #     
        #     print('For four dataframes: ')
        #     print(f'{m}:',df[m].isna().sum())
        
    final_df = pd.concat(dflist,ignore_index=True)
    
    print(final_df.shape)
    print(dflist[1].columns)
   
    # for y in final_df.columns:
        # 
        # print(f'Columns: {y}')
        # print('For final dataframe')
        # print(final_df[y].isna().sum())

    final_df.fillna(-999,inplace=True)

    final_df_to_tables = Table.from_pandas(final_df)
    

    primary = fits.PrimaryHDU(header=primary_header_main)
    binary = fits.BinTableHDU(data=final_df_to_tables,header=binary_header_main)

    final_cat = fits.HDUList([primary,binary])
    final_cat.writeto(output_file)

    
    
merge_all('/Users/swagat98/MAIN_UVIT_cat/','/Users/swagat98/Documents/UVIT-cat/UVIT-cat_merged.fits')
