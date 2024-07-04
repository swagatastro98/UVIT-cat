
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table
import pandas as pd
from glob import glob
import os
import subprocess
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.coordinates import Angle


global path_to_fuv
global path_to_nuv
global same_field_folder_path

path_to_nuv = '/Users/swagat98/Documents/Combine_cat_2/NUV-cat/' # path for the FUV files
path_to_fuv = '/Users/swagat98/Documents/Combine_cat_2/FUV-cat/' # path for the NUV files
same_field_folder_path = '/Users/swagat98/Documents/Combine_cat_2/Same_field/'  # path where we should create a new same field file
same_field_path_fuv = '/Users/swagat98/Documents/Combine_cat_2/Same_field/FUV_same/'
same_field_path_nuv = '/Users/swagat98/Documents/Combine_cat_2/Same_field/NUV_same/'

path_to_fuv_2 = '/Users/swagat98/Documents/Combine_cat_2/FUV-cat_2/'




# create same field folder
# check for similar pointing images
    # in fuv folder first, and thn nuv folder
# copy to same field folder
    # for remaining images, create a catalog



def merged():

    if os.path.exists(same_field_folder_path):
        print('FOlder exists. Skipping..')
    else:
        print("Folder doesn't exist. Making one")
        os.mkdir(same_field_folder_path)

    if os.path.exists(same_field_path_fuv):
        print('Folder exists. SKipping..')
    else:
        print("Folder doesn't exists. Making one")
        os.mkdir(same_field_path_fuv)

    if os.path.exists(same_field_path_nuv):
        print('Folder exists. SKipping..')
    else:
        print("Folder doesn't exists. Making one")
        os.mkdir(same_field_path_nuv)


    fuv_1 = glob(f'{path_to_fuv}*.fits')


    dist = Angle(14,u.arcmin)
    print(dist)


    nuv_1 = glob(f'/Users/swagat98/Documents/Combine_cat_2/NUV-cat/*.fits')


    # print(type( fuv_2 ))

    # fuv_2.remove(f"{path_to_fuv_2}{'SMC2-F1-G08_006T02_9000001730-catalog.fits'}")
    # print(len(fuv_2))

    obj_names_to_deleted = []

    for x in fuv_1:



        img_fuv = fits.open(x)

        header_fuv = img_fuv[1].header

        ra_pnt = header_fuv['RA_PNT']
        dec_pnt = header_fuv['DEC_PNT']

        coord_fuv = SkyCoord(ra=ra_pnt,dec=dec_pnt,frame='icrs',unit='deg')
        # print(coord_fuv)
        obj = x.split('/')[-1]

        print('\n')
        print('Object File which is searching: ',obj)

        x = f'cp -r {path_to_fuv_2}{obj} /Users/swagat98/Documents/Combine_cat_2/test/' #copying the main files to test folder
        subprocess.run(x,shell=True,check=True)
        print('Copying main file')

# remove the file so that while searching the program won't find the file in the second folder
        y = f'rm -r {path_to_fuv_2}{obj}'
        subprocess.run(y,shell=True,check=True)
        print('Deleting main file')

        fuv_2 = glob(f'{path_to_fuv_2}*.fits') # keep it inside the loop so that the list returns back to its original state 

        # print(fuv_2)
        # The file which is searching will not be there in this list, and therefore, in the next loop should
        # not be considered


        for inside_fuv in range(len(fuv_2)):
            # print(len(fuv_2))
            # print('Objects inside the 2nd folder: ',fuv_2)
            # print(len(fuv_2))
            img_inside = fits.open(fuv_2[inside_fuv])
            # img_inside.info()

            header_inside_fuv = img_inside[1].header

            ra_in_fuv = header_inside_fuv['RA_PNT']
            dec_in_fuv = header_inside_fuv['DEC_PNT']
            obj_inside = fuv_2[inside_fuv].split('/')[-1].split('/')[0]

            coord_fuv_in = SkyCoord(ra=ra_in_fuv,dec=dec_in_fuv,frame='icrs',unit='deg')

            separation = coord_fuv_in.separation(coord_fuv)
            separation_in_deg = separation.to(u.arcmin)

            # # print(separation_in_deg)
            if separation_in_deg < dist:
                print('\t\tMatch found')
                print('Object: ',)
                print(separation_in_deg)

                print('Object searching inside',obj_inside)
                

                copy = f'cp -r {fuv_2[inside_fuv]} {same_field_path_fuv}' #similar images are separated out
                subprocess.run(copy,shell=True,check=True)
                print('Copying matched files')

                copy_2 = f'cp -r {path_to_fuv}{obj} {same_field_path_fuv}'
                subprocess.run(copy_2,shell=True,check=True)
                print('Copying original matching file')

                print('\n \n')

                obj_names_to_deleted.append(obj)

        
        for inside_nuv in range(len(nuv_1)):
            
            nuv_file = fits.open(nuv_1[inside_nuv])
            
            nuv_hdr = nuv_file[1].header

            ra_nuv = nuv_hdr['RA_PNT']
            dec_nuv = nuv_hdr['DEC_PNT']
            obj_nuv = nuv_1[inside_nuv].split('/')[-1]

            nuv_coord = SkyCoord(ra=ra_nuv,dec=dec_nuv,frame='icrs',unit='deg')
            nuv_sep = nuv_coord.separation(coord_fuv)
            nuv_sep_in_deg = nuv_sep.to(u.arcmin)

            if nuv_sep_in_deg < dist:
                print('\t\tNUV match found')
                print(f'Object: {obj_nuv}')
                print(nuv_sep_in_deg)

                copy_nuv = f'cp -r {nuv_1[inside_nuv]} {same_field_path_nuv}'
                subprocess.run(copy_nuv,shell=True,check=True)
                print('Copying NUV object')


 
    print(len(obj_names_to_deleted))

    # for m in obj_names_to_deleted:

        # remove from
        # print(m)
    
    path_to_test = '/Users/swagat98/Documents/Combine_cat_2/Same_field/FUV_same/'
    path_to_nuv = '/Users/swagat98/Documents/Combine_cat_2/Same_field/NUV_same/'

    test_files = glob(f'{path_to_test}*.fits')
    nuv_files = glob(f'{path_to_nuv}*.fits')

    # print(len(test_files))
    for m in test_files:

        test_obj = m.split('/')[-1]
        # print(test_obj)

        remove_cmd = f'rm -r /Users/swagat98/Documents/Combine_cat_2/test/{test_obj}'
        subprocess.run(remove_cmd,shell=True,check=True)

    for n in nuv_files:

        nuv_obj = n.split('/')[-1]

        remove_nuv = f'rm -r /Users/swagat98/Documents/Combine_cat_2/NUV-cat/{nuv_obj}'
        subprocess.run(remove_nuv,check=True,shell=True)
    

    
    




# merged()



## make filterwise catalog
## run merging script to find duplicate sources
##  one by one run the second layer of merging script between filters and add in columns



def remove_99(df):   #function to remove 99 values from the dataframe

    indexlist = []
    for i in range(len(df)):
        magnitude = df['MAG_AUTO'].iloc[i]
        if magnitude == 99:
            indexlist.append(i)
    
    df.drop(indexlist,axis='rows',inplace=True)
    return df



def dataframe_build_1():


    same_field_files = glob(f'{same_field_path_fuv}*.fits')
    # print(len( same_field_files ))


    f1_df = []
    f2_df = []
    f3_df = []
    f5_df = []
    f7_df = []






    for x in same_field_files:

        img = fits.open(x)

        pri_data = img[0].data
        pri_header = img[0].header

        data = img[1].data
        header = img[1].header
        obj = header['OBJECT']

        filter_name = header['FILTER']

        datatable = Table(data)
        df = datatable.to_pandas()

        if filter_name == 'F1':

            f1_remove_99 = remove_99(df)

            f1_remove_99.rename(columns={'MAG_AUTO':'F1_M_A',  # mag auto
                                               'MAGERR_AUTO':'F1_MER_A', #mag error auto 
                                               'MAG_ISO':'F1_M_I',  # mag iso
                                               'MAGERR_ISO':'F1_MER_I', #mag err iso
                                               'FLUX_AUTO':'F1_F_A', #
                                               'FLUXERR_AUTO':'F1_FER_A'},inplace=True)
            # print(f1_remove_99.columns)         
            f1_df.append(f1_remove_99)


        elif filter_name == 'F2':

            f2_remove_99 = remove_99(df)

            f2_remove_99.rename(columns={'MAG_AUTO':'F2_M_A',
                                               'MAGERR_AUTO':'F2_MER_A',
                                               'MAG_ISO':'F2_M_I',
                                               'MAGERR_ISO':'F2_MER_I',
                                               'FLUX_AUTO':'F2_F_A',
                                               'FLUXERR_AUTO':'F2_FER_A'},
                                               inplace=True)
            
            f2_df.append(f2_remove_99)
        elif filter_name == 'F3':

            f3_remove_99 = remove_99(df)

            f3_remove_99.rename(columns={'MAG_AUTO':'F3_M_A',
                                               'MAGERR_AUTO':'F3_MER_A',
                                               'MAG_ISO':'F3_M_I',
                                               'MAGERR_ISO':'F3_MER_I',
                                               'FLUX_AUTO':'F3_F_A',
                                               'FLUXERR_AUTO':'F3_FER_A'},
                                               inplace=True)
            
            f3_df.append(f3_remove_99)


        elif filter_name == 'F5':
        
            f5_remove_99 = remove_99(df)

            f5_remove_99.rename(columns={'MAG_AUTO':'F5_M_A',
                                               'MAGERR_AUTO':'F5_MER_A',
                                               'MAG_ISO':'F5_M_I',
                                               'MAGERR_ISO':'F5_MER_I',
                                               'FLUX_AUTO':'F5_F_A',
                                               'FLUXERR_AUTO':'F5_FER_A'},
                                               inplace=True)
            
            f5_df.append(f5_remove_99)

        elif filter_name == 'F7':
            
            f7_remove_99 = remove_99(df)

            f7_remove_99.rename(columns={'MAG_AUTO':'F7_M_A',
                                               'MAGERR_AUTO':'F7_MER_A',
                                               'MAG_ISO':'F7_M_I',
                                               'MAGERR_ISO':'F7_MER_I',
                                               'FLUX_AUTO':'F7_F_A',
                                               'FLUXERR_AUTO':'F7_FER_A'},
                                               inplace=True)
            
            f7_df.append(f7_remove_99)
        else:
            print('None error')

    print(len(f1_df)) #lengths of dataframes per filter FUV
    print(len(f2_df))
    print(len(f3_df))
    print(len(f5_df))
    print(len(f7_df))

    f1_conc = pd.concat(f1_df,ignore_index=True)
    # # print(f1_conc)
    f2_conc = pd.concat(f2_df,ignore_index=True)
    f3_conc = pd.concat(f3_df,ignore_index=True)
    f5_conc = pd.concat(f5_df,ignore_index=True)
    f7_conc = pd.concat(f7_df,ignore_index=True)

    # print(f1_conc.columns)

    f1tabledf = Table.from_pandas(f1_conc)
    f2tabledf = Table.from_pandas(f2_conc)
    f3tabledf = Table.from_pandas(f3_conc)
    f5tabledf = Table.from_pandas(f5_conc)
    f7tabledf = Table.from_pandas(f7_conc)

    
    
    
    
    
    
    ###### Now merge the files filterwise

# merging dist

    rad = Angle(3,u.arcsec)




    # print(f1_conc.columns)
    # f1_conc

    # f1_coord_matches = []

    # for f1_index in range(len(f1_conc)):

    #     f1_ra = f1_conc['ALPHA_J2000'].iloc[f1_index]
    #     f1_dec = f1_conc['DELTA_J2000'].iloc[f1_index]

    #     coord_f1_out = SkyCoord(ra=f1_ra,dec=f1_dec,unit='deg',frame='icrs') 

    #     print(f'Doing for coord: {coord_f1_out}')

    #     index_of_matched_coord_f1 = []

    #     for f1_index_inside in range(len(f1_conc)):
    #         ra_inside_f1 = f1_conc['ALPHA_J2000'].iloc[f1_index_inside]
    #         dec_inside_f1 = f1_conc['DELTA_J2000'].iloc[f1_index_inside]

    #         coord_f1_inside = SkyCoord(ra=ra_inside_f1,dec=dec_inside_f1,frame='icrs',unit='deg')

    #         separation_f1 = coord_f1_inside.separation(coord_f1_out)
    #         separation_f1_deg = separation_f1.to(u.arcsec)

    #         if separation_f1_deg < rad:

    #             if separation_f1_deg == 0:
    #                 pass
    #             else:
    #                 print(f'The separation is {separation_f1_deg}')

    #                 f1_coord_matches.append(coord_f1_inside)

    #                 if coord_f1_out in f1_coord_matches:
    #                     pass
    #                 else:
    #                     index_of_matched_coord_f1.append(f1_index_inside)
            

    #     f1_conc.drop(index_of_matched_coord_f1,axis='rows',inplace=True)
    #     f1_conc.reset_index(drop=True,inplace=True)
    #     print(len(f1_conc))

    
    # final_f1 = Table.from_pandas(f1_conc)

    # binary_f1 = fits.BinTableHDU(data=final_f1,header=header)
    # primary_f1 = fits.PrimaryHDU(data=pri_data,header=pri_header)

    # final_f1_df = fits.HDUList([primary_f1,binary_f1])
    # final_f1_df.writeto('/Users/swagat98/Documents/Combine_cat_2/f1_merged.fits',overwrite=True)




# dataframe_build_1()

def dataframe_build_2():

### make another folder of same contents
## check for proximity of one folder
# if match found
    # find for sources
    # if filter is same find the mean
    # if filter is differeent shift the positions
    path_to_same = '/Users/swagat98/Documents/Combine_cat_2/Same_field/FUV_same/'
    path_to_same_2 = '/Users/swagat98/Documents/Combine_cat_2/Same_field/FUV_same_2/'
    
    
    # copy

    # copy_cmd = f'cp -r {path_to_same} /Users/swagat98/Documents/Combine_cat_2/Same_field/FUV_same_2'
    # subprocess.run(copy_cmd,shell=True,check=True)
    # print('Copy folder created')


    fuv_files = glob(f'{path_to_same}*.fits')
    # print(fuv_files)
    

    dist = Angle(14,u.arcmin)

    for x in fuv_files:



        img_fuv = fits.open(x)

        header_fuv = img_fuv[1].header
        data_fuv = img_fuv[1].data

        filter_out = header_fuv['FILTER']

        out_tabledf = Table(data_fuv)
        outtabledf = out_tabledf.to_pandas()

        outdf = remove_99(outtabledf) # dataframe of one searching

        if filter_out == 'F1':
            outdf.rename(columns={'MAG_AUTO':'F1_M_A',  # mag auto
                                               'MAGERR_AUTO':'F1_MER_A', #mag error auto 
                                               'MAG_ISO':'F1_M_I',  # mag iso
                                               'MAGERR_ISO':'F1_MER_I', #mag err iso
                                               'FLUX_AUTO':'F1_F_A', #
                                               'FLUXERR_AUTO':'F1_FER_A'},inplace=True)
        elif filter_out == 'F2':
            outdf.rename(columns={'MAG_AUTO':'F2_M_A',
                                               'MAGERR_AUTO':'F2_MER_A',
                                               'MAG_ISO':'F2_M_I',
                                               'MAGERR_ISO':'F2_MER_I',
                                               'FLUX_AUTO':'F2_F_A',
                                               'FLUXERR_AUTO':'F2_FER_A'},
                                               inplace=True) 
        elif filter_out == 'F3':
            outdf.rename(columns={'MAG_AUTO':'F3_M_A',
                                               'MAGERR_AUTO':'F3_MER_A',
                                               'MAG_ISO':'F3_M_I',
                                               'MAGERR_ISO':'F3_MER_I',
                                               'FLUX_AUTO':'F3_F_A',
                                               'FLUXERR_AUTO':'F3_FER_A'},
                                               inplace=True)
        elif filter_out == 'F5':
            outdf.rename(columns={'MAG_AUTO':'F5_M_A',
                                               'MAGERR_AUTO':'F5_MER_A',
                                               'MAG_ISO':'F5_M_I',
                                               'MAGERR_ISO':'F5_MER_I',
                                               'FLUX_AUTO':'F5_F_A',
                                               'FLUXERR_AUTO':'F5_FER_A'},
                                               inplace=True)
            
        elif filter_out == 'F7':
            outdf.rename(columns={'MAG_AUTO':'F7_M_A',
                                               'MAGERR_AUTO':'F7_MER_A',
                                               'MAG_ISO':'F7_M_I',
                                               'MAGERR_ISO':'F7_MER_I',
                                               'FLUX_AUTO':'F7_F_A',
                                               'FLUXERR_AUTO':'F7_FER_A'},
                                               inplace=True)



        ra_pnt = header_fuv['RA_PNT']
        dec_pnt = header_fuv['DEC_PNT']

        coord_fuv = SkyCoord(ra=ra_pnt,dec=dec_pnt,frame='icrs',unit='deg')
        # print(coord_fuv)
        obj = x.split('/')[-1]

        print('\n')
        print('Object File which is searching: ',obj)

        try:
            x = f'cp -r {path_to_same_2}{obj} /Users/swagat98/Documents/Combine_cat_2/Same_field/N/' #copying the main files to test folder
            subprocess.run(x,shell=True,check=True)
            print('Copying main file')
        except:
            print('Matched file copied')

# remove the file so that while searching the program won't find the file in the second folder
        try:
            y = f'rm -r {path_to_same_2}{obj}'
            subprocess.run(y,shell=True,check=True)
            print('Deleting main file')
        except:
            print('Matched file deleted')


        fuv_2 = glob(f'{path_to_same_2}*.fits') # keep it inside the loop so that the list returns back to its original state 

        # print(fuv_2)
        # The file which is searching will not be there in this list, and therefore, in the next loop should
        # not be considered

        matched_dataframe_list = []
        matched_dataframe_list.append(outdf) #put the dataframe of the one searching
        print('Dataframe of the one searching: ',matched_dataframe_list)
        print(outdf.shape)



        for inside_fuv in range(len(fuv_2)):




            
            img_inside = fits.open(fuv_2[inside_fuv])
            # img_inside.info()

            primary_header = img_inside[0].header
            primary_data = img_inside[0].data

            header_inside_fuv = img_inside[1].header
            data_inside_fuv = img_inside[1].data

            ra_in_fuv = header_inside_fuv['RA_PNT']
            dec_in_fuv = header_inside_fuv['DEC_PNT']
            obj_inside = fuv_2[inside_fuv].split('/')[-1].split('/')[0]

            filter_in = header_inside_fuv['FILTER']

            coord_fuv_in = SkyCoord(ra=ra_in_fuv,dec=dec_in_fuv,frame='icrs',unit='deg')

            separation = coord_fuv_in.separation(coord_fuv)
            separation_in_deg = separation.to(u.arcmin)

            if separation_in_deg < dist:

                print('\t\tMatch found')
                print(f'Object found: {obj_inside}')

# delete the matched object so that it doesn't search again

                del_match = f'rm -r {path_to_same_2}{obj_inside}'
                subprocess.run(del_match,shell=True,check=True)
                print('Deleted matched object')






                print('Separation in deg',separation_in_deg)

                print('Computing the dataframe')

                dftable_inside = Table(data_inside_fuv)
                df_inside_pre = dftable_inside.to_pandas()

                df_inside = remove_99(df_inside_pre)

                if filter_in == 'F1':
                    df_inside.rename(columns={'MAG_AUTO':'F1_M_A',  # mag auto
                                               'MAGERR_AUTO':'F1_MER_A', #mag error auto 
                                               'MAG_ISO':'F1_M_I',  # mag iso
                                               'MAGERR_ISO':'F1_MER_I', #mag err iso
                                               'FLUX_AUTO':'F1_F_A', #
                                               'FLUXERR_AUTO':'F1_FER_A'},inplace=True)
                elif filter_in == 'F2':
                    df_inside.rename(columns={'MAG_AUTO':'F2_M_A',
                                               'MAGERR_AUTO':'F2_MER_A',
                                               'MAG_ISO':'F2_M_I',
                                               'MAGERR_ISO':'F2_MER_I',
                                               'FLUX_AUTO':'F2_F_A',
                                               'FLUXERR_AUTO':'F2_FER_A'},
                                               inplace=True)
                elif filter_in == 'F3':
                    df_inside.rename(columns={'MAG_AUTO':'F3_M_A',
                                               'MAGERR_AUTO':'F3_MER_A',
                                               'MAG_ISO':'F3_M_I',
                                               'MAGERR_ISO':'F3_MER_I',
                                               'FLUX_AUTO':'F3_F_A',
                                               'FLUXERR_AUTO':'F3_FER_A'},
                                               inplace=True)
                elif filter_in == 'F5':
                    df_inside.rename(columns={'MAG_AUTO':'F5_M_A',
                                               'MAGERR_AUTO':'F5_MER_A',
                                               'MAG_ISO':'F5_M_I',
                                               'MAGERR_ISO':'F5_MER_I',
                                               'FLUX_AUTO':'F5_F_A',
                                               'FLUXERR_AUTO':'F5_FER_A'},
                                               inplace=True)
                elif filter_in == 'F7':
                    df_inside.rename(columns={'MAG_AUTO':'F7_M_A',
                                               'MAGERR_AUTO':'F7_MER_A',
                                               'MAG_ISO':'F7_M_I',
                                               'MAGERR_ISO':'F7_MER_I',
                                               'FLUX_AUTO':'F7_F_A',
                                               'FLUXERR_AUTO':'F7_FER_A'},
                                               inplace=True)

                matched_dataframe_list.append(df_inside)
                print(df_inside.shape)

                # print(df_inside.shape)

                # dataframe of matched points

        # print('Matched dataframe list: ',matched_dataframe_list[1:])

        
        
        if len(matched_dataframe_list) == 1:
            print('Only main file present. No matches')
        elif len(matched_dataframe_list) > 1:
            matched_concat = pd.concat(matched_dataframe_list,ignore_index=True)
            print(matched_concat.columns)
            print(matched_concat.shape)

            tabledf = Table.from_pandas(matched_concat)


            binary = fits.BinTableHDU(data=tabledf,header=header_fuv)
            primary = fits.PrimaryHDU(data=primary_data,header=primary_header)
            final = fits.HDUList([primary,binary])
            final.writeto(f'/Users/swagat98/Documents/Combine_cat_2/Same_field/Catalog_combined/{obj}',overwrite=True)



        

# dataframe_build_2()



def merge_sources():
    
    path_combined = '/Users/swagat98/Documents/Combine_cat_2/Same_field/Catalog_combined/'

    files = glob(f'{path_combined}*.fits')



    path_merged = '/Users/swagat98/Documents/Combine_cat_2/Same_field/Merged_folder/'

    if os.path.exists(path_merged):
        print('Path exists')
    else:
        print('Creating folder')
        os.mkdir(path_merged)


    # print(files)

    # img = fits.open(files[0])
    # img.info()
    rad = Angle(3,u.arcsec)

    for x in range(len(files)):
        
        img = fits.open(files[x])
        img.info()

        primary_header= img[0].header
        primary_data = img[0].data

        header = img[1].header
        data = img[1].data

        obj = files[x].split('/')[-1]
        print(obj)

        tabledf = Table(data)
        df = tabledf.to_pandas()

        # m = df.iloc[5].to_frame().T

        print(df.columns)

        df_copy = df.copy() ##df copy

        # print(df.reset_index())
        df.reset_index(inplace=True)

        # index_l = [190,191,192]
        # df.drop(index_l,axis='rows',inplace=True)
        # df.reset_index(inplace=True)
        # print(df)
        

        ra_dec_out = SkyCoord(ra=df['ALPHA_J2000'],dec=df['DELTA_J2000'],unit='deg',frame='icrs')
        

#put magnitude list here
# put searching index here
    #unique liat it out

        list_of_matches = []

        coordinates_of_matches = []
        for i in df.index:
            
            index_of_matches = []


            try:
                ra_out = df['ALPHA_J2000'].loc[i]
                dec_out = df['DELTA_J2000'].loc[i]
                # index_of_matches.append(i)

            except KeyError as e:
                print('Error',e)
            out_cord = SkyCoord(ra=ra_out,dec=dec_out,frame='icrs',unit='deg')

            # out_cord = ra_dec_out[i]
            print(f'Doing for coord number: {i}')
            print(f'Doing for coord number {out_cord}')

            # print(out_cord.ra.deg)
            print('Finding the row')
            # print(df['F1_M_A'].loc[df['ALPHA_J2000'] == out_cord.ra]) #column having outcord
            # print(df.loc[df['ALPHA_J2000'] == out_cord.ra])

            #uncomment the above two lines

            # print(df)Z

     #248.73578264, 38.00153502       

            if out_cord in coordinates_of_matches:
                print('Pass')
                continue
            else: 
                for j in df.index:
                    
                    ra_in = df['ALPHA_J2000'].loc[j]
                    dec_in = df['DELTA_J2000'].loc[j]

                    in_coord = SkyCoord(ra=ra_in,dec=dec_in,frame='icrs',unit='deg')

                    separation = in_coord.separation(out_cord)
                    sep_in_deg = separation.to(u.arcsec)

                    if sep_in_deg < rad:
                        if sep_in_deg == 0:
                            pass
                        else:
                            print('\n\t\tMATCH FOUND')
                            print(f'Coordinate number: {j}')
                            print(f'Coordinate found: {in_coord}')
                            coordinates_of_matches.append(in_coord)

                            index_of_matches.append(j)



                
                print('Index of matches', index_of_matches)
                
                if len(index_of_matches) > 0:
                    drop_rows = df.iloc[index_of_matches].reset_index(drop=True)

                    print('Drop rows:\n',drop_rows)


                    mag_auto_f1_sum = 0
                    mag_auto_err_f1_sum = 0
                    mag_iso_f1_sum = 0
                    mag_err_iso_f1_sum = 0
                    flux_auto_f1_sum = 0
                    flux_err_auto_f1_sum = 0

                    mag_auto_f2_sum = 0
                    mag_auto_err_f2_sum = 0
                    mag_iso_f2_sum = 0
                    mag_err_iso_f2_sum = 0
                    flux_auto_f2_sum = 0
                    flux_err_auto_f2_sum = 0

                    mag_auto_f3_sum = 0
                    mag_auto_err_f3_sum = 0
                    mag_iso_f3_sum = 0
                    mag_err_iso_f3_sum = 0
                    flux_auto_f3_sum = 0
                    flux_err_auto_f3_sum = 0

                    mag_auto_f5_sum = 0
                    mag_auto_err_f5_sum = 0
                    mag_iso_f5_sum = 0
                    mag_err_iso_f5_sum = 0
                    flux_auto_f5_sum = 0
                    flux_err_auto_f5_sum = 0

                    mag_auto_f7_sum = 0
                    mag_auto_err_f7_sum = 0
                    mag_iso_f7_sum = 0
                    mag_err_iso_f7_sum = 0
                    flux_auto_f7_sum = 0
                    flux_err_auto_f7_sum = 0


                    count_f1 = 0
                    count_f2 = 0
                    count_f3 = 0
                    count_f5 = 0
                    count_f7 = 0
                    # print(drop_rows.loc[index_of_matches,'F1_M_A'])

                    for val in range(len(drop_rows)):
                        print('val: ',val)

                        # print('Rows: ',drop_rows.loc[val,['F1_M_A','F2_M_A','F7_M_A']])

                        # print('F7 column values',drop_rows['F7_M_A'].iloc[val] )
                    
                        # print(count_f7)                        

                        # mag_f1 = drop_rows['']

                        try:
                            if drop_rows['F1_M_A'].iloc[val] > 0:
                                mag = drop_rows['F1_M_A'].iloc[val]
                                mag_auto_f1_sum += mag
                                count_f1 += 1
                                # print('Drop rows mag F1',mag)
                                print('F1',drop_rows['F1_MER_A'].iloc[val])


                                mag_err = drop_rows['F1_MER_A'].iloc[val]
                                mag_auto_err_f1_sum += mag_err
                                print('F1 err:',drop_rows['F1_MER_A'].iloc[val])

                                mag_iso = drop_rows['F1_M_I'].iloc[val]
                                mag_iso_f1_sum += mag_iso

                                mag_iso_err = drop_rows['F1_MER_I'].iloc[val]
                                mag_err_iso_f1_sum += mag_iso_err

                                flux_auto = drop_rows['F1_F_A'].iloc[val]
                                flux_auto_f1_sum += flux_auto

                                flux_err_auto = drop_rows['F1_FER_A'].iloc[val]
                                flux_err_auto_f1_sum += flux_err_auto

                            else:
                                print('F1 nan value')
                            # print('F7',drop_rows['F7_M_A'].iloc[val])

                            # if drop_rows['F1_ME_A']
                        
                        except:
                            print('F1 not present')
                        
                        try:
                            if drop_rows['F2_M_A'].iloc[val]>0:
                                mag = drop_rows['F2_M_A'].iloc[val]
                                mag_auto_f2_sum += mag
                                count_f2 += 1
                                # print('Drop rows mag F2',mag)
                                print('F2',drop_rows['F2_M_A'].iloc[val])

                                mag_err = drop_rows['F2_MER_A'].iloc[val]
                                mag_auto_err_f2_sum += mag_err
                                print('F2 err:',drop_rows['F2_MER_A'].iloc[val])

                                mag_iso = drop_rows['F2_M_I'].iloc[val]
                                mag_iso_f2_sum += mag_iso

                                mag_iso_err = drop_rows['F2_MER_I'].iloc[val]
                                mag_err_iso_f2_sum += mag_iso_err

                                flux_auto = drop_rows['F2_F_A'].iloc[val]
                                flux_auto_f2_sum += flux_auto

                                flux_err_auto = drop_rows['F2_FER_A'].iloc[val]
                                flux_err_auto_f2_sum += flux_err_auto

                            else:
                                print('F2 nan value')

                        except:
                            print('F2 not present')

                        try:
                            if drop_rows['F3_M_A'].iloc[val]>0:
                                mag = drop_rows['F3_M_A'].iloc[val]
                                mag_auto_f3_sum += mag
                                count_f3 += 1
                                # print('Drop rows mag F3:',mag)

                                print('F3',drop_rows['F3_M_A'].iloc[val])

                                mag_err = drop_rows['F3_MER_A'].iloc[val]
                                mag_auto_err_f3_sum += mag_err
                                print('F3 err:', drop_rows['F3_MER_A'].iloc[val])

                                mag_iso = drop_rows['F3_M_I'].iloc[val]
                                mag_iso_f3_sum += mag_iso

                                mag_iso_err = drop_rows['F3_MER_I'].iloc[val]
                                mag_err_iso_f3_sum += mag_iso_err

                                flux_auto = drop_rows['F3_F_A'].iloc[val]
                                flux_auto_f3_sum += flux_auto

                                flux_err_auto = drop_rows['F3_FER_A'].iloc[val]
                                flux_err_auto_f3_sum += flux_err_auto

                            else:
                                print('F3 nan value')
                        except:
                            print('F3 not present')

                        try:
                            if drop_rows['F5_M_A'].iloc[val]>0:
                                mag = drop_rows['F5_M_A'].iloc[val]
                                mag_auto_f5_sum += mag
                                count_f5 += 1

                                print('F5',drop_rows['F5_M_A'].iloc[val])


                                mag_err = drop_rows['F5_MER_A'].iloc[val]
                                mag_auto_err_f5_sum += mag_err
                                print('F5 err:',drop_rows['F5_MER_A'].iloc[val])

                                mag_iso = drop_rows['F5_M_I'].iloc[val]
                                mag_iso_f5_sum += mag_iso

                                mag_iso_err = drop_rows['F5_MER_I'].iloc[val]
                                mag_err_iso_f5_sum += mag_iso_err

                                flux_auto = drop_rows['F5_F_A'].iloc[val]
                                flux_auto_f5_sum += flux_auto

                                flux_err_auto = drop_rows['F5_FER_A'].iloc[val]
                                flux_err_auto_f5_sum += flux_err_auto


                            else:
                                print('F5 nan value')

                        except:
                            print('F5 not present')

                        try:
                            if drop_rows['F7_M_A'].iloc[val] > 0:
                                mag = drop_rows['F7_M_A'].iloc[val]
                                mag_auto_f7_sum += mag
                                count_f7 += 1
                                # print('Drop rows mag F7: ',mag)

                                print('F7',drop_rows['F7_M_A'].iloc[val])

                                mag_err = drop_rows['F7_MER_A'].iloc[val]
                                mag_auto_err_f7_sum += mag_err
                                print('F7 err:',drop_rows['F7_MER_A'].iloc[val])

                                mag_iso = drop_rows['F7_M_I'].iloc[val]
                                mag_iso_f7_sum += mag_iso

                                mag_iso_err = drop_rows['F7_MER_I'].iloc[val]
                                mag_err_iso_f7_sum += mag_iso_err

                                flux_auto = drop_rows['F7_F_A'].iloc[val]
                                flux_auto_f7_sum += flux_auto

                                flux_err_auto = drop_rows['F7_FER_A'].iloc[val]
                                flux_err_auto_f7_sum += flux_err_auto

                            else:
                                print('F7 nan value')
                        # print('F7 Loc val',df.loc[val,'F7_M_A'])
                        except:
                            print('F7 not present')


                        





                    print('mean_mag_sum_f1:', mag_auto_f1_sum)
                    print('mean_mag_sum_f3:', mag_auto_f3_sum)
                    print('mean_mag_sum_f2:', mag_auto_f2_sum)
                    print('mean_mag_sum_f5:', mag_auto_f5_sum)
                    print('mean_mag_sum_f7:', mag_auto_f7_sum)

                    print('Count F1:', count_f1)
                    print('Count F2:', count_f2)
                    print('Count F3:', count_f2)
                    print('Count F5:', count_f2)
                    print('Count F7:', count_f7)
                    # print('Count F7:', count_f7)
                    # print('Count F7:', count_f7)


                    try:
                        magnitude_replace = mag_auto_f1_sum/count_f1
                        df['F1_M_A'].loc[df['ALPHA_J2000'] == out_cord.ra] = magnitude_replace
                        # df.loc[]
                        # print(df['F1_M_A'].loc[df['ALPHA_J2000'] == out_cord.ra])


                        mag_auto_err_replace = mag_auto_err_f1_sum/count_f1
                        df['F1_MER_A'].loc[df['ALPHA_J2000'] == out_cord.ra] = mag_auto_err_replace
                        # print(df['F1_MER_A'].loc[df['ALPHA_J2000']==out_cord.ra])

                        mag_iso_replace = mag_iso_f1_sum/count_f1
                        df['F1_M_I'].loc[df['ALPHA_J2000'] == out_cord.ra] = mag_iso_replace

                        mag_err_iso_replace = mag_err_iso_f1_sum / count_f1
                        df['F1_MER_I'].loc[df['ALPHA_J2000'] == out_cord.ra] = mag_err_iso_replace

                        flux_auto_replace = flux_auto_f1_sum / count_f1
                        df['F1_F_A'].loc[df['ALPHA_J2000'] == out_cord.ra] = flux_auto_replace

                        flux_err_replace = flux_err_auto_f1_sum / count_f1
                        df['F1_FER_A'].loc[df['ALPHA_J2000'] == out_cord.ra] = flux_err_replace

                    except:
                        print('F1 data not present')

                    try:
                        magnitude_replace = mag_auto_f2_sum/count_f2

                        df['F2_M_A'].loc[df['ALPHA_J2000'] == out_cord.ra] = magnitude_replace

                        mag_auto_err_replace = mag_auto_err_f2_sum/count_f2
                        df['F2_MER_A'].loc[df['ALPHA_J2000'] == out_cord.ra] = mag_auto_err_replace

                        mag_iso_replace = mag_iso_f2_sum / count_f2
                        df['F2_M_I'].loc[df['ALPHA_J2000'] == out_cord.ra] = mag_iso_replace

                        mag_err_iso_replace = mag_err_iso_f2_sum / count_f2
                        df['F2_MER_I'].loc[df['ALPHA_J2000'] == out_cord.ra] = mag_err_iso_replace

                        flux_auto_replace = flux_auto_f2_sum / count_f2
                        df['F2_F_A'].loc[df['ALPHA_J2000'] == out_cord.ra] = flux_auto_replace

                        flux_err_replace = flux_err_auto_f2_sum / count_f2
                        df['F2_FER_A'].loc[df['ALPHA_J2000'] == out_cord.ra] = flux_err_replace 


                    except:
                        print('F2 data not present')
                    
                    try:
                        magnitude_replace = mag_auto_f3_sum/count_f3
                        df['F3_M_A'].loc[df['ALPHA_J2000'] == out_cord.ra] = magnitude_replace

                        mag_auto_err_replace = mag_auto_err_f3_sum / count_f3
                        df['F3_MER_A'].loc[df['ALPHA_J2000'] == out_cord.ra] = mag_auto_err_replace

                        mag_iso_replace = mag_iso_f3_sum / count_f3
                        df['F3_M_I'].loc[df['ALPHA_J2000'] == out_cord.ra] = mag_iso_replace

                        mag_err_iso_replace = mag_err_iso_f3_sum / count_f3
                        df['F3_MER_I'].loc[df['ALPHA_J2000'] == out_cord.ra]  = mag_err_iso_replace

                        flux_auto_replace = flux_auto_f3_sum / count_f3
                        df['F3_F_A'].loc[df['ALPHA_J2000'] == out_cord.ra] = flux_auto_replace

                        flux_err_replace = flux_err_auto_f3_sum / count_f3
                        df['F3_FER_A'].loc[df['ALPHA_J2000'] == out_cord.ra] = flux_err_replace

                    except:
                        print('F3 data not present')
                    
                    try:
                        magnitude_replace = mag_auto_f5_sum/count_f5
                        df['F5_M_A'].loc[df['ALPHA_J2000']==out_cord.ra] = magnitude_replace

                        mag_auto_err_replace = mag_auto_err_f5_sum / count_f5
                        df['F5_MER_A'].loc[df['ALPHA_J2000'] == out_cord.ra] = mag_auto_err_replace

                        mag_iso_replace = mag_iso_f5_sum / count_f5
                        df['F5_M_I'].loc[df['ALPHA_J2000'] == out_cord.ra] = mag_iso_replace

                        mag_err_iso_replace = mag_err_iso_f5_sum / count_f5
                        df['F5_MER_I'].loc[df['ALPHA_J2000'] == out_cord.ra] = mag_err_iso_replace

                        flux_auto_replace = flux_auto_f5_sum / count_f5
                        df['F5_F_A'].loc[df['ALPHA_J2000'] == out_cord.ra] = flux_auto_replace

                        flux_err_replace = flux_err_auto_f5_sum / count_f5
                        df['F5_FER_A'].loc[df['ALPHA_J2000'] == out_cord.ra]  = flux_err_replace

                    except:
                        print('F5 data not present')



                    try:
                        magnitude_replace=mag_auto_f7_sum/count_f7
                        df['F7_M_A'].loc[df['ALPHA_J2000'] == out_cord.ra] = magnitude_replace

                        mag_auto_err_replace = mag_auto_err_f7_sum / count_f7
                        df['F7_MER_A'].loc[df['ALPHA_J2000'] == out_cord.ra] = mag_auto_err_replace

                        mag_iso_replace = mag_iso_f7_sum / count_f7
                        df['F7_M_I'].loc[df['ALPHA_J2000'] == out_cord.ra] = mag_iso_replace

                        mag_err_iso_replace = mag_err_iso_f7_sum / count_f7
                        df['F7_MER_I'].loc[df['ALPHA_J2000'] == out_cord.ra] = mag_err_iso_replace

                        flux_auto_replace = flux_auto_f7_sum / count_f7
                        df['F7_F_A'].loc[df['ALPHA_J2000'] == out_cord.ra] = flux_auto_replace

                        flux_err_replace = flux_err_auto_f7_sum / count_f7
                        df['F7_FER_A'].loc[df['ALPHA_J2000'] == out_cord.ra] = flux_err_replace

                    except:
                        print('F7 data not present')

                    
                        print('\n')

                df.drop(index_of_matches,inplace=True,axis='rows')
                df.reset_index(drop=True,inplace=True)
                
                try:
                    df.drop('index',inplace=True)
                except:
                    print('Index not present')
                
                    
     
                # print('shape',df.shape) # after dropping the matched location, the shape
                
                print('\n\n')
                

            




                # df.drop(index_list,axis='rows',inplace=True)
                # df.reset_index(drop=True,inplace=True)

        tabledf = Table.from_pandas(df)              

        binary = fits.BinTableHDU(data=tabledf,header=header)
        primary = fits.PrimaryHDU(header=primary_header)

        final = fits.HDUList([primary,binary])
        final.writeto(f'{path_merged}{obj}',overwrite=True)
        




















                #     # print(df_copy)
                #     print('\n')
# initial df doesn't have index


# apply next filter by searching for len(index) >1 or ==1

                # print('Count')

                # try:
                #     print('Mean magnitude count',mean_num_f1/count_f1)

                #     f1_mag_auto = mean_num_f1/count_f1 
                #     drop_rows.loc[0,'F1_M_A'] = f1_mag_auto
                # except:
                #     print('No f1 values')

                
                # print(drop_rows)
            # replace_rows = drop_rows.iloc[0].to_frame().T

                # final = pd.concat([df,replace_rows],ignore_index=True)
                # print(final.shape)

        # binary = fits.BinTableHDU(data=df)
        # primary = fits.PrimaryHDU(data=primary_data,header=primary_header)

        # final = fits.HDUList([primary,binary])
        # final.writeto(f'{path_merged}{obj}',overwrite=True)

                        



                





        



# merge_sources()




def dataframe_build3():

    nuv_path = '/Users/swagat98/Documents/Combine_cat_2/Same_field/NUV_same/'
    test_folder = '/Users/swagat98/Documents/Combine_cat_2/test/'
    test_folder_2 = '/Users/swagat98/Documents/Combine_cat_2/Same_field/FUV_same/'

    nuv_files = glob(f'{nuv_path}*.fits')
    fuv_files = glob(f'{test_folder}*.fits')
    fuv_files_2 = glob(f'{test_folder_2}*.fits')
    
    dist = Angle(14,u.arcmin)


    nuv_table = []
    fuv_table = []
    
    obj_list = []


    for nuv_len in nuv_files:

        nuv = fits.open(nuv_len)
        nuv.info()

        primary_hdr = nuv[0].header
        nuv_data = nuv[1].data
        nuv_header= nuv[1].header

        nuv_filter = nuv_header['FILTER']

        nuvtabledf = Table(nuv_data)
        nuvdf = nuvtabledf.to_pandas()

        print(nuvdf.columns)

        if nuv_filter == 'F1':
            nuvdf.rename(columns={'MAG_AUTO':'N1_M_A',
                                  'MAGERR_AUTO':'N1_MER_A',
                                  'MAG_ISO':'N1_M_I',
                                  'MAGERR_ISO':'N1_MER_I',
                                  'FLUX_AUTO':'N1_F_A',
                                  'FLUXERR_AUTO':'N1_FER_A'},inplace=True)
        if nuv_filter == 'F2':
            nuvdf.rename(columns={'MAG_AUTO':'N2_M_A',
                                  'MAGERR_AUTO':'N2_MER_A',
                                  'MAG_ISO':'N2_M_I',
                                  'MAGERR_ISO':'N2_MER_I',
                                  'FLUX_AUTO':'N2_F_A',
                                  'FLUXERR_AUTO':'N2_FER_A'},inplace=True)

            
        if nuv_filter == 'F3':
            nuvdf.rename(columns={'MAG_AUTO':'N3_M_A',
                                  'MAGERR_AUTO':'N3_MER_A',
                                  'MAG_ISO':'N3_M_I',
                                  'MAGERR_ISO':'N3_MER_I',
                                  'FLUX_AUTO':'N3_F_A',
                                  'FLUXERR_AUTO':'N3_FER_A'},inplace=True)
            
        if nuv_filter == 'F5':
            nuvdf.rename(columns={'MAG_AUTO':'N5_M_A',
                                  'MAGERR_AUTO':'N5_MER_A',
                                  'MAG_ISO':'N5_M_I',
                                  'MAGERR_ISO':'N5_MER_I',
                                  'FLUX_AUTO':'N5_F_A',
                                  'FLUXERR_AUTO':'N5_FER_A'},inplace=True)
            
        if nuv_filter == 'F6':
            nuvdf.rename(columns={'MAG_AUTO':'N6_M_A',
                                  'MAGERR_AUTO':'N6_MER_A',
                                  'MAG_ISO':'N6_M_I',
                                  'MAGERR_ISO':'N6_MER_I',
                                  'FLUX_AUTO':'N6_F_A',
                                  'FLUXERR_AUTO':'N6_FER_A'},inplace=True)
            

        
        ra_nuv = nuv_header['RA_PNT']
        dec_nuv = nuv_header['DEC_PNT']
        
        
        coord_nuv = SkyCoord(ra=ra_nuv,dec=dec_nuv,frame='icrs',unit='deg')

        for fuv_len in fuv_files:

            fuv = fits.open(fuv_len)      

            fuv_data = fuv[1].data
            fuv_hdr = fuv[1].header

            fuv_filter = fuv_hdr['FILTER']

            fuvtabledf = Table(fuv_data)
            fuvdf = fuvtabledf.to_pandas()
            
            obj_fuv = fuv_len.split('/')[-1]

            if fuv_filter == 'F1':
                fuvdf.rename(columns={'MAG_AUTO':'F1_M_A',
                                  'MAGERR_AUTO':'F1_MER_A',
                                  'MAG_ISO':'F1_M_I',
                                  'MAGERR_ISO':'F1_MER_I',
                                  'FLUX_AUTO':'F1_F_A',
                                  'FLUXERR_AUTO':'F1_FER_A'},inplace=True)

                
            if fuv_filter == 'F2':
                fuvdf.rename(columns={'MAG_AUTO':'F2_M_A',
                                  'MAGERR_AUTO':'F2_MER_A',
                                  'MAG_ISO':'F2_M_I',
                                  'MAGERR_ISO':'F2_MER_I',
                                  'FLUX_AUTO':'F2_F_A',
                                  'FLUXERR_AUTO':'F2_FER_A'},inplace=True)


            if fuv_filter == 'F3':
                fuvdf.rename(columns={'MAG_AUTO':'F3_M_A',
                                  'MAGERR_AUTO':'F3_MER_A',
                                  'MAG_ISO':'F3_M_I',
                                  'MAGERR_ISO':'F3_MER_I',
                                  'FLUX_AUTO':'F3_F_A',
                                  'FLUXERR_AUTO':'F3_FER_A'},inplace=True)
                

            if fuv_filter == 'F5':
                fuvdf.rename(columns={'MAG_AUTO':'F5_M_A',
                                  'MAGERR_AUTO':'F5_MER_A',
                                  'MAG_ISO':'F5_M_I',
                                  'MAGERR_ISO':'F5_MER_I',
                                  'FLUX_AUTO':'F5_F_A',
                                  'FLUXERR_AUTO':'F5_FER_A'},inplace=True)
                

            if fuv_filter == 'F7':
                fuvdf.rename(columns={'MAG_AUTO':'F7_M_A',
                                  'MAGERR_AUTO':'F7_MER_A',
                                  'MAG_ISO':'F7_M_I',
                                  'MAGERR_ISO':'F7_MER_I',
                                  'FLUX_AUTO':'F7_F_A',
                                  'FLUXERR_AUTO':'F7_FER_A'},inplace=True)

                
            
            
            ra_fuv = fuv_hdr['RA_PNT']
            dec_fuv = fuv_hdr['DEC_PNT']
            
            coord_fuv = SkyCoord(ra=ra_fuv,dec=dec_fuv,frame='icrs',unit='deg')
            
            sep = coord_fuv.separation(coord_nuv)
            sep_in_deg = sep.to(u.arcmin)
            
            if sep_in_deg < dist:
                
                print('\t\tMATCH FOUND')
                print('\tObject found: ',obj_fuv)
                
                
                fuv_table.append(fuvdf)
                nuv_table.append(nuvdf)
                
                obj_list.append(obj_fuv)
                
                
                # copy to FUV_after_NUV_check
                # remove from
                
                copy_cmd = f'cp -r {test_folder}{obj_fuv} /Users/swagat98/Documents/Combine_cat_2/FUV_after_NUV_check/'
                subprocess.run(copy_cmd,shell=True,check=True)
                print('Files copied')
                
                
    print(len(fuv_table))
    print(len(nuv_table))
    
    nuv_check_path = '/Users/swagat98/Documents/Combine_cat_2/FUV_after_NUV_check/'
    files_in_nuv_check_folder = glob(f'{nuv_check_path}*.fits')
    
    for t in files_in_nuv_check_folder:
        
        obj_names = t.split('/')[-1]
        print(obj_names)
        
        del_cmd = f'rm -r /Users/swagat98/Documents/Combine_cat_2/test/{obj_names}'
        subprocess.run(del_cmd,shell=True,check=True)
        print('Files deleted')

    
#     combined_dataframe_list = []
    
#     for m in range(len(fuv_table)):
        
#         fuv_t = fuv_table[m]
#         nuv_t = nuv_table[m]
        
# #         print(fuv_t.shape)
# #         print(nuv_t.shape)
        
#         concat_df = pd.concat([nuv_t,fuv_t])
        
#         combined_dataframe_list.append(concat_df)
        
#     print(combined_dataframe_list[0].columns)
    
    
    rad = Angle(6,u.arcsec)
    
    for m in range(len(fuv_table)):
        
        fuv_t = fuv_table[m]
        
        nuv_t = nuv_table[m]
        
        
        
        nuv_t['F1_M_A'] = -999
        nuv_t['F1_MER_A'] = -999
        nuv_t['F1_M_I']=-999
        nuv_t['F1_MER_I']=-999
        nuv_t['F1_F_A']=-999
        nuv_t['F1_FER_A']=-999
        
        nuv_t['F2_M_A'] = -999
        nuv_t['F2_MER_A'] = -999
        nuv_t['F2_M_I']=-999
        nuv_t['F2_MER_I']=-999
        nuv_t['F2_F_A']=-999
        nuv_t['F2_FER_A']=-999
        
        nuv_t['F3_M_A'] = -999
        nuv_t['F3_MER_A'] = -999
        nuv_t['F3_M_I']=-999
        nuv_t['F3_MER_I']=-999
        nuv_t['F3_F_A']=-999
        nuv_t['F3_FER_A']=-999
        
        nuv_t['F5_M_A'] = -999
        nuv_t['F5_MER_A'] = -999
        nuv_t['F5_M_I']=-999
        nuv_t['F5_MER_I']=-999
        nuv_t['F5_F_A']=-999
        nuv_t['F5_FER_A']=-999
        
        nuv_t['F7_M_A'] = -999
        nuv_t['F7_MER_A'] = -999
        nuv_t['F7_M_I']=-999
        nuv_t['F7_MER_I']=-999
        nuv_t['F7_F_A']=-999
        nuv_t['F7_FER_A']=-999
        
        
        print(nuv_t)
        
        
        
        
        
        nuv_coord = SkyCoord(ra=nuv_t['ALPHA_J2000'],dec=nuv_t['DELTA_J2000'],frame='icrs',unit='deg')
        
        for count_0 in range(len(nuv_t)):
            
#             print(nuv_coord[count_0])
            print('Doing for coord number: ',count_0)
    
    
            index_to_be_deleted = []
        
            for count_1 in range(len(fuv_t)):
                
                ra_in_fuv = fuv_t.loc[count_1,'ALPHA_J2000']
                dec_in_fuv = fuv_t.loc[count_1,'DELTA_J2000']
                
                fuv_coord = SkyCoord(ra=ra_in_fuv,dec=dec_in_fuv,frame='icrs',unit='deg')
                
                sep_t = nuv_coord[count_0].separation(fuv_coord)
                sep_t_deg = sep_t.to(u.arcsec)
                
#                 print(sep_t_deg)
                if sep_t_deg < rad:
                    
                    
                    print('\t\tMATCH FOUND')
                    print(f'Coordinate match found: {fuv_coord}')
                    print(f'Index found: {count_1}')
                    print(f'Fuv index: {count_0}')
                    
                    # below will only work if there is a 
                    
                    try:
                        print(f'Magnitude Before: {nuv_t["F1_M_A"].iloc[count_0]}')

                        nuv_t['F1_M_A'].iloc[count_0] = fuv_t['F1_M_A'].iloc[count_1]
                        nuv_t['F1_MER_A'].iloc[count_0] = fuv_t['F1_MER_A'].iloc[count_1]
                        nuv_t['F1_M_I'].iloc[count_0] = fuv_t['F1_M_I'].iloc[count_1]
                        nuv_t['F1_MER_I'].iloc[count_0] = fuv_t['F1_MER_I'].iloc[count_1]
                        nuv_t['F1_F_A'].iloc[count_0] = fuv_t['F1_F_A'].iloc[count_1]
                        nuv_t['F1_FER_A'].iloc[count_0] = fuv_t['F1_FER_A'].iloc[count_1]


                        print(f'Magnitude After: {nuv_t["F1_M_A"].iloc[count_0]}')

                        index_to_be_deleted.append(count_1)

                        print('\n\n')
                    except:
                        print('F1 not present')
                        print('\n')
                        

                    try:
                        print(f'Magnitude Before: {nuv_t["F2_M_A"].iloc[count_0]}')

                        nuv_t['F2_M_A'].iloc[count_0] = fuv_t['F2_M_A'].iloc[count_1]
                        nuv_t['F2_MER_A'].iloc[count_0] = fuv_t['F2_MER_A'].iloc[count_1]
                        nuv_t['F2_M_I'].iloc[count_0] = fuv_t['F2_M_I'].iloc[count_1]
                        nuv_t['F2_MER_I'].iloc[count_0] = fuv_t['F2_MER_I'].iloc[count_1]
                        nuv_t['F2_F_A'].iloc[count_0] = fuv_t['F2_F_A'].iloc[count_1]
                        nuv_t['F2_FER_A'].iloc[count_0] = fuv_t['F2_FER_A'].iloc[count_1]


                        print(f'Magnitude After: {nuv_t["F2_M_A"].iloc[count_0]}')

                        index_to_be_deleted.append(count_1)

                        print('\n\n')
                    except:
                        print('F2 not present')
                        print('\n')
                    

                    try:
                        print(f'Magnitude Before: {nuv_t["F3_M_A"].iloc[count_0]}')

                        nuv_t['F3_M_A'].iloc[count_0] = fuv_t['F3_M_A'].iloc[count_1]
                        nuv_t['F3_MER_A'].iloc[count_0] = fuv_t['F3_MER_A'].iloc[count_1]
                        nuv_t['F3_M_I'].iloc[count_0] = fuv_t['F3_M_I'].iloc[count_1]
                        nuv_t['F3_MER_I'].iloc[count_0] = fuv_t['F3_MER_I'].iloc[count_1]
                        nuv_t['F3_F_A'].iloc[count_0] = fuv_t['F3_F_A'].iloc[count_1]
                        nuv_t['F3_FER_A'].iloc[count_0] = fuv_t['F3_FER_A'].iloc[count_1]


                        print(f'Magnitude After: {nuv_t["F3_M_A"].iloc[count_0]}')

                        index_to_be_deleted.append(count_1)

                        print('\n\n')
                    except:
                        print('F3 not present')
                        print('\n')
                        
                        

                    try:
                        print(f'Magnitude Before: {nuv_t["F5_M_A"].iloc[count_0]}')

                        nuv_t['F5_M_A'].iloc[count_0] = fuv_t['F5_M_A'].iloc[count_1]
                        nuv_t['F5_MER_A'].iloc[count_0] = fuv_t['F5_MER_A'].iloc[count_1]
                        nuv_t['F5_M_I'].iloc[count_0] = fuv_t['F5_M_I'].iloc[count_1]
                        nuv_t['F5_MER_I'].iloc[count_0] = fuv_t['F5_MER_I'].iloc[count_1]
                        nuv_t['F5_F_A'].iloc[count_0] = fuv_t['F5_F_A'].iloc[count_1]
                        nuv_t['F5_FER_A'].iloc[count_0] = fuv_t['F5_FER_A'].iloc[count_1]


                        print(f'Magnitude After: {nuv_t["F5_M_A"].iloc[count_0]}')

                        index_to_be_deleted.append(count_1)

                        print('\n\n')
                    except:
                        print('F5 not present')
                        print('\n')
                        
                        

                    try:
                        print(f'Magnitude Before: {nuv_t["F7_M_A"].iloc[count_0]}')

                        nuv_t['F7_M_A'].iloc[count_0] = fuv_t['F7_M_A'].iloc[count_1]
                        nuv_t['F7_MER_A'].iloc[count_0] = fuv_t['F7_MER_A'].iloc[count_1]
                        nuv_t['F7_M_I'].iloc[count_0] = fuv_t['F7_M_I'].iloc[count_1]
                        nuv_t['F7_MER_I'].iloc[count_0] = fuv_t['F7_MER_I'].iloc[count_1]
                        nuv_t['F7_F_A'].iloc[count_0] = fuv_t['F7_F_A'].iloc[count_1]
                        nuv_t['F7_FER_A'].iloc[count_0] = fuv_t['F7_FER_A'].iloc[count_1]


                        print(f'Magnitude After: {nuv_t["F7_M_A"].iloc[count_0]}')

                        index_to_be_deleted.append(count_1)

                        print('\n\n')
                    except:
                        print('F7 not present')
                        print('\n')
                        
                        
                        
                        
            
            print('Shape before drop axis: ',fuv_t.shape)
            fuv_t.drop(index_to_be_deleted,inplace=True,axis='rows')
            fuv_t.reset_index(drop=True,inplace=True)
            print('Shape after drop axis:',fuv_t.shape)
            
        #concat the remaining fuv_t and nuv_t
            
        
        final = pd.concat([nuv_t,fuv_t])        
        print('final shape: ',final.shape)
        
        binarytable = Table.from_pandas(final)
        
        binary = fits.BinTableHDU(data=binarytable)
        primary = fits.PrimaryHDU()
        
        final_hdu = fits.HDUList([primary,binary])
        final_hdu.writeto(f'/Users/swagat98/Documents/Combine_cat_2/Same_field/NUV_output/{obj_list[m]}.fits',overwrite=True)
    
#     del_cmd = f'rm -r {test_folder}{obj_fuv}'
#     subprocess.run(del_cmd,shell=True,check=True)
#     print('Files deleted')
        

# dataframe_build3()


def dataframe_build4():

    path_to_same_nuv = '/Users/swagat98/Documents/Combine_cat_2/Same_field/NUV_same/'
    path_to_merged_fuv = '/Users/swagat98/Documents/Combine_cat_2/Same_field/Merged_folder/'

    nuv = glob(f'{path_to_same_nuv}*.fits')
    fuv_merged = glob(f'{path_to_merged_fuv}*.fits')


    dist = Angle(14,u.arcmin)

    fuv_table = []
    nuv_table = []
    obj_list = []

    print('Copying NUV data to FUV merged...')
    print('\n')

    for x in nuv:
        
        nuv_img = fits.open(x)
        
        nuv_hdr = nuv_img[1].header
        nuv_data = nuv_img[1].data

        ra_nuv = nuv_hdr['RA_PNT']
        dec_nuv= nuv_hdr['DEC_PNT']

        nuv_filter = nuv_hdr['FILTER']

        coord_nuv = SkyCoord(ra=ra_nuv,dec=dec_nuv,frame='icrs',unit='deg')

        obj = x.split('/')[-1]
        print(f'Doing for object (NUV): {obj}')


        nuvtabledf = Table(nuv_data)
        nuvdf = nuvtabledf.to_pandas()

        if nuv_filter == 'F1':
            nuvdf.rename(columns={'MAG_AUTO':'N1_M_A',
                                  'MAGERR_AUTO':'N1_MER_A',
                                  'MAG_ISO':'N1_M_I',
                                  'MAGERR_ISO':'N1_MER_I',
                                  'FLUX_AUTO':'N1_F_A',
                                  'FLUXERR_AUTO':'N1_FER_A'},inplace=True)
        if nuv_filter == 'F2':
            nuvdf.rename(columns={'MAG_AUTO':'N2_M_A',
                                  'MAGERR_AUTO':'N2_MER_A',
                                  'MAG_ISO':'N2_M_I',
                                  'MAGERR_ISO':'N2_MER_I',
                                  'FLUX_AUTO':'N2_F_A',
                                  'FLUXERR_AUTO':'N2_FER_A'},inplace=True)

            
        if nuv_filter == 'F3':
            nuvdf.rename(columns={'MAG_AUTO':'N3_M_A',
                                  'MAGERR_AUTO':'N3_MER_A',
                                  'MAG_ISO':'N3_M_I',
                                  'MAGERR_ISO':'N3_MER_I',
                                  'FLUX_AUTO':'N3_F_A',
                                  'FLUXERR_AUTO':'N3_FER_A'},inplace=True)
            
        if nuv_filter == 'F5':
            nuvdf.rename(columns={'MAG_AUTO':'N5_M_A',
                                  'MAGERR_AUTO':'N5_MER_A',
                                  'MAG_ISO':'N5_M_I',
                                  'MAGERR_ISO':'N5_MER_I',
                                  'FLUX_AUTO':'N5_F_A',
                                  'FLUXERR_AUTO':'N5_FER_A'},inplace=True)
            
        if nuv_filter == 'F6':
            nuvdf.rename(columns={'MAG_AUTO':'N6_M_A',
                                  'MAGERR_AUTO':'N6_MER_A',
                                  'MAG_ISO':'N6_M_I',
                                  'MAGERR_ISO':'N6_MER_I',
                                  'FLUX_AUTO':'N6_F_A',
                                  'FLUXERR_AUTO':'N6_FER_A'},inplace=True)
            


        for y in fuv_merged:
            fuv_img = fits.open(y)

            fuv_hdr = fuv_img[1].header
            fuv_data = fuv_img[1].data

            ra_fuv = fuv_hdr['RA_PNT']
            dec_fuv= fuv_hdr['DEC_PNT']

            coord_fuv = SkyCoord(ra=ra_fuv,dec=dec_fuv,frame='icrs',unit='deg')

            sep = coord_nuv.separation(coord_fuv)
            sep_in_deg = sep.to(u.arcmin) 

            obj_fuv = y.split('/')[-1]

            fuvtabledf = Table(fuv_data)
            fuvdf = fuvtabledf.to_pandas()


            if sep_in_deg < dist:

                print('\t\tMATCH FOUND')

                print(f'Match found in object (FUV) : {obj_fuv}')
                print('\n')

                obj_list.append(obj_fuv) #append obj name

                fuv_table.append(fuvdf)
                nuv_table.append(nuvdf)

    # print(fuv_table[0])
    # print(nuv_table[0])

    rad = Angle(1.5,u.arcsec)


    
    for m in range(len(fuv_table)):

        fuv_t = fuv_table[m]
        nuv_t = nuv_table[m]

        fuv_t['N1_M_A'] = -999
        fuv_t['N1_MER_A'] = -999
        fuv_t['N1_M_I'] = -999
        fuv_t['N1_MER_I'] = -999
        fuv_t['N1_F_A'] = -999
        fuv_t['N1_FER_A'] = -999

        fuv_t['N2_M_A'] = -999
        fuv_t['N2_MER_A'] = -999
        fuv_t['N2_M_I'] = -999
        fuv_t['N2_MER_I'] = -999
        fuv_t['N2_F_A'] = -999
        fuv_t['N2_FER_A'] = -999

        fuv_t['N3_M_A'] = -999
        fuv_t['N3_MER_A'] = -999
        fuv_t['N3_M_I'] = -999
        fuv_t['N3_MER_I'] = -999
        fuv_t['N3_F_A'] = -999
        fuv_t['N3_FER_A'] = -999

        fuv_t['N5_M_A'] = -999
        fuv_t['N5_MER_A'] = -999
        fuv_t['N5_M_I'] = -999
        fuv_t['N5_MER_I'] = -999
        fuv_t['N5_F_A'] = -999
        fuv_t['N5_FER_A'] = -999

        fuv_t['N6_M_A'] = -999
        fuv_t['N6_MER_A'] = -999
        fuv_t['N6_M_I'] = -999
        fuv_t['N6_MER_I'] = -999
        fuv_t['N6_F_A'] = -999
        fuv_t['N6_FER_A'] = -999

        print(fuv_t.columns)
        for count0 in range(len(fuv_t)):

            ra_f = fuv_t['ALPHA_J2000'].iloc[count0]
            dec_f = fuv_t['DELTA_J2000'].iloc[count0]

            fuvcoord = SkyCoord(ra=ra_f,dec=dec_f,frame='icrs',unit='deg')

            print(f'Doing for coord : {fuvcoord}')

            index_of_matches = []

            for count1 in range(len(nuv_t)):
                
                # print('Doing for number',count1)# 
                ra_n = nuv_t['ALPHA_J2000'].iloc[count1]
                dec_n = nuv_t['DELTA_J2000'].iloc[count1]

                nuvcoord = SkyCoord(ra=ra_n,dec=dec_n,frame='icrs',unit='deg')

                separation = nuvcoord.separation(fuvcoord)
                separation_in_deg = separation.to(u.arcsec)

                if separation_in_deg < rad:

                    
                    print('MATCH FOUND')
                    print(f'Match found in coord: {nuvcoord}')
                    print(f'Index number of match: {count1}')


                    index_of_matches.append(count1)
                    # print(nuv_t.columns)

                    # print(f'Before {fuv_t["N2_M_A"].iloc[count0]}')
                    # fuv_t['N2_M_A'].iloc[count0] = nuv_t['N2_M_A'].iloc[count1]
                    # print(f'After {fuv_t["N2_M_A"].iloc[count0]}')

                    try:
                        fuv_t["N1_M_A"].iloc[count0] = nuv_t["N1_M_A"].iloc[count1]
                        fuv_t["N1_MER_A"].iloc[count0] = nuv_t['N1_MER_A'].iloc[count1]
                        fuv_t['N1_M_I'].iloc[count0] = nuv_t['N1_M_I'].iloc[count1]
                        fuv_t['N1_MER_I'].iloc[count0] = nuv_t['N1_MER_I'].iloc[count1]
                        fuv_t['N1_F_A'].iloc[count0] = nuv_t['N1_F_A'].iloc[count1]
                        fuv_t['N1_FER_A'].iloc[count0] = nuv_t['N1_FER_A'].iloc[count1]
                    except:
                        print('N1 not present')

                    try:
                        fuv_t["N2_M_A"].iloc[count0] = nuv_t["N2_M_A"].iloc[count1]
                        fuv_t["N2_MER_A"].iloc[count0] = nuv_t['N2_MER_A'].iloc[count1]
                        fuv_t['N2_M_I'].iloc[count0] = nuv_t['N2_M_I'].iloc[count1]
                        fuv_t['N2_MER_I'].iloc[count0] = nuv_t['N2_MER_I'].iloc[count1]
                        fuv_t['N2_F_A'].iloc[count0] = nuv_t['N2_F_A'].iloc[count1]
                        fuv_t['N2_FER_A'].iloc[count0] = nuv_t['N2_FER_A'].iloc[count1]
                    except:
                        print('N2 not present')

                    try:
                        fuv_t["N3_M_A"].iloc[count0] = nuv_t["N3_M_A"].iloc[count1]
                        fuv_t["N3_MER_A"].iloc[count0] = nuv_t['N3_MER_A'].iloc[count1]
                        fuv_t['N3_M_I'].iloc[count0] = nuv_t['N3_M_I'].iloc[count1]
                        fuv_t['N3_MER_I'].iloc[count0] = nuv_t['N3_MER_I'].iloc[count1]
                        fuv_t['N3_F_A'].iloc[count0] = nuv_t['N3_F_A'].iloc[count1]
                        fuv_t['N3_FER_A'].iloc[count0] = nuv_t['N3_FER_A'].iloc[count1]
                    except:
                        print('N3 not present')

                    try:
                        fuv_t["N5_M_A"].iloc[count0] = nuv_t["N5_M_A"].iloc[count1]
                        fuv_t["N5_MER_A"].iloc[count0] = nuv_t['N5_MER_A'].iloc[count1]
                        fuv_t['N5_M_I'].iloc[count0] = nuv_t['N5_M_I'].iloc[count1]
                        fuv_t['N5_MER_I'].iloc[count0] = nuv_t['N5_MER_I'].iloc[count1]
                        fuv_t['N5_F_A'].iloc[count0] = nuv_t['N5_F_A'].iloc[count1]
                        fuv_t['N5_FER_A'].iloc[count0] = nuv_t['N5_FER_A'].iloc[count1]
                    except:
                        print('N5 not present')

                    try:
                        fuv_t["N6_M_A"].iloc[count0] = nuv_t["N6_M_A"].iloc[count1]
                        fuv_t["N6_MER_A"].iloc[count0] = nuv_t['N6_MER_A'].iloc[count1]
                        fuv_t['N6_M_I'].iloc[count0] = nuv_t['N6_M_I'].iloc[count1]
                        fuv_t['N6_MER_I'].iloc[count0] = nuv_t['N6_MER_I'].iloc[count1]
                        fuv_t['N6_F_A'].iloc[count0] = nuv_t['N6_F_A'].iloc[count1]
                        fuv_t['N6_FER_A'].iloc[count0] = nuv_t['N6_FER_A'].iloc[count1]
                    except:
                        print('N6 not present')

                    print('\n')
            
            print(f'Index of matches: {index_of_matches}')

            # remove matched index:

            # print('Before dropping',nuv_t.shape)
            nuv_t.drop(index_of_matches,axis='rows',inplace=True)
            nuv_t.reset_index(drop=True,inplace=True)

            print('\n')
            # print(nuv_t.shape)

        concat_df = pd.concat([fuv_t,nuv_t])
        tabledf = Table.from_pandas(concat_df)

        binary = fits.BinTableHDU(data=tabledf)
        primary = fits.PrimaryHDU()

        final_hdu = fits.HDUList([primary,binary])
        final_hdu.writeto(f'/Users/swagat98/Documents/Combine_cat_2/Same_field/Merged_folder/{obj_list[m]}',overwrite=True)






# dataframe_build4()

def only_fuv_nuv():

    path_to_fuv = '/Users/swagat98/Documents/Combine_cat_2/test/'
    path_to_nuv = '/Users/swagat98/Documents/Combine_cat_2/NUV-cat/'

    fuv = glob(f'{ path_to_fuv }*.fits')
    nuv = glob(f'{ path_to_nuv }*.fits')

    f1 = []
    f2 = []
    f3 = []
    f5 = []
    f7 = []
    
    for x in fuv:

        filex = fits.open(x)

        header = filex[1].header
        data = filex[1].data

        filtername = header['FILTER']

        tabledf = Table(data)
        df = tabledf.to_pandas()

        df99remove = remove_99(df)

        if filtername == 'F1':
            f1.append(df99remove)
        
        if filtername == 'F2':
            f2.append(df99remove)

        if filtername == 'F3':
            f3.append(df99remove)

        if filtername == 'F5':
            f5.append(df99remove)

        if filtername == 'F7':
            f7.append(df99remove)

        
    final_list = []

    try:
        f1_concat = pd.concat(f1)
        f1_concat.rename(columns={'MAG_AUTO':'F1_M_A',
                                  'MAGERR_AUTO':'F1_MER_A',
                                  'MAG_ISO':'F1_M_I',
                                  'MAGERR_ISO':'F1_MER_I',
                                  'FLUX_AUTO':'F1_F_A',
                                  'FLUXERR_AUTO':'F1_FER_A'},inplace=True)
        
        final_list.append(f1_concat)
    except:
        print('No f1 files present')

    try:
        f2_concat = pd.concat(f2)
        f2_concat.rename(columns={'MAG_AUTO':'F2_M_A',
                                  'MAGERR_AUTO':'F2_MER_A',
                                  'MAG_ISO':'F2_M_I',
                                  'MAGERR_ISO':'F2_MER_I',
                                  'FLUX_AUTO':'F2_F_A',
                                  'FLUXERR_AUTO':'F2_FER_A'},inplace=True)
        final_list.append(f2_concat)
    except:
        print('No f2 files present')

    try:
        f3_concat = pd.concat(f3)
        f3_concat.rename(columns={'MAG_AUTO':'F3_M_A',
                                  'MAGERR_AUTO':'F3_MER_A',
                                  'MAG_ISO':'F3_M_I',
                                  'MAGERR_ISO':'F3_MER_I',
                                  'FLUX_AUTO':'F3_F_A',
                                  'FLUXERR_AUTO':'F3_FER_A'},inplace=True)
        final_list.append(f3_concat)
    except:
        print('No f3 files present')
    try:
        f5_concat = pd.concat(f5)
        f5_concat.rename(columns={'MAG_AUTO':'F5_M_A',
                                  'MAGERR_AUTO':'F5_MER_A',
                                  'MAG_ISO':'F5_M_I',
                                  'MAGERR_ISO':'F5_MER_I',
                                  'FLUX_AUTO':'F5_F_A',
                                  'FLUXERR_AUTO':'F5_FER_A'},inplace=True)
        final_list.append(f5_concat)
    except:
        print('No f5 files present')
    try:
        f7_concat = pd.concat(f7)
        f7_concat.rename(columns={'MAG_AUTO':'F7_M_A',
                                  'MAGERR_AUTO':'F7_MER_A',
                                  'MAG_ISO':'F7_M_I',
                                  'MAGERR_ISO':'F7_MER_I',
                                  'FLUX_AUTO':'F7_F_A',
                                  'FLUXERR_AUTO':'F7_FER_A'},inplace=True)
        final_list.append(f7_concat)
    except:
        print('No f7 files present')

    final_concat= pd.concat(final_list)
    print(final_concat.columns)
    print(final_concat.shape)

    
    finaldf = Table.from_pandas(final_concat)

    binary = fits.BinTableHDU(data=finaldf)
    binary.writeto('/Users/swagat98/Documents/Combine_cat_2/RESULT/fuv_merge.fits',overwrite=True)



    f1_n = []
    f2_n = []
    f3_n = []
    f5_n = []
    f6_n = []
    

    for y in nuv:
        # print(y)
        filey = fits.open(y)

        headery = filey[1].header
        datay = filey[1].data

        filtername_y = headery['FILTER']

        tabledf_y = Table(datay)
        dfy = tabledf_y.to_pandas()

        df99remove_y = remove_99(dfy)

        if filtername_y == 'F1':
            f1_n.append(df99remove_y)
        
        if filtername_y == 'F2':
            f2_n.append(df99remove_y)

        if filtername_y == 'F3':
            f3_n.append(df99remove_y)

        if filtername_y == 'F5':
            f5_n.append(df99remove_y)

        if filtername_y == 'F6':
            f6_n.append(df99remove_y)

    final_list_y = []
    
    try: 
        n1_concat = pd.concat(f1_n)

        n1_concat.rename(columns={'MAG_AUTO':'N1_M_A',
                                  'MAGERR_AUTO':'N1_MER_A',
                                  'MAG_ISO':'N1_M_I',
                                  'MAGERR_ISO':'N1_MER_I',
                                  'FLUX_AUTO':'N1_F_A',
                                  'FLUXERR_AUTO':'N1_FER_A'},inplace=True)
        
        final_list_y.append(n1_concat)

    except:
        print('No N1 files present')

    try: 
        n2_concat = pd.concat(f2_n)

        n2_concat.rename(columns={'MAG_AUTO':'N2_M_A',
                                  'MAGERR_AUTO':'N2_MER_A',
                                  'MAG_ISO':'N2_M_I',
                                  'MAGERR_ISO':'N2_MER_I',
                                  'FLUX_AUTO':'N2_F_A',
                                  'FLUXERR_AUTO':'N2_FER_A'},inplace=True)
        
        final_list_y.append(n2_concat)

    except:
        print('No N2 files present')


    try: 
        n3_concat = pd.concat(f3_n)

        n3_concat.rename(columns={'MAG_AUTO':'N3_M_A',
                                  'MAGERR_AUTO':'N3_MER_A',
                                  'MAG_ISO':'N3_M_I',
                                  'MAGERR_ISO':'N3_MER_I',
                                  'FLUX_AUTO':'N3_F_A',
                                  'FLUXERR_AUTO':'N3_FER_A'},inplace=True)
        
        final_list_y.append(n3_concat)

    except:
        print('No N3 files present')

    try: 
        n3_concat = pd.concat(f3_n)

        n3_concat.rename(columns={'MAG_AUTO':'N3_M_A',
                                  'MAGERR_AUTO':'N3_MER_A',
                                  'MAG_ISO':'N3_M_I',
                                  'MAGERR_ISO':'N3_MER_I',
                                  'FLUX_AUTO':'N3_F_A',
                                  'FLUXERR_AUTO':'N3_FER_A'},inplace=True)
        
        final_list_y.append(n3_concat)

    except:
        print('No N3 files present')

    try: 
        n5_concat = pd.concat(f5_n)

        n5_concat.rename(columns={'MAG_AUTO':'N5_M_A',
                                  'MAGERR_AUTO':'N5_MER_A',
                                  'MAG_ISO':'N5_M_I',
                                  'MAGERR_ISO':'N5_MER_I',
                                  'FLUX_AUTO':'N5_F_A',
                                  'FLUXERR_AUTO':'N5_FER_A'},inplace=True)
        
        final_list_y.append(n5_concat)

    except:
        print('No N5 files present')

    try: 
        n6_concat = pd.concat(f6_n)

        n6_concat.rename(columns={'MAG_AUTO':'N6_M_A',
                                  'MAGERR_AUTO':'N6_MER_A',
                                  'MAG_ISO':'N6_M_I',
                                  'MAGERR_ISO':'N6_MER_I',
                                  'FLUX_AUTO':'N6_F_A',
                                  'FLUXERR_AUTO':'N6_FER_A'},inplace=True)
        
        final_list_y.append(n6_concat)

    except:
        print('No N6 files present')

    
    final_df_y = pd.concat(final_list_y)

    print(final_df_y.columns)
    print(final_df_y)

    binarydf_y = Table.from_pandas(final_df_y)

    binary_n = fits.BinTableHDU(data=binarydf_y)
    primary_n = fits.PrimaryHDU()

    final_hdu_n = fits.HDUList([primary_n,binary_n])
    final_hdu_n.writeto('/Users/swagat98/Documents/Combine_cat_2/RESULT/nuv_merge.fits',overwrite=True)

# only_fuv_nuv()

def merged_and_nuv_output():
    
    path_to_merged = '/Users/swagat98/Documents/Combine_cat_2/Same_field/Merged_folder/'
    path_to_nuv_output = '/Users/swagat98/Documents/Combine_cat_2/Same_field/NUV_output/'

    merged = glob(f'{path_to_merged}*.fits')
    nuv_o = glob(f'{path_to_nuv_output}*.fits')

    merged_list=[]

    for x in merged:
        
        img0 = fits.open(x)

        header0 = img0[1].header
        data0 = img0[1].data

        tabledf0 = Table(data0)
        df0 = tabledf0.to_pandas()

        merged_list.append(df0)

    
    conc0 = pd.concat(merged_list)
    print(conc0)

    binaryconc0 = Table.from_pandas(conc0)

    binary0 = fits.BinTableHDU(data=binaryconc0)
    primary0 = fits.PrimaryHDU()

    final_hdu0 = fits.HDUList([primary0,binary0])

    final_hdu0.writeto('/Users/swagat98/Documents/Combine_cat_2/RESULT/merged_merge.fits',overwrite=True)


    nuv_merged_list = []
    for y in nuv_o:

        img1 = fits.open(y)

        data1 = img1[1].data

        tabledf1 = Table(data1)
        data1 = tabledf1.to_pandas()

        nuv_merged_list.append(data1)

    conc1 = pd.concat(nuv_merged_list)

    binaryconc1 = Table.from_pandas(conc1)

    binary1 = fits.BinTableHDU(data=binaryconc1)
    primary1 = fits.PrimaryHDU()

    final_hdu1 = fits.HDUList([primary1,binary1])
    final_hdu1.writeto('/Users/swagat98/Documents/Combine_cat_2/RESULT/nuv_merged_merge.fits',overwrite=True)



# merged_and_nuv_output()

def result_folder():

    path_to_result = '/Users/swagat98/Documents/Combine_cat_2/RESULT/'

    files = glob(f'{path_to_result}*.fits')

    final_list = []
    for x in files:

        img = fits.open(x)

        data = img[1].data

        tabledf = Table(data)
        df = tabledf.to_pandas()

        print(df.shape)
        final_list.append(df)

    
    final_conc = pd.concat(final_list)

    final_conc.drop('index',axis='columns',inplace=True)
    print(final_conc)
    print('Final shape:',final_conc.shape)
    print(final_conc.columns)

    final_conc.fillna(-999,inplace=True)

    finaltable = Table.from_pandas(final_conc)

    binary = fits.BinTableHDU(data=finaltable)
    primary = fits.PrimaryHDU()

    final = fits.HDUList([primary,binary])
    final.writeto('/Users/swagat98/Documents/Combine_cat_2/UVIT-cat.fits',overwrite=True)






result_folder()