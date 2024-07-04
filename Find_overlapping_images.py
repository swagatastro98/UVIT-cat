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

path_to_nuv = '/Users/swagat98/Documents/Combine_cat/NUV-cat/' # path for the FUV files
path_to_fuv = '/Users/swagat98/Documents/Combine_cat/FUV-cat/' # path for the NUV files
same_field_folder_path = '/Users/swagat98/Documents/Combine_cat/Same_field/'  # path where we should create a new same field file


def separation_in_filters():
    
    fuv = glob(f'{ path_to_fuv }*.fits')

    print(fuv)

    f1_folder = '/Users/swagat98/Documents/Combine_cat/F1-FUV/'
    f2_folder = '/Users/swagat98/Documents/Combine_cat/F2-FUV/'
    f3_folder = '/Users/swagat98/Documents/Combine_cat/F3-FUV/'
    f5_folder = '/Users/swagat98/Documents/Combine_cat/F5-FUV/'
    f7_folder = '/Users/swagat98/Documents/Combine_cat/F7-FUV/'

    if os.path.exists(f'{f1_folder}'):
        print('Folder exists')

    else:
        print('Creating one')
        os.mkdir(f'{f1_folder}')


    if os.path.exists(f'{f2_folder}'):
        print('Folder exists')

    else:
        print('Creating one')
        os.mkdir(f'{f2_folder}')
               

    if os.path.exists(f'{f3_folder}'):
        print('Folder exists')

    else:
        print('Creating one')
        os.mkdir(f'{f3_folder}')

    if os.path.exists(f'{f5_folder}'):
        print('Folder exists')

    else:
        print('Creating one')
        os.mkdir(f'{f5_folder}')


    if os.path.exists(f'{f7_folder}'):
        print('Folder exists')

    else:
        print('Creating one')
        os.mkdir(f'{f7_folder}')



    for x in fuv:

        img = fits.open(x)

        header = img[1].header

        filter = header['FILTER']


        if filter == 'F1':

            print('Copying to F1')

            copy_f1 = f'cp -r {x} {f1_folder}'
            subprocess.run(copy_f1,shell=True,check=True)

        if filter == 'F2':

            print('Copying to F2')

            copy_f2 = f'cp -r {x} {f2_folder}'
            subprocess.run(copy_f2,shell=True,check=True)
        if filter == 'F3':

            print('Copying to F3')

            copy_f3 = f'cp -r {x} {f3_folder}'
            subprocess.run(copy_f3,shell=True,check=True)
        if filter == 'F5':

            print('Copying to F5')

            copy_f5 = f'cp -r {x} {f5_folder}'
            subprocess.run(copy_f5,shell=True,check=True)
        if filter == 'F7':

            print('Copying to F7')

            copy_f7 = f'cp -r {x} {f7_folder}'
            subprocess.run(copy_f7,shell=True,check=True)





separation_in_filters()



def check_merged():

    path_f1 = '/Users/swagat98/Documents/Combine_cat/F1-FUV/'
    path_f2 = '/Users/swagat98/Documents/Combine_cat/F2-FUV/'
    path_f3 = '/Users/swagat98/Documents/Combine_cat/F3-FUV/'



    fuv = glob(f'{path_f2}*.fits')
    # print(len(fuv))
    nuv = glob(f'{path_to_nuv}*.fits')


    close_point_folder_f1 = '/Users/swagat98/Documents/Combine_cat/close_folder_f1/'
    close_point_folder_f2 = '/Users/swagat98/Documents/Combine_cat/close_folder_f2/'
    close_point_folder_f3 = '/Users/swagat98/Documents/Combine_cat/close_folder_f3/'

    if os.path.exists(f'{close_point_folder_f1}'):
        print('Folder exists. Skipping')

    else:
        print('Folder not found creating one')
        os.mkdir(f'{close_point_folder_f1}')


    if os.path.exists(f'{close_point_folder_f2}'):
        print('Folder exists. Skipping')

    else:
        print('Folder not found creating one')
        os.mkdir(f'{close_point_folder_f2}')

    if os.path.exists(f'{close_point_folder_f3}'):
        print('Folder exists. Skipping')

    else:
        print('Folder not found creating one')
        os.mkdir(f'{close_point_folder_f3}')

    field_dist = 14
    field_dist_in_deg = Angle(field_dist, unit=u.arcmin)
    print(f'Searching for fields within {field_dist_in_deg}')





    fuv_pointings_list = []
    out_object_list = []

    for x in fuv:

        fuv_img = fits.open(x)

        header = fuv_img[1].header
        data  = fuv_img[1].data

        newdf = Table(data)
        df  = newdf.to_pandas()        

        ra =  header['RA_PNT']
        dec = header['DEC_PNT']
        obj = header['OBJECT']

        
        pointing = SkyCoord(ra=ra,dec=dec,unit='deg',frame='icrs')
        # print(pointing)

        fuv_pointings_list.append(pointing)
        out_object_list.append(obj)

    # print(len(fuv_pointings_list))

    same_pointing_images = []
    coordinates_of_matches = []

    for out_point in range( len( fuv_pointings_list ) ):

        # pointing that checks

        print('Doing for pointing : ',out_object_list[out_point])
        
        for in_point in range(len(fuv)):

            # print(f'{ in_point }')
            in_img = fits.open(fuv[in_point])

            in_header = in_img[1].header
            in_data = in_img[1].data

            in_ra = in_header['RA_PNT']
            in_dec = in_header['DEC_PNT']
            in_object = in_header['OBJECT']
            #
            print('\tChecking object: ',in_object)


            in_pointing = SkyCoord(ra=in_ra,dec=in_dec,unit='deg',frame='icrs')

            separation_of_fields = fuv_pointings_list[out_point].separation(in_pointing)
            separation_in_deg = separation_of_fields.to(u.arcmin)
            #
            print('\t\t Separation in degrees is: ',separation_in_deg)

            if separation_in_deg < field_dist_in_deg:

                if separation_in_deg == 0:
                    
                    pass

                else:

                    print(f'The matched coordinate is: {in_pointing}')
                    print(f'The separation in arcmin is:{separation_in_deg} ')
                    coordinates_of_matches.append(in_pointing)

                    # if fuv_pointings_list[out_point] in coordinates_of_matches:
                        # pass
                    # else:
                    same_pointing_images.append(in_object)




    print(f'Same pointing images are: {same_pointing_images}')
    #copy the files:

    for obj_names in same_pointing_images:

        try:
            print('Copying close images')
    #         
            print(obj_names)
            fuv_copy = f'cp -r {path_f2}{obj_names}* {close_point_folder_f2}'
            x = subprocess.run(fuv_copy,shell=True,check=True)
            print('File copied')
        #
        except subprocess.CalledProcessError as e1:
            print('Error:',e1)
    #
    #

    for obj_names in same_pointing_images:
        try:
    #         
            fuv_del = f'rm -r {path_f2}{obj_names}*'
            y = subprocess.run(fuv_del,shell=True,check=True)
            print('File deleted')
        except:
            print('Error')



    # copy close pointing 

    # for the same pointings found, copy those files to a new folder
    # for the rest combine to make the UVIT-catalog,

    # if separation is not 0;
        #move to a different folder
    # if separation is 0;
        #

    # in that fuv folder, combine to make a single catalog





# check_merged()
