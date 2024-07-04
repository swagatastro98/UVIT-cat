
from astropy.io import fits
import numpy as np
import pandas as pd
from astropy import units as u
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy.coordinates import Angle


# Run script to make a merged catalog for al the images found from Final_overlapping_images.py (filterwise)
# Run separately for each folder, made filterwise 
# Run Merge_all_files.py




def merged_cat(path):

    img = fits.open(path)
    img.info()


    data = img[1].data
    header = img[1].header
    pri_header = img[0].header
    pri_data = img[0].data

    table = Table(data)
    df = table.to_pandas()

    print(df)
    try:
        f1 = df[df['F1_M_A']>0] 
    except:
        print('No f1 data')
    try:
        f2 = df[df['F2_M_A']>0] 

    except:
        print('No f2 data')
    try:
        f3 = df[df['F3_M_A']>0] 
    except:
        print('No f3 data')

    try:
        f5 = df[df['F5_M_A']>0] 
    except:
        print('No f5 data')
    try:
        f7 = df[df['F7_M_A']>0] 
    except:
        print('No f7 data')

    #length of the dataframes
    # print(len(f1))
    print(len(f2))
    # print(len(f3))
    # print(len(f5))
    # print(len(f7))

    #for f7


    ra_dec_f1 = list(SkyCoord(ra=f1['RA'],dec=f1['DEC'],unit='deg',frame='icrs'))
    ra_dec_f2 = list(SkyCoord(ra=f2['RA'],dec=f2['DEC'],unit='deg',frame='icrs'))
    ra_dec_f3 = list(SkyCoord(ra=f3['RA'],dec=f3['DEC'],unit='deg',frame='icrs'))
    ra_dec_f5 = list(SkyCoord(ra=f5['RA'],dec=f5['DEC'],unit='deg',frame='icrs'))
    ra_dec_f7 = list(SkyCoord(ra=f7['RA'],dec=f7['DEC'],unit='deg',frame='icrs'))


    radius_of_merging = 3
    angle_in_degree = Angle(radius_of_merging,unit=u.arcsec)


    print('\nDoing for filter F1\n')

    coordinates_of_matches_f1 = []

    for x1 in range(len(f1)):

        index_of_match_coords_f1 = [] # containts the matched coordinates

        print(f'Doing for coord number: {x1}')
        print(f'Coordinate: {ra_dec_f1[x1]}')

        if ra_dec_f1[x1] in coordinates_of_matches_f1:
            pass
        else:

            for x1inside in range(len(f1)):
                ra_inside_f1 = f1['RA'].iloc[x1inside]
                dec_inside_f1 = f1['DEC'].iloc[x1inside]
                coord_inside_f1 = SkyCoord(ra=ra_inside_f1,dec=dec_inside_f1,frame='icrs',unit='deg')
                separation_f1 = ra_dec_f1[x1].separation(coord_inside_f1)
                separation_in_deg_f1 = separation_f1.to(u.arcsec)

                if separation_in_deg_f1 < angle_in_degree:

                    if separation_in_deg_f1 == 0:
                        pass
                    else:
                        print(f'The separation is {separation_in_deg_f1}')
                        print(f'The matched coordinate is: {coord_inside_f1}')
                        coordinates_of_matches_f1.append(coord_inside_f1)

                        if ra_dec_f1[x1] in coordinates_of_matches_f1:
                            pass
                        else:
                            index_of_match_coords_f1.append(x1inside)

            print(f'Matched coordinate list: {coordinates_of_matches_f1}')

            print(index_of_match_coords_f1)

            f1.drop(index_of_match_coords_f1,axis='rows',inplace=True)
            f1.reset_index(drop=True,inplace=True)
            print(len(f1))

            print('\n\n')

    finaldf_f1 = Table.from_pandas(f1)

    binary_f1 = fits.BinTableHDU(data=finaldf_f1,header=header)
    primary_f1 = fits.PrimaryHDU(data=pri_data,header=pri_header)

    final_f1 = fits.HDUList([primary_f1,binary_f1])
    final_f1.writeto('/Users/swagat98/Documents/UVIT-cat/f1mergedtable.fits',overwrite=True)


        # For filter F2

    print('\nDoing for Filter F2\n')

    coordinates_of_matches_f2 = []

    for x2 in range(len(f2)):

        index_of_match_coords_f2 = []

        print(f'Doing for coord number: {x2}')
        print(f'Coordinate: {ra_dec_f2[x2]}')

        if ra_dec_f2[x2] in coordinates_of_matches_f2:
            pass
        else:
            for x2inside in range(len(f2)):
                ra_inside_f2 = f2['RA'].iloc[x2inside]
                dec_inside_f2 = f2['DEC'].iloc[x2inside]
                coord_inside_f2 = SkyCoord(ra=ra_inside_f2,dec=dec_inside_f2,frame='icrs',unit='deg')
                separation_f2 = ra_dec_f2[x2].separation(coord_inside_f2)
                separation_in_deg_f2 = separation_f2.to(u.arcsec)


                if separation_in_deg_f2 < angle_in_degree:

                    if separation_in_deg_f2 == 0:
                        pass
                    else:
                        print(f'The separation is {separation_in_deg_f2}')
                        print(f'The matched coordinate is: {coord_inside_f2}')
                        coordinates_of_matches_f2.append(coord_inside_f2)

                        if ra_dec_f2[x2] in coordinates_of_matches_f2:
                            pass
                        else:
                            index_of_match_coords_f2.append(x2inside)


            print(f'Matched coordinate list: {coordinates_of_matches_f2}')

            print(index_of_match_coords_f2)

            f2.drop(index_of_match_coords_f2,axis='rows',inplace=True)
            f2.reset_index(drop=True,inplace=True)
            print(len(f2))

            print('\n\n')

    finaldf_f2 = Table.from_pandas(f2)

    binary_f2 = fits.BinTableHDU(data=finaldf_f2,header=header)
    primary_f2 = fits.PrimaryHDU(data=pri_data,header=pri_header)

    final_f2 = fits.HDUList([primary_f2,binary_f2])
    final_f2.writeto('/Users/swagat98/Documents/UVIT-cat/f2mergedtable.fits',overwrite=True)


    print('\nDoing for filter F3\n')

    coordinates_of_matches_f3 = []

    for x3 in range(len(f3)):

        index_of_match_coords_f3 = []

        print(f'Doing for coord number: {x3}')
        print(f'Coordinate: {ra_dec_f3[x3]}')

        if ra_dec_f3[x3] in coordinates_of_matches_f3:
            pass
        else:

            for x3inside in range(len(f3)):
                ra_inside_f3 = f3['RA'].iloc[x3inside]
                dec_inside_f3 = f3['DEC'].iloc[x3inside]
                coord_inside_f3 = SkyCoord(ra=ra_inside_f3,dec=dec_inside_f3,frame='icrs',unit='deg')
                separation_f3 = ra_dec_f3[x3].separation(coord_inside_f3)
                separation_in_deg_f3 = separation_f3.to(u.arcsec)

                if separation_in_deg_f3 < angle_in_degree:

                    if separation_in_deg_f3 == 0:
                        pass
                    else:
                        print(f'The separation is {separation_in_deg_f3}')
                        print(f'The matched coordinate is: {coord_inside_f3}')
                        coordinates_of_matches_f3.append(coord_inside_f3)

                        if ra_dec_f3[x3] in coordinates_of_matches_f3:
                            pass

                        else:
                            index_of_match_coords_f3.append(x3inside)


            print(f'Matched coordinate list: {coordinates_of_matches_f3}')
            print(index_of_match_coords_f3)

            f3.drop(index_of_match_coords_f3,axis='rows',inplace=True)
            f3.reset_index(drop=True,inplace=True)
            print(len(f3))

            print('\n\n')

    finaldf_f3 = Table.from_pandas(f3)

    binary_f3 = fits.BinTableHDU(data=finaldf_f3,header=header)
    primary_f3 = fits.PrimaryHDU(data=pri_data,header=pri_header)

    final_f3 = fits.HDUList([primary_f3,binary_f3])
    final_f3.writeto('/Users/swagat98/Documents/UVIT-cat/f333rgedtable.fits',overwrite=True)




    print('\nDoing for filter F5\n')

    coordinates_of_matches_f5 = []

    for x5 in range(len(f5)):
        
        index_of_match_coords_f5= []

        print(f'Doing for coord number: {x5}')
        print(f'Coordinate: {ra_dec_f5[x5]}')

        if ra_dec_f5[x5] in coordinates_of_matches_f5:
            pass
        else:

            for x5inside in range(len(f5)):
                ra_inside_f5 = f5['RA'].iloc[x5inside]
                dec_inside_f5 = f5['DEC'].iloc[x5inside]
                coord_inside_f5 = SkyCoord(ra=ra_inside_f5,dec=dec_inside_f5,frame='icrs',unit='deg')
                separation_f5 = ra_dec_f5[x5].separation(coord_inside_f5)
                separation_in_deg_f5 = separation_f5.to(u.arcsec)

                if separation_in_deg_f5< angle_in_degree:

                    if separation_in_deg_f5 == 0:
                        pass

                    else:
                        print(f'The separation is {separation_in_deg_f5}')
                        print(f'The matched coordinate is: {coord_inside_f5}')
                        coordinates_of_matches_f5.append(coord_inside_f5)

                        if ra_dec_f5[x5] in coordinates_of_matches_f5:
                            pass
                        else:
                            index_of_match_coords_f5.append(x5inside)


            print(f'Matched coordinate list: {coordinates_of_matches_f5}')
            print(index_of_match_coords_f5)

            f5.drop(index_of_match_coords_f5,axis='rows',inplace=True)
            f5.reset_index(drop=True,inplace=True)
            print(len(f5))
            print('\n\n')

    finaldf_f5 = Table.from_pandas(f5)

    binary_f5 = fits.BinTableHDU(data=finaldf_f5,header=header)
    primary_f5 = fits.PrimaryHDU(data=pri_data,header=pri_header)

    final_f5 = fits.HDUList([primary_f5,binary_f5])
    final_f5.writeto('/Users/swagat98/Documents/UVIT-cat/f5mergedtable.fits',overwrite=True)



    print('\nDoing for Filter F7\n')

    coordinates_of_matches_f7 = []

    for x7 in range(len(f7)):
        
        index_of_match_coords_f7 = [] #this contains the matched coordinates

        print(f'Doing for coord number: {x7}')
        print(f'Coordinate: {ra_dec_f7[x7]}')

        if ra_dec_f7[x7] in coordinates_of_matches_f7:
            pass
        else:

            for x7inside in range(len(f7)):
                ra_inside_f7 = f7['RA'].iloc[x7inside]
                dec_inside_f7 = f7['DEC'].iloc[x7inside]
                coord_inside_f7 = SkyCoord(ra=ra_inside_f7,dec=dec_inside_f7,frame='icrs',unit='deg')
                separation_f7 = ra_dec_f7[x7].separation(coord_inside_f7)
                separation_in_deg_f7 = separation_f7.to(u.arcsec)
                # print(separation_in_deg_f7)

                # if match is found store it in a list the index values 

                if separation_in_deg_f7 < angle_in_degree:
                    
                    if separation_in_deg_f7 == 0:
                        pass
                    else:
                        print(f'The separation is {separation_in_deg_f7}')
                        print(f'The matched coordinate is: {coord_inside_f7}')
                        coordinates_of_matches_f7.append(coord_inside_f7)

                        if ra_dec_f7[x7] in coordinates_of_matches_f7: #if the searching coordinate is already wihtin the list of the matched coordinate before ignore it
                            pass
                        else:
                            index_of_match_coords_f7.append(x7inside)
            
            print(f'Matched coordinate list:{coordinates_of_matches_f7}')

    # ignore all the coordinates that are found in the first match


                # ask what to do next
            print(index_of_match_coords_f7)

            f7.drop(index_of_match_coords_f7,axis='rows',inplace=True)
            f7.reset_index(drop=True,inplace=True)
            print(len(f7))


            print('\n\n')
        # else:
        #     print('No match found')

    finaldf_f7 = Table.from_pandas(f7)

    binary_f7 = fits.BinTableHDU(data=finaldf_f7,header=header)
    primary_f7 = fits.PrimaryHDU(data=pri_data,header=pri_header)

    final_f7 = fits.HDUList([primary_f7,binary_f7])
    final_f7.writeto('/Users/swagat98/Documents/UVIT-cat/f7mergedtable.fits',overwrite=True)


merged_cat('/Users/swagat98/MAIN_UVIT_cat/UVIT-cat_f2.fits')
# merged_cat('/Users/swagat98/Documents/testdf_f1.fits')
# merged_cat('/Users/swagat98/Documents/testdf_f2.fits')
# merged_cat('/Users/swagat98/Documents/testdf_f3.fits')
# merged_cat('/Users/swagat98/Documents/testdf_f5.fits')
# merged_cat('/Users/swagat98/Documents/testdf_f7.fits')

