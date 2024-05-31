import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table
import pandas as pd
from glob import glob
import os
import subprocess


global path_to_fuv
global path_to_nuv
global same_field_folder_path

path_to_nuv = '/Users/swagat98/Documents/Combine_cat/NUV-cat/' # path for the FUV files
path_to_fuv = '/Users/swagat98/Documents/Combine_cat/FUV-cat/' # path for the NUV files
same_field_folder_path = '/Users/swagat98/Documents/Combine_cat/Same_field/'  # path where we should create a new same field file



def check():

    fuv_files = glob(f'{path_to_fuv}*.fits')
    nuv_files = glob(f'{path_to_nuv}*.fits')


    final_fuv = []
    final_nuv = []

    f1 = []
    f2 = []
    f3 = []
    f5 = []
    f7 = []

    count0 = 0 

    while count0 < len(fuv_files):

        files = fits.open(fuv_files[count0])
        
        data = files[1].data
        header = files[1].header

        filters = header['FILTER']

        dfnew = Table(data)
        df = dfnew.to_pandas()

        final_fuv.append(df)

        if filters == 'F1':
            f1.append(df)
        
        if filters == 'F2':
            f2.append(df)

        if filters == 'F3':
            f3.append(df)

        if filters == 'F5':
            f5.append(df)

        if filters == 'F7':
            f7.append(df)


        count0 += 1

    f1_com = pd.concat(f1)
    f2_com = pd.concat(f2)
    f3_com = pd.concat(f3)
    f5_com = pd.concat(f5)
    # f7_com = pd.concat(f7)

    print(f1_com.shape)
    print(f2_com.shape)
    print(f3_com.shape)
    print(f5_com.shape)
    # print(f7.shape)


    count1 = 0


    n1 = []
    n2 = []
    n3 = []
    n5 = []
    n6 = []

    while count1 < len(nuv_files):
        
        files = fits.open(nuv_files[count1])

        data = files[1].data
        header = files[1].header

        filters = header['FILTER']

        dfnew = Table(data)
        df = dfnew.to_pandas()

        final_nuv.append(df)


        if filters == 'F1':
            n1.append(df)

        if filters == 'F2':
            n2.append(df)
        
        if filters == 'F3':
            n3.append(df)
        if filters == 'F5':
            n5.append(df)
        if filters == 'F6':
            n6.append(df)

        count1 += 1

    n1_com = pd.concat(n1)
    n2_com = pd.concat(n2)
    n3_com = pd.concat(n3)
    n5_com = pd.concat(n5)
    n6_com = pd.concat(n6)

    print(n1_com.shape)
    print(n2_com.shape)
    print(n3_com.shape)
    print(n5_com.shape)
    print(n6_com.shape)
    
    fuv_com = pd.concat(final_fuv)
    nuv_com = pd.concat(final_nuv)

    com = pd.concat([fuv_com,nuv_com])
    print(com.shape) #total number of detections raw 
    # this is just for checking


# check()





def separate_same_field_NUV_FUV():

    
    #we find the pointing information from both the folders and compare them both. 
    # We compare the NUV images first as they are of lower number than the FUV, and find the similar pointing images
    
    fuv_files = glob(f'{path_to_fuv}*.fits')
    nuv_files = glob(f'{path_to_nuv}*.fits')
    
    nuv_pointing = []
    nuv_obj = []
    
    count0 = 0
    while count0 < len(nuv_files):
        
        #open the nuv_files
        
        nuv_image = fits.open(nuv_files[count0])
#         nuv_image.info()
        
        # store the pointing information
        nuv_header = nuv_image[1].header
        
        ra_nuv = nuv_header['RA_PNT']
        dec_nuv = nuv_header['DEC_PNT']      
        nuv_pointing.append((ra_nuv,dec_nuv))
        
        obj_nuv = nuv_header['OBJECT'] #the stored information of the object name
        nuv_obj.append(obj_nuv)
        
        count0 += 1
    
#above are the nuv pointings and objects
    nuv_folder = list(zip(nuv_obj,nuv_pointing))
    print('\n','Images in NUVfolder: ',nuv_folder,'\n')



    fuv_pointing = []
    fuv_obj = []

    count1 = 0
    while count1 < len(fuv_files):
        
        fuv_image = fits.open(fuv_files[count1])

        fuv_header = fuv_image[1].header

        ra_fuv = fuv_header['RA_PNT']
        dec_fuv = fuv_header['DEC_PNT']
        fuv_pointing.append((ra_fuv,dec_fuv))
        # print(fuv_header['OBJECT'])

        fuv_obj.append(fuv_header['OBJECT'])

        count1 += 1
    
    fuv_folder = list(zip(fuv_obj,fuv_pointing))
    print('Images in FUV Folder: ',fuv_folder,'\n')



    same_field_files = []

    for i in range(len(nuv_pointing)):
        
        if nuv_pointing[i] in fuv_pointing:
            print(f'For NUV Pointing {nuv_folder[i]}: ')
            print(f'\tFUV pointing is present in FUV folder\n')
            same_field_files.append(nuv_obj[i])
            
        else:
            print(f'For NUV Pointing {nuv_folder[i]}: ')
            print('\tFUV pointing is not present in FUV folder\n')

    print('The similar field images are: ',same_field_files,'\n')


    


    if os.path.exists(f'{same_field_folder_path}'):
        print('Same field folder already exists. \n\tSkipping...')

    else:
        print('Same field folder not found! \n\tCreating one... \n')
        os.mkdir(f'{same_field_folder_path}')   #create folder named same_field
    
    ### Now copy the images to the same field folder, and delete from the original folder
        
        ## inside this has to be nuv and fuv folders again of same field, to use it later. 
        

    fuv_folder_inside_same_field_folder = f'{same_field_folder_path}/FUV_same_field/'  # creates new FUV folder inside the same field folder for the images with same pointings
    nuv_folder_inside_same_field_folder = f'{same_field_folder_path}/NUV_same_field/'  # creates new NUV folder inside the same field folder for the images with same pointings

    #creates FUV and NUV folders within Same Field Folder:
        #for FUV folder
    if os.path.exists(f'{fuv_folder_inside_same_field_folder}'):
        print('FUV folder already exists. \n\tSkipping...')
    else:
        print('FUV Folder not found! \n\tCreating one...\n')
        os.mkdir(f'{fuv_folder_inside_same_field_folder}')  # creates FUV folder inside the same field folder

        #for NUV folder
    if os.path.exists(f'{nuv_folder_inside_same_field_folder}'):
        print('NUV folder already exists. \n\tSkipping...\n')
    else:
        print('NUV Folder not found! \n\tCreating one...\n')
        os.mkdir(f'{nuv_folder_inside_same_field_folder}') #create NUV folder inside the same field folder



    # copy fuv files:
    for objnames in same_field_files:
        try:
            print(f'IMAGE NAME: {objnames}')  #image name
            
            fuv_copytosamefield = f'cp -r {path_to_fuv}{objnames}* {fuv_folder_inside_same_field_folder}'
            x = subprocess.run(fuv_copytosamefield,shell=True,check=True)
            print('\tFUV file copied successfully!')
            print('\tSubprocess Return Code',x.returncode,'\n')

            nuv_copytosamefield = f'cp -r {path_to_nuv}{objnames}* {nuv_folder_inside_same_field_folder}'
            y = subprocess.run(nuv_copytosamefield,shell=True,check=True)
            print('\tNUV file copied successfully!')
            print('\tSubprocess Return Code: ',y.returncode,'\n')

        except subprocess.CalledProcessError as e1:
            print('Error:', e1)
    

    #delete from the main folder
            
    for objnames in same_field_files:
        try:
            print(f'IMAGE NAME: {objnames}') #image name deleteing 
            
            fuv_delete_same_field = f'rm -r {path_to_fuv}{objnames}*'
            m = subprocess.run(fuv_delete_same_field,shell=True,check=True)
            print('\tFUV file deleted successfully from the main folder! ')
            print('\tSubprocess Return Code: ',m.returncode,'\n')

            nuv_delete_same_field = f'rm -r {path_to_nuv}{objnames}*'
            n = subprocess.run(nuv_delete_same_field,shell=True,check=True)
            print('\tNUV file deleted successfully from the main folder! ')
            print('\tSubprocess Return Code: ',n.returncode,'\n')

        except subprocess.CalledProcessError as e2:
            print('Error:', e2)

    # combine the fuv field images to a single catalog







    # sort the dataframes according to the filters
    # replace the existing present single catalog with the removed 99 ones
    # or take the complete datas in a list, and then use global variable to use the same in the next method
    # in the next method, we clean the catalog around a particular bright object, or any extended sources. 
    # (Before the previous step check images and whether the very bright sources has spreaded out detections)

    # the next step will be merging the datas and removing the duplicate objects after crossmatching (should have the same ra/dec of the field images)
            

   

    
# separate_same_field_NUV_FUV()




def remove_99(df):   #function to remove 99 values from the dataframe


    indexlist = []
    for i in range(len(df)):
        magnitude = df['MAG_AUTO'].iloc[i]
        if magnitude == 99:
            indexlist.append(i)
    
    df.drop(indexlist,axis='rows',inplace=True)
    return df
    

        
def sort_filters():

        new_fuvfolder = glob(f'{path_to_fuv}*')
        new_nuvfolder = glob(f'{path_to_nuv}*')

        F1_df = []
        F2_df = []
        F3_df = []
        F5_df = []
        F7_df = []

        N1_df = []
        N2_df = []
        N3_df = []
        N5_df = []
        N6_df = []
        
        
        combined_fuv = []

        count_0 = 0
        while count_0 < len(new_fuvfolder):
            

            img = fits.open(new_fuvfolder[count_0])
            # img.info()

            data = img[1].data
            header = img[1].header 

            

            image_filter = header['FILTER']
            objectname = header['OBJECT']
            
            imagedf = Table(data)
            df = imagedf.to_pandas()
            # print(df.keys())

            combined_fuv.append(df)

            if image_filter == 'F1':
                f1_99remove_df = remove_99(df)

                # print(f1_99remove_df.keys())

                #now change the column names and then append it to the F1_df
                f1_99remove_df.rename(columns={'MAG_AUTO':'F1_M_A',  # mag auto
                                               'MAGERR_AUTO':'F1_MER_A', #mag error auto 
                                               'MAG_ISO':'F1_M_I',  # mag iso
                                               'MAGERR_ISO':'F1_MER_I', #mag err iso
                                               'FLUX_AUTO':'F1_F_A', #
                                               'FLUXERR_AUTO':'F1_FER_A'},inplace=True)
                # print(f1_99remove_df.keys())
                F1_df.append(f1_99remove_df)


            if image_filter == 'F2':
                f2_99remove_df = remove_99(df)

                f2_99remove_df.rename(columns={'MAG_AUTO':'F2_M_A',
                                               'MAGERR_AUTO':'F2_MER_A',
                                               'MAG_ISO':'F2_M_I',
                                               'MAGERR_ISO':'F2_MER_I',
                                               'FLUX_AUTO':'F2_F_A',
                                               'FLUXERR_AUTO':'F2_FER_A'},
                                               inplace=True)
                
                # print(f2_99remove_df.keys())

                F2_df.append(f2_99remove_df)

            if image_filter == 'F3':
                f3_99remove_df = remove_99(df)

                f3_99remove_df.rename(columns={'MAG_AUTO':'F3_M_A',
                                               'MAGERR_AUTO':'F3_MER_A',
                                               'MAG_ISO':'F3_M_I',
                                               'MAGERR_ISO':'F3_MER_I',
                                               'FLUX_AUTO':'F3_F_A',
                                               'FLUXERR_AUTO':'F3_FER_A'},
                                               inplace=True)
                F3_df.append(f3_99remove_df)
            
            if image_filter == 'F5':
                f5_99remove_df = remove_99(df)

                f5_99remove_df.rename(columns={'MAG_AUTO':'F5_M_A',
                                               'MAGERR_AUTO':'F5_MER_A',
                                               'MAG_ISO':'F5_M_I',
                                               'MAGERR_ISO':'F5_MER_I',
                                               'FLUX_AUTO':'F5_F_A',
                                               'FLUXERR_AUTO':'F5_FER_A'},
                                               inplace=True)
                F5_df.append(f5_99remove_df)

                print(f5_99remove_df.keys())

            if image_filter == 'F7':
                f7_99remove_df = remove_99(df)

                f7_99remove_df.rename(columns={'MAG_AUTO':'F7_M_A',
                                               'MAGERR_AUTO':'F7_MER_A',
                                               'MAG_ISO':'F7_M_I',
                                               'MAGERR_ISO':'F7_MER_I',
                                               'FLUX_AUTO':'F7_F_A',
                                               'FLUXERR_AUTO':'F7_FER_A'},
                                               inplace=True)
                F7_df.append(f7_99remove_df)

                

            count_0 += 1

        x = pd.concat(combined_fuv)
        print(x)
# now for NUV 


        count_1 = 0

        while count_1 < len(new_nuvfolder):

            
            img = fits.open(new_nuvfolder[count_1])
            # img.info()

            data = img[1].data
            header = img[1].header 

            

            image_filter = header['FILTER']
            objectname = header['OBJECT']
            
            imagedf = Table(data)
            df = imagedf.to_pandas()
            # print(df.keys())

            if image_filter == 'F1':
                n1_99remove_df = remove_99(df)

                # print(f1_99remove_df.keys())

                #now change the column names and then append it to the F1_df
                n1_99remove_df.rename(columns={'MAG_AUTO':'N1_M_A',
                                               'MAGERR_AUTO':'N1_MER_A',
                                               'MAG_ISO':'N1_M_I',
                                               'MAGERR_ISO':'N1_MER_I',
                                               'FLUX_AUTO':'N1_F_A',
                                               'FLUXERR_AUTO':'N1_FER_A'},inplace=True)
                # print(f1_99remove_df.keys())
                N1_df.append(n1_99remove_df)
            
            if image_filter == 'F2':
                n2_99remove_df = remove_99(df)

                # print(f1_99remove_df.keys())

                #now change the column names and then append it to the F1_df
                n2_99remove_df.rename(columns={'MAG_AUTO':'N2_M_A',
                                               'MAGERR_AUTO':'N2_MER_A',
                                               'MAG_ISO':'N2_M_I',
                                               'MAGERR_ISO':'N2_MER_I',
                                               'FLUX_AUTO':'N2_F_A',
                                               'FLUXERR_AUTO':'N2_FER_A'},inplace=True)
                # print(f1_99remove_df.keys())
                N2_df.append(n2_99remove_df)


            if image_filter == 'F3':
                n3_99remove_df = remove_99(df)

                # print(f1_99remove_df.keys())

                #now change the column names and then append it to the F1_df
                n3_99remove_df.rename(columns={'MAG_AUTO':'N3_M_A',
                                               'MAGERR_AUTO':'N3_MER_A',
                                               'MAG_ISO':'N3_M_I',
                                               'MAGERR_ISO':'N3_MER_I',
                                               'FLUX_AUTO':'N3_F_A',
                                               'FLUXERR_AUTO':'N3_FER_A'},inplace=True)
                # print(f1_99remove_df.keys())
                N3_df.append(n3_99remove_df)


            if image_filter == 'F5':
                n5_99remove_df = remove_99(df)

                # print(f1_99remove_df.keys())

                #now change the column names and then append it to the F1_df
                n5_99remove_df.rename(columns={'MAG_AUTO':'N5_M_A',
                                               'MAGERR_AUTO':'N5_MER_A',
                                               'MAG_ISO':'N5_M_I',
                                               'MAGERR_ISO':'N5_MER_I',
                                               'FLUX_AUTO':'N5_F_A',
                                               'FLUXERR_AUTO':'N5_FER_A'},inplace=True)
                # print(f1_99remove_df.keys())
                N5_df.append(n5_99remove_df)


            if image_filter == 'F6':   #NUV has F6 instead of F7
                n6_99remove_df = remove_99(df)

                # print(f1_99remove_df.keys())

                #now change the column names and then append it to the F1_df1
                n6_99remove_df.rename(columns={'MAG_AUTO':'N6_M_A',
                                               'MAGERR_AUTO':'N6_MER_A',
                                               'MAG_ISO':'N6_M_I',
                                               'MAGERR_ISO':'N6_MER_I',
                                               'FLUX_AUTO':'N6_F_A',
                                               'FLUXERR_AUTO':'N6_FER_A'},inplace=True)
                # print(f1_99remove_df.keys())
                N6_df.append(n6_99remove_df)

            count_1 += 1


        print('Number of images in F1 FUV filter: ',len(F1_df))
        print('Number of images in F2 FUV filter: ',len(F2_df))
        print('Number of images in F3 FUV filter: ',len(F3_df))
        print('Number of images in F5 FUV filter: ',len(F5_df))
        print('Number of images in F7 FUV filter: ',len(F7_df))
        

        print('')
        print('Number of images in F1 NUV filter: ',len(N1_df))
        print('Number of images in F2 NUV filter: ',len(N2_df))
        print('Number of images in F3 NUV filter: ',len(N3_df))
        print('Number of images in F5 NUV filter: ',len(N5_df))
        print('Number of images in F6 NUV filter: ',len(N6_df))

        #concat the dataframes according to each filter
        try:
            f1_combined = pd.concat(F1_df)
            print(f'F1 dataframe: \n{f1_combined}')
            print(f'F1 dataframe shape: {f1_combined.shape}')
        except:         
            print('\nF1 filter\'s data not present')


        try:
            f2_combined = pd.concat(F2_df)
            print(f'F2 dataframe: \n{f2_combined}')
            print(f'F2 dataframe shape: {f2_combined.shape}')
        except:
            print('\nF2 filter\'s data not present')


        try:
            f3_combined = pd.concat(F3_df)
            print(f'F3 dataframe: \n{f3_combined}')
            print(f'F3 dataframe shape: {f3_combined.shape}')
        except:
            print('\nF3 filter\'s data not present')


        try:
            f5_combined = pd.concat(F5_df)
            print(f'F5 dataframe: \n{f5_combined}')
            print(f'F5 dataframe shape: {f5_combined.shape}')
        except:
            print('\nF5 filter\'s data not present' )


        try:
            f7_combined = pd.concat(F7_df)
            print(f'F7 dataframe: \n{f7_combined}')
            print(f'F7 dataframe shape: {f7_combined.shape}')
        except:
            print('\nF7 filter\'s data not present')


# combine dataframes for NUV filters
        
        try:
            n1_combined = pd.concat(N1_df)
            print(f'N1 dataframe: \n{n1_combined}')
            print(f'N1 dataframe shape: {n1_combined.shape}')
        except:
            print('\nN1 filter\'s data not present')


        try:
            n2_combined = pd.concat(N2_df)
            print(f'N2 dataframe: \n{n2_combined}')
            print(f'N2 dataframe shape: {n2_combined.shape}')
        except:
            print('\nN2 filter\'s data not present')


        try:
            n3_combined = pd.concat(N3_df)
            print(f'N3 dataframe: \n{n3_combined}')
            print(f'N3 dataframe shape: {n3_combined.shape}')
        except:
            print('\nN3 filter\'s data not present')

        try: 
            n5_combined = pd.concat(N5_df)
            print(f'N5 dataframe: \n{n5_combined}')
            print(f'N5 dataframe shape: {n5_combined.shape}')
        except:
            print('N5 filter\'s data not present')



        try:
            n6_combined = pd.concat(N6_df)
            print(f'N6 dataframe: \n{n6_combined}')
            print(f'N6 dataframe shape: {n6_combined.shape}')
        except:
            print('\nN6 filter\'s data not present')



# now rename the main number column:
            
        # Make a new dataframe where the columns are present according to the filters.
        #change the number to a different ID. like photoextractID from GUVcat.
                # we do the ID -  <Filter_name>-<>  
                    # this should be just before combining all to a single dataframe


# For F1:

        

        try:
            f1_combined.drop(columns=['NUMBER'],inplace=True)
            # print(f1_combined)
            
            df_length = len(f1_combined)
            id_number = np.arange(0,df_length)
            
            f1_id = np.vectorize(lambda x: 'F1-' + str(x))(id_number) #the ID to be created
            
            f1_combined['Source ID'] = f1_id

            

            first_column = f1_combined.pop('Source ID')
            f1_combined.insert(0,'Source ID',first_column)
            
            print('F1 dataframe after Source ID: \n')
            print(f1_combined)
        except:
            f1_combined = pd.DataFrame()
            print(f1_combined)
            print('F1 data not available')
        finally:
            print('ID created for the dataframe\n')

        
            
# For F2: 
       
        try: 
            f2_combined.drop(columns=['NUMBER'],inplace=True)
            # print(f2_combined)
            
            df_length = len(f2_combined)
            id_number = np.arange(0,df_length)
            
            f2_id = np.vectorize(lambda x: 'F2-' + str(x))(id_number) #the ID to be created
            
            f2_combined['Source ID'] = f2_id

            first_column = f2_combined.pop('Source ID')
            f2_combined.insert(0,'Source ID',first_column)


            print('F2 dataframe after Source ID: \n')
            print(f2_combined)
        except:
            f2_combined = pd.DataFrame()
            print(f2_combined)
            print('F2 data not available')
        finally:
            print('ID created for the dataframe\n')

# For F3: 
  
        try:
            f3_combined.drop(columns=['NUMBER'],inplace=True)

            df_length = len(f3_combined)
            id_number = np.arange(0,df_length)

            f3_id = np.vectorize(lambda x: 'F3-' + str(x))(id_number) # the ID to be created for F3

            f3_combined['Source ID'] = f3_id

            first_column = f3_combined.pop('Source ID')
            f3_combined.insert(0,'Source ID',first_column)

            print('F3 dataframe after Source ID: \n')
            print(f3_combined)
        except:
            f3_combined = pd.DataFrame()
            print(f3_combined)
            print('F3 data not available')
        finally:
            print('ID created for the dataframe\n')


# For F5: 
        
        try:
            f5_combined.drop(columns=['NUMBER'],inplace=True)

            df_length = len(f5_combined)
            id_number = np.arange(0,df_length)

            f5_id = np.vectorize(lambda x: 'F5-'+str(x))(id_number)

            f5_combined['Source ID'] = f5_id

            first_column = f5_combined.pop('Source ID')
            f5_combined.insert(0,'Source ID',first_column)

            print('F5 dataframe after Source ID: \n')
            print(f5_combined)        
        except:
            f5_combined = pd.DataFrame()
            print(f5_combined)
            print('F5 data not available')
        finally:
            print('ID created for the dataframe\n')

# for F7:
        
        try:
            f7_combined.drop(columns=['NUMBER'],inplace=True)

            df_length = len(f7_combined)
            id_number = np.arange(0,df_length)

            f7_id = np.vectorize(lambda x: 'F7-' + str(x))(id_number)

            f7_combined['Source ID'] = f7_id

            first_column = f7_combined.pop('Source ID')
            f7_combined.insert(0,'Source ID',first_column)
            
            print('F7 dataframe after Source ID: \n')
            print(f7_combined)
        except:
            f7_combined = pd.DataFrame()
            print(f7_combined)
            print('F7 data not available')
        finally:
            print('ID created for the dataframe\n')




#for N1:
        
        try:
            n1_combined.drop(columns=['NUMBER'],inplace=True)

            df_length = len(n1_combined)
            id_number = np.arange(0,df_length)

            n1_id = np.vectorize(lambda x: 'N1-' + str(x))(id_number)

            n1_combined['Source ID'] = n1_id

            first_column = n1_combined.pop('Source ID')
            n1_combined.insert(0,'Source ID',first_column)

            print('N1 dataframe after Source ID: \n')
            print(n1_combined)
        except:
            n1_combined = pd.DataFrame()
            print(n1_combined)
            print('N1 data not available')
        finally:
            print('ID created for the dataframe\n')


#For N2

        try:
            n2_combined.drop(columns=['NUMBER'],inplace=True)

            df_length = len(n2_combined)
            id_number = np.arange(0,df_length)

            n2_id = np.vectorize(lambda x: 'N2-' + str(x))(id_number)

            n2_combined['Source ID'] = n2_id

            first_column = n2_combined.pop('Source ID')
            n2_combined.insert(0,'Source ID',first_column)

            print('N2 dataframe after Source ID: \n')
            print(n2_combined)
        except:
            n2_combined = pd.DataFrame()
            print(n2_combined)
            print('N2 data not available')
        finally:
            print('ID created for the dataframe\n')
            

#for N3:
            
        try:
            n3_combined.drop(columns=['NUMBER'],inplace=True)

            df_length = len(n3_combined)
            id_number = np.arange(0,df_length)

            n3_id = np.vectorize(lambda x: 'N3-' + str(x))(id_number)

            n3_combined['Source ID'] = n3_id

            first_column = n3_combined.pop('Source ID')
            n3_combined.insert(0,'Source ID',first_column)

            print('N3 dataframe after Source ID: \n')
            print(n3_combined)
        except:
            n3_combined = pd.DataFrame()
            print(n3_combined)
            print('N3 data not available')
        finally:
            print('ID created for the dataframe\n')


# For N5:
        try:
            n5_combined.drop(columns=['NUMBER'],inplace=True)

            df_length = len(n5_combined)
            id_number = np.arange(0,df_length)

            n5_id = np.vectorize(lambda x: 'N5-' + str(x))(id_number)

            n5_combined['Source ID'] = n5_id

            first_column = f5_combined.pop('Source ID')
            f5_combined.insert(0,'Source ID',first_column)

            print('N5 dataframe after Source ID: \n')
            print(n5_combined)
        except:
            n5_combined = pd.DataFrame()
            print(n5_combined)
            print('N5 data not available')
        finally:
            print('ID created for the dataframe\n')

# For N6:
        try:
            n6_combined.drop(columns=['NUMBER'],inplace=True)

            df_length = len(n6_combined)
            id_number = np.arange(0,df_length)

            n6_id = np.vectorize(lambda x: 'N6-' + str(x))(id_number)

            n6_combined['Source ID'] = n6_id

            first_column = n6_combined.pop('Source ID')
            n6_combined.insert(0,'Source ID',first_column)

            print('N6 dataframe after Source ID: \n')
            print(n6_combined)
        except:
            n6_combined = pd.DataFrame()
            print(n6_combined)
            print('N6 data not available')
        finally:
            print('ID created for the dataframe\n')


        print('F1 filter source: ',len(f1_combined))
        print('F2 filter source: ',len(f2_combined))
        print('F3 filter source: ',len(f3_combined))
        print('F5 filter source: ',len(f5_combined))
        print('F7 filter source: ',len(f7_combined))

        print('NUV F1 filter source: ',len(n1_combined))
        print('NUV F2 filter source: ',len(n2_combined))
        print('NUV F3 filter source: ',len(n3_combined))
        print('NUV F5 filter source: ',len(n5_combined))
        print('NUV F6 filter source: ',len(n6_combined))


        fuv = len(f1_combined)+len(f2_combined)+len(f3_combined)+len(f5_combined)+len(f7_combined)
        nuv = len(n1_combined)+len(n2_combined)+len(n3_combined)+len(n5_combined)+len(n6_combined)

        print(f'Total Number of FUV sources: {fuv}')

        print(f'Total Number of NUV sources: {nuv}')

        print(f'Total number of sources: {fuv+nuv}') 
        # print(sum([len(f1_combined),len(f2_combined),len(f3_combined),len(f5_combined),len(f7_combined),len(n1_combined),len(n2_combined),len(n3_combined),len(n5_combined),len(n6_combined)]))
        dataframes_list = []

        

        if len(f1_combined) != 0:
            dataframes_list.append(f1_combined)
        else:
            print('F1 dataframe not available')

        if len(f2_combined) != 0:
            dataframes_list.append(f2_combined)
        else:
            print('F2 dataframe not available')

        if len(f3_combined) != 0:
            dataframes_list.append(f3_combined)
        else:
            print('F3 dataframe not available')

        if len(f5_combined) != 0:
            dataframes_list.append(f5_combined)
        else:
            print('F5 dataframe not available')

        if len(f7_combined) != 0:
            dataframes_list.append(f7_combined)
        else:
            print('F7 dataframe not available')
#nuv
        if len(n1_combined) != 0:
            dataframes_list.append(n1_combined)
        else:
            print('N1 dataframe not available')

        if len(n2_combined) != 0:
            dataframes_list.append(n2_combined)
        else:
            print('N2 dataframe not available')

        if len(n3_combined) != 0:
            dataframes_list.append(n3_combined)
        else:
            print('N3 dataframe not available')

        if len(n5_combined) != 0:
            dataframes_list.append(n5_combined)
        else:
            print('N5 dataframe not available')

        if len(n6_combined) != 0:
            dataframes_list.append(n6_combined)
        else:
            print('N6 dataframe not available')

        final_dataframe = pd.concat(dataframes_list, ignore_index=True)


        renamed_columns = {'Source ID':'SRC_ID',
                           'ALPHA_J2000':'RA',
                           'DELTA_J2000':'DEC',
                           'ERRA_IMAGE':'ERRA_IMG',
                           'ERRB_IMAGE':'ERRB_IMG',
                           'ERRA_WORLD':'ERRA_WLD',
                           'ERRB_WORLD':'ERRB_WLD',
                           'FWHM_WORLD':'FWHM_WLD',
                           'FWHM_IMAGE':'FWHM_IMG',
                           'CLASS_STAR':'CLS_STR',
                           'ELLIPTICITY':'ELLIP',
                           'THETA_J2000':'THEJ2000',
                           'KRON_RADIUS':'KRON_RAD',
                           'Distance from FOV Centre':'DIST_FOV'}
        final_dataframe.rename(columns=renamed_columns,inplace=True)



        print(final_dataframe.keys())
        print(final_dataframe.shape)


        final_binary_table = Table.from_pandas(final_dataframe)
        final_binary_table.write('/Users/swagat98/UVIT-cat.fits',format='fits',overwrite=True)


def add_license_removenan(path):

    img = fits.open(path)

    data = img[1].data
    header = img[1].header
    pri_header = img[0].header
    pri_data = img[0].data

    # print(header.keys)

    newdf = Table(data)
    df = newdf.to_pandas()

    # print(df.info())

    nanremovedf = df.fillna(-999)
    
    nanremovetable = Table.from_pandas(nanremovedf)
    # print(type(nanremovedf))
    img.close()
    # fits.writeto('/Users/swagat98/newtest.fits',data,header,overwrite=True)
    primary = fits.PrimaryHDU(data=pri_data,header=pri_header)
    binary = fits.BinTableHDU(data=nanremovetable,header=header)

    finalfits = fits.HDUList([primary,binary])
    finalfits.writeto('/Users/swagat98/UVIT-cat.fits',overwrite=True)




def add_keyword_comment():
    


    img = fits.open( '/Users/swagat98/UVIT-cat.fits',mode='update')
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
    final.writeto('/Users/swagat98/UVIT-cat.fits',overwrite=True)

        
def add_license(path):

    path_to_license = '/Users/swagat98/Documents/Copyright 2024 Swagat Bordoloi'
    
    license_list = []
    with open(path_to_license,'r') as file:

        for line in file:
            print(type( line.strip()) )
            license_list.append(line.strip())

    img= fits.open(path,mode='update')
    img.info()

    primary_header = img[0].header
    primary_data = img[0].header
    binary_header = img[1].header
    binary_data = img[1].data

    # print(primary_header.keys)

    # primary_header.add_blank()
    
    primary_header[''] = '31-05-24'
    primary_header[''] = 'Version 1.0'

    for x in license_list:
        primary_header[''] = x 

    img.flush()
    print(primary_header.keys)
    # print(primary_header.keys)
    # primary = fits.PrimaryHDU(primary_data,primary_header)
    # primary = fits.PrimaryHDU(primary_data,primary_header)
    # binary = fits.BinTableHDU(binary_data,header=binary_header)

    # final = fits.HDUList([primary,binary])
    # final.writeto('/Users/swagat98/UVIT-cat.fits',overwrite=True)

sort_filters()
add_license_removenan('/Users/swagat98/UVIT-cat.fits')
add_keyword_comment()
add_license('/Users/swagat98/UVIT-cat.fits')



# Also make for NUV field images.
# maybe do the NUV files again?

