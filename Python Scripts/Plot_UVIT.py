import astropy 
import numpy as np
import pandas as pd
from astropy.io import fits
from glob import glob
import matplotlib.pyplot as plt
import statistics as stats 
from astropy.coordinates import SkyCoord, SkyCoordInfo
from astropy import units as u
from astropy.table import Table
from astropy import units as u
import seaborn as sns
import os





def count_vs_magnitude_distribution():
    
    global path
    global path2
    path = '/Users/swagat98/Documents/Combine_cat/FUV-cat/*'
    path2= '/Users/swagat98/Documents/Combine_cat/NUV-cat/*'

    path_fits = '/Users/swagat98/UVIT-cat.fits'

    images_fuv = glob(path)
    # print(images)
    
    images_nuv = glob(path2)
    # print(images0)
    
    
    count_f1 = 0
    count_f2 = 0
    count_f3 = 0
    count_f4 = 0
    count_f5 = 0
    count_f6 = 0
    count_f7 = 0 
    
    nuv_f1 = 0
    nuv_f2 = 0
    nuv_f3 = 0
    nuv_f5 = 0
    nuv_f6 = 0
    
    # combined_magnitude = []
    combined_ab_magnitude_F1 = []
    combined_ab_magnitude_F2 = []
    combined_ab_magnitude_F3 = []
    combined_ab_magnitude_F4 = []
    combined_ab_magnitude_F5 = []
    combined_ab_magnitude_F6 = []
    combined_ab_magnitude_F7 = []
    
    combined_nuv_ab_mag_f1 = []
    combined_nuv_ab_mag_f2 = []
    combined_nuv_ab_mag_f3 = []
    combined_nuv_ab_mag_f5 = []
    combined_nuv_ab_mag_f6 = []
    
    new_df = pd.DataFrame()
    nuv_df = pd.DataFrame()


    for j in range(len(images_fuv)):
        
        # print(f'{images[j]}\n')
        
        variable = fits.open(images_fuv[j])
        # print(variable)
        
        image = variable[1].data
        header = variable[1].header
        header0 = variable[0].header

        filters = header['FILTER']
        
        mag = image['MAG_AUTO']

        if filters == 'F1':

            dataframef1 = pd.DataFrame(mag)
            combined_ab_magnitude_F1.append(dataframef1)
        
        if filters == 'F2':

            dataframef2 = pd.DataFrame(mag)
            combined_ab_magnitude_F2.append(dataframef2)
        
        if filters == 'F3':

            dataframef3 = pd.DataFrame(mag)
            combined_ab_magnitude_F3.append(dataframef3)
        
        if filters == 'F5':

            dataframef5 = pd.DataFrame(mag)
            combined_ab_magnitude_F5.append(dataframef5)

        if filters == 'F7':

            dataframef7 = pd.DataFrame(mag)
            combined_ab_magnitude_F7.append(dataframef7)
            # print(dataframef1)

    

    for k in range(len(images_nuv)):

        variable_nuv = fits.open(images_nuv[k])
        # variable.info()

        image_nuv = variable_nuv[1].data
        header_nuv = variable_nuv[1].header

        filters = header_nuv['FILTER']

        mag_nuv = image_nuv['MAG_AUTO']

        if filters == 'F1':

            dataframe_nuvf1 = pd.DataFrame(mag_nuv)
            combined_nuv_ab_mag_f1.append(dataframe_nuvf1)
        
        if filters == 'F2':

            dataframe_nuvf2 = pd.DataFrame(mag_nuv)
            combined_nuv_ab_mag_f2.append(dataframe_nuvf2)
        if filters == 'F3':

            dataframe_nuvf3 = pd.DataFrame(mag_nuv)
            combined_nuv_ab_mag_f3.append(dataframe_nuvf3)
        if filters == 'F5':

            dataframe_nuvf5 = pd.DataFrame(mag_nuv)
            combined_nuv_ab_mag_f5.append(dataframe_nuvf5)
        
        if filters == 'F6':

            dataframe_nuvf6 = pd.DataFrame(mag_nuv)
            combined_nuv_ab_mag_f6.append(dataframe_nuvf6)


    
    #20.775


    # plt.figure(1)
    # print(len(combined_ab_magnitude_F1))
    if len(combined_ab_magnitude_F1) == 0:

        print('No image in this filter')
    else:
        dataframef1_99 = pd.concat(combined_ab_magnitude_F1)
        

        print(dataframef1_99)

        dataframe_F1 = dataframef1_99[dataframef1_99[0] != 99]
        print(dataframe_F1)

        number1,bins1,patches1 = plt.hist(dataframe_F1,bins=50,histtype='step',color='orange')
        # ##### to find the bin with the maximum value
        # bin_with_max_value_1 = bins1[np.argmax(number1)] 
        # print(f'Bin with the maximum source count in F1: {bin_with_max_value_1}')



    if len(combined_ab_magnitude_F2) == 0:

        print('No image in this filter')
    else:
        dataframef2_99 = pd.concat(combined_ab_magnitude_F2)
        

        print(dataframef2_99)

        dataframe_F2 = dataframef2_99[dataframef2_99[0] != 99]
        print(dataframe_F2)

        number2,bins2,patches2 = plt.hist(dataframe_F2,bins=50,histtype='step',color='gray')
        # ##### to find the bin with the maximum value
        # bin_with_max_value_2 = bins2[np.argmax(number2)] 
        # print(f'Bin with the maximum source count in F2: {bin_with_max_value_2}')

        


    if len(combined_ab_magnitude_F3) == 0:

        print('No image in this filter')
    else:
        dataframef3_99 = pd.concat(combined_ab_magnitude_F3)
        

        print(dataframef3_99)

        dataframe_F3 = dataframef3_99[dataframef3_99[0] != 99]
        print(dataframe_F3)

        number3,bins3,patches3 = plt.hist(dataframe_F3,bins=50,histtype='step',color='olive')
        # # ##### to find the bin with the maximum value
        # bin_with_max_value_3 = bins3[np.argmax(number3)] 
        # print(f'Bin with the maximum source count in F3: {bin_with_max_value_3}')


    
    if len(combined_ab_magnitude_F5) == 0:

        print('No image in this filter')
    else:
        dataframef5_99 = pd.concat(combined_ab_magnitude_F5)
        

        print(dataframef5_99)

        dataframe_F5 = dataframef5_99[dataframef5_99[0] != 99]
        print(dataframe_F5)

        number5,bins5,patches5 = plt.hist(dataframe_F5,bins=50,histtype='step',color='purple')
        # # ##### to find the bin with the maximum value
        # bin_with_max_value_5 = bins5[np.argmax(number5)] 
        # print(f'Bin with the maximum source count in F5: {bin_with_max_value_5}')




    if len(combined_ab_magnitude_F7) == 0:

        print('No image in this filter')
    else:
        dataframef7_99 = pd.concat(combined_ab_magnitude_F7)
        

        print(dataframef7_99)

        dataframe_F7 = dataframef7_99[dataframef7_99[0] != 99]
        print(dataframe_F7) # after removing 99 values

        number7,bins7,patches7 = plt.hist(dataframe_F7,bins=50,histtype='step',color='brown')
        # # ##### to find the bin with the maximum value
        # bin_with_max_value_7 = bins7[np.argmax(number7)] 
        # print(f'Bin with the maximum source count in F7: {bin_with_max_value_7}')
    


    
    # plt.xlabel('AB Magnitude',fontsize=12)
    # plt.ylabel('Source count',fontsize=12)
    # labels = ['F148W; CaF2-1','F154W; BaF2','F169M; Sapphire','F172M; Silica','F148Wa; CaF2-2']
    # plt.legend(labels,loc='upper left',fontsize='7')
    # plt.yscale('log') #log transformation
    # plt.show()



    #for nuv
        

    plt.figure(2)
    if len(combined_nuv_ab_mag_f1) == 0:
        print('No image in this filter')
    else:
        dataframe_nuvf1_99 = pd.concat(combined_nuv_ab_mag_f1)

        print(dataframe_nuvf1_99)

        dataframen1 = dataframe_nuvf1_99[dataframe_nuvf1_99[0] != 99]
        print(dataframen1)

        # nuvnumber1,nuvbins1,nuvpatches1 = plt.hist(dataframen1,bins=50,histtype='step',color='orange')
        # # ##### to find the bin with the maximum value
        # bin_with_max_value_n1 = nuvbins1[np.argmax(nuvnumber1)] 
        # print(f'Bin with the maximum source count in N1: {bin_with_max_value_n1}')


        # cutoff_point = nuvbins1[np.argmax(nuvnumber1)+5]
        # next_cutoffpoint = nuvbins1[np.argmax(nuvnumber1)+6]
        # final_cutoff = (cutoff_point+next_cutoffpoint)/2

        # print(final_cutoff)
        # plt.axvline(final_cutoff,c='black',linestyle='-')


    if len(combined_nuv_ab_mag_f2) == 0:
        print('No image in this filter')
    else:
        dataframe_nuvf2_99 = pd.concat(combined_nuv_ab_mag_f2)

        print(dataframe_nuvf2_99)


        dataframen2 = dataframe_nuvf2_99[dataframe_nuvf2_99[0] != 99]
        print(dataframen2)

        # nuvnumber2,nuvbins2,nuvpatches2 = plt.hist(dataframen2,bins=50,histtype='step',color='green')
        # # ##### to find the bin with the maximum value
        # bin_with_max_value_n2 = nuvbins2[np.argmax(nuvnumber2)] 
        # print(f'Bin with the maximum source count in N2: {bin_with_max_value_n2}')

    if len(combined_nuv_ab_mag_f3) == 0:
        print('No image in this filter')
    else:
        dataframe_nuvf3_99 = pd.concat(combined_nuv_ab_mag_f3)

        print(dataframe_nuvf3_99)

        dataframen3 = dataframe_nuvf3_99[dataframe_nuvf3_99[0] != 99]
        print(dataframen3)

        # nuvnumber3,nuvbins3,nuvpatches3 = plt.hist(dataframen3,bins=50,histtype='step',color='red')
        # # ##### to find the bin with the maximum value
        # bin_with_max_value_n3 = nuvbins3[np.argmax(nuvnumber3)] 
        # print(f'Bin with the maximum source count in N3: {bin_with_max_value_n3}')


    if len(combined_nuv_ab_mag_f5) == 0:
        print('No image in this filter')
    else:
        dataframe_nuvf5_99 = pd.concat(combined_nuv_ab_mag_f5)

        print(dataframe_nuvf5_99)

    
        dataframen5 = dataframe_nuvf5_99[dataframe_nuvf5_99[0] != 99]
        print(dataframen5)

        # nuvnumber5,nuvbins5,nuvpatches5 = plt.hist(dataframen5,bins=50,histtype='step',color='black')
        # # ##### to find the bin with the maximum value
        # bin_with_max_value_n5 = nuvbins5[np.argmax(nuvnumber5)] 
        # print(f'Bin with the maximum source count in N5: {bin_with_max_value_n5}')        

    if len(combined_nuv_ab_mag_f6) == 0:
        print('No image in this filter')
    else:
        dataframe_nuvf6_99 = pd.concat(combined_nuv_ab_mag_f6)

        print(dataframe_nuvf6_99)

        dataframen6 = dataframe_nuvf6_99[dataframe_nuvf6_99[0] != 99]
        print(dataframen6)

        # nuvnumber6,nuvbins6,nuvpatches6 = plt.hist(dataframen6,bins=50,histtype='step',color='blue')
        # # ##### to find the bin with the maximum value
        # bin_with_max_value_n6 = nuvbins6[np.argmax(nuvnumber6)] 
        # print(f'Bin with the maximum source count in N6: {bin_with_max_value_n6}')


    # print(combined_nuv_ab_mag_f5)
    # print(combined_nuv_ab_mag_f6)


    plt.xlabel('AB Magnitude',fontsize=12)
    plt.ylabel('Source count',fontsize=12)
    labels = ['N242W; Silica-1','N219M; NUVB15','N245M; NUVB13','N263M; NUVB4','N279N; NUVN2']
    plt.legend(labels,loc='upper left',fontsize='7')
    # plt.yscale('log') #log transformation
    # plt.show()

    # print(len( combined_ab_magnitude_F1 ))
    # print(len( combined_ab_magnitude_F2 ))
    # print(len( combined_ab_magnitude_F3 ))
    # print(len( combined_ab_magnitude_F5 ))
    # print(len( combined_ab_magnitude_F7 ))
        

    


    
global path_to_fuv
global path_to_nuv
global same_field_folder_path

path_to_nuv = '/Users/swagat98/Documents/Combine_cat/NUV-cat/' # path for the FUV files
path_to_fuv = '/Users/swagat98/Documents/Combine_cat/FUV-cat/' # path for the NUV files
same_field_folder_path = '/Users/swagat98/Documents/Combine_cat/Same_field/'  # path where we should create a new same field file
    
            

   

    

    
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



def source_dist_plot():
    path = '/Users/swagat98/Documents/Combine_cat/FUV-cat/'
    path2= '/Users/swagat98/Documents/Combine_cat/NUV-cat/'
    
    folder = glob(f'{path}*')
    print(len(folder))

    folder2 = glob(f'{path2}*')
    print(len(folder2))
    
    folder3 = glob(f'{same_field_folder_path}FUV_same_field/*')
    print(len(folder3))

    # print(folder3)
    dfs = []
    dfs2 = []
    dfs3 = []

    for files in range(len(folder)):
        openfiles = fits.open(folder[files])
        # openfiles.info()
        data = openfiles[1].data
        
        header = openfiles[1].header
        # print(header.keys)
        
        
        df = pd.DataFrame(data)
        df_new = pd.DataFrame()

    
        ra = df['ALPHA_J2000'].astype('float64')
        dec= df['DELTA_J2000'].astype('float64')
    
        gal_ra = df['GAL_RA'].astype('float64')
        gal_dec= df['GAL_DEC'].astype('float64')
    
        df_new['RA']=ra
        df_new['DEC']=dec
        df_new['GAL_RA']=gal_ra
        df_new['GAL_DEC']=gal_dec
        # print(df_new)
        
        dfs.append(df_new)
    
    for files2 in range(len(folder2)):
        openfiles2 = fits.open(folder2[files2])
        # openfiles2.info()
        data2 = openfiles2[1].data
        
        header2 = openfiles2[1].header
        # print(header.keys)
        
        
        df2 = pd.DataFrame(data2)
        df_new2 = pd.DataFrame()

    
        ra2 = df2['ALPHA_J2000'].astype('float64')
        dec2= df2['DELTA_J2000'].astype('float64')
    
        gal_ra2 = df2['GAL_RA'].astype('float64')
        gal_dec2= df2['GAL_DEC'].astype('float64')
    
        df_new2['RA']=ra2
        df_new2['DEC']=dec2
        df_new2['GAL_RA']=gal_ra2
        df_new2['GAL_DEC']=gal_dec2
        # print(df_new)
        
        dfs2.append(df_new2)
   
    for files3 in range(len(folder3)):

        openfiles3 = fits.open(folder3[files3])
        openfiles3.info()
        data3 = openfiles3[1].data
        header3 = openfiles3[1].header
        
        df3 = pd.DataFrame(data3)
        df_new3 = pd.DataFrame()

        ra3 = df3['ALPHA_J2000'].astype('float64')
        dec3 = df3['DELTA_J2000'].astype('float64')

        gal_ra3 = df3['GAL_RA'].astype('float64')
        gal_dec3 = df3['GAL_DEC'].astype('float64')

        df_new3['RA'] = ra3
        df_new3['DEC'] = dec3
        df_new3['GAL_RA'] = gal_ra3
        df_new3['GAL_DEC'] = gal_dec3

        dfs3.append(df_new3)


    com1 = pd.concat(dfs)
    com2 = pd.concat(dfs2)
    com3 = pd.concat(dfs3)
    
    combined_df = pd.concat([com1,com2])
    # print(combined_df)
        
    coords_1 = SkyCoord(combined_df['RA'],combined_df['DEC'],frame='icrs',unit=u.deg)
    # print(coords_1)
    # changed_to_gal = coords_1.transform_to('galactic')
    
 
    coords_2_1 = SkyCoord(com1['GAL_RA'],com1['GAL_DEC'],frame='galactic',unit=u.deg)
    coords_2_2 = SkyCoord(com2['GAL_RA'],com2['GAL_DEC'],frame='galactic',unit=u.deg)
    coords_2_3 = SkyCoord(com3['GAL_RA'],com3['GAL_DEC'],frame='galactic',unit=u.deg)

    # print(coords_2)
    # changed_to_eq = coords_2.transform_to('icrs')
    
    
    
    random_ra = np.random.randint(0,360,360)
    ran_dec   = np.zeros(360)
    cords  = list(zip(random_ra,ran_dec))
    # print(cords[:])
    galactic = SkyCoord(random_ra,ran_dec,unit=u.deg,frame='galactic')
    # print(galactic)

    equatorial = galactic.transform_to('icrs')
    # print('equatorial',equatorial)

    
    #####################
    
    plt.figure()
    plt.subplot(111, projection='aitoff')
    plt.grid(True)
    plt.xlabel('Galactic Longitude (l)',fontsize=11)
    plt.ylabel('Galactic Latitude (b)',fontsize=11)
    
    
    #these are in galactic l,b projection, galactic plane along 0deg
    # plt.scatter(galactic.l.wrap_at('180d').radian,galactic.b.radian,s=5)
    plt.scatter(coords_2_1.l.wrap_at('180d').radian,coords_2_1.b.radian,s=7,marker='+',color='blue')
    plt.scatter(coords_2_2.l.wrap_at('180d').radian,coords_2_2.b.radian,s=10,marker='+',color='green')
    plt.scatter(coords_2_3.l.wrap_at('180d').radian,coords_2_3.b.radian,s=7,marker='+',color='red')
    
    
    #these are in ra_dec projection, galactic plane along the curve.
    # plt.scatter(equatorial.ra.wrap_at('180d').radian,equatorial.dec.radian,s=5,c='red')
    # plt.scatter(coords_1.ra.wrap_at('180d').radian,coords_1.dec.radian,s=5,marker='+') 
    #both are projection of galactic plane. One is in icrs format and other is in galactic format

    # plt.show()
    




def mag_vs_magerr():
    
    path = '/Users/swagat98/Documents/Combine_cat/FUV-cat/*'
    image2 = glob(path)
    # print(image2)
    
    path2='/Users/swagat98/Documents/Combine_cat/NUV-cat/*'
    # image2 = glob(path2)
    
    
    
    combined_magnitude_and_error = []
    
            
    mag_and_error_f1 = []
    mag_and_error_f2 = []
    mag_and_error_f3 = []
    mag_and_error_f4 = []
    mag_and_error_f5 = []
    mag_and_error_f6 = []
    mag_and_error_f7 = []
        
     
    test_object = []
  
       
    for i in range(len(image2)):
        img2 = fits.open(image2[i])
        # print(img2)
        
        data = img2[1].data
        header = img2[1].header
   
   
        object_name = header['OBJECT']
   
        filters = header['FILTER']
        print(filters)     
        # print(header.keys)
        ra_point = [{keys:values} for keys,values in header.items() if keys == 'RA_PNT']
        dec_point = [{keys:values} for keys,values in header.items() if keys == 'DEC_PNT']

        ra_pointing = list(ra_point[0].values())[0]
        dec_pointing = list(dec_point[0].values())[0]
        
        print(f'Ra, Dec of the field is: {ra_pointing},{dec_pointing}')

        if filters == 'F1':
            mag = data['MAG_AUTO'].astype('float64')
            magerr = data['MAGERR_AUTO'].astype('float64')
 
        
            df1 = pd.DataFrame()
        
            df1['Magnitude'] = mag
            df1['Magnitude Error'] = magerr
        
            # print('dataframe',df1)
            df_new = df1
            
            mag_and_error_f1.append(df1)
            
        if filters == 'F2':
            mag = data['MAG_AUTO'].astype('float64')
            magerr = data['MAGERR_AUTO'].astype('float64')
 
        
            df2 = pd.DataFrame()
        
            df2['Magnitude'] = mag
            df2['Magnitude Error'] = magerr
        
            mag_and_error_f2.append(df2)

        if filters == 'F3':
            mag = data['MAG_AUTO'].astype('float64')
            magerr = data['MAGERR_AUTO'].astype('float64')
 

            df3 = pd.DataFrame()
        
            df3['Magnitude'] = mag
            df3['Magnitude Error'] = magerr
            
            
            mag_and_error_f3.append(df3)
            
           
            
            
        if filters == 'F4':
            mag = data['MAG_AUTO'].astype('float64')
            magerr = data['MAGERR_AUTO'].astype('float64')
 
        
            df4 = pd.DataFrame()
        
            df4['Magnitude'] = mag
            df4['Magnitude Error'] = magerr
        
            mag_and_error_f4.append(df4)
            
            
        if filters == 'F5':

            mag = data['MAG_AUTO'].astype('float64')
            magerr = data['MAGERR_AUTO'].astype('float64')
 
        
            df5 = pd.DataFrame()
        
            df5['Magnitude'] = mag
            df5['Magnitude Error'] = magerr
        
            mag_and_error_f5.append(df5)
                
        if filters == 'F6':
            mag = data['MAG_AUTO'].astype('float64')
            magerr = data['MAGERR_AUTO'].astype('float64')
 
        
            df6 = pd.DataFrame()
        
            df6['Magnitude'] = mag
            df6['Magnitude Error'] = magerr
        
            mag_and_error_f6.append(df6)

            
        if filters == 'F7':
            mag = data['MAG_AUTO'].astype('float64')
            magerr = data['MAGERR_AUTO'].astype('float64')
 
        
            df7 = pd.DataFrame()
        
            df7['Magnitude'] = mag
            df7['Magnitude Error'] = magerr
        
            mag_and_error_f7.append(df7)

    print(len(mag_and_error_f1))

    
        # labels = [,, ,]
    
    if len(mag_and_error_f1) == 0:
        print('No f1 info. ')
    else:
        df1_new = pd.concat(mag_and_error_f1)

        df1_list_sorted = df1_new.sort_values(by='Magnitude')
        print(df1_list_sorted)
        
        plt.figure(1)    

        df1_list_sorted[df1_list_sorted['Magnitude'] != 99]

        df1_mag = df1_list_sorted['Magnitude']
        df1_magerr = df1_list_sorted['Magnitude Error']
        number1 = plt.scatter(df1_mag,df1_magerr,color='green',s=1)
        plt.ylim(bottom=-0.1,top=0.5)
        plt.xlim(left=12,right=27)
        plt.xlabel('Magnitude',fontsize=12)
        plt.ylabel('Magnitude Error',fontsize=12)
        plt.legend(['N242W; Silica-1'],loc='upper left')
        
        num_of_f1_sources = len(df1_mag)
        print(f'Number of sources for F1 filter: {num_of_f1_sources}')
        
    if len(mag_and_error_f2) == 0:
        print('No f2 info. ')
    else:
        df2_new = pd.concat(mag_and_error_f2)

        df2_list_sorted = df2_new.sort_values(by='Magnitude')
        print(df2_list_sorted)
        
        plt.figure(2)
        df2_mag = df2_list_sorted['Magnitude']
        df2_maggerr = df2_list_sorted['Magnitude Error']
        number2 = plt.scatter(df2_mag,df2_maggerr,color='red',s=1)
        plt.ylim(bottom=-0.1,top=0.5)
        plt.xlim(left=12,right=27)
        plt.xlabel('Magnitude',fontsize=12)
        plt.ylabel('Magnitude Error',fontsize=12)
        plt.legend(['N219M; NUVB15'],loc='upper left')
       
        num_of_f2_sources = len(df2_mag)
        print(f'Number of sources for F2 filter: {num_of_f2_sources}')
        
         
    if len(mag_and_error_f3) == 0:
        print('No f3 info. ')
    else:
        df3_new = pd.concat(mag_and_error_f3)

        df3_list_sorted = df3_new.sort_values(by='Magnitude')
        print(df3_list_sorted)
        
        plt.figure(3)
        df3_mag = df3_list_sorted['Magnitude']
        df3_magerr = df3_list_sorted['Magnitude Error']
        number3 = plt.scatter(df3_mag,df3_magerr,color='blue',s=1)
        plt.ylim(bottom=-0.1,top=0.5)
        plt.xlim(left=12,right=27)
        plt.xlabel('Magnitude',fontsize=12)
        plt.ylabel('Magnitude Error',fontsize=12)
        plt.legend(['N245M; NUVB13'],loc='upper left')
        
        num_of_f3_sources = len(df3_mag)
        print(f'Number of sources for F3 filter: {num_of_f3_sources}')
        
        
    if len(mag_and_error_f5) == 0:
        print('No f5 info. ')
    else:
        df5_new = pd.concat(mag_and_error_f5)

        df5_list_sorted = df5_new.sort_values(by='Magnitude')
        print(df5_list_sorted)

        plt.figure(5)
        df5_mag = df5_list_sorted['Magnitude']
        df5_magerr = df5_list_sorted['Magnitude Error']
        number5 = plt.scatter(df5_mag,df5_magerr,color='black',s=1)
        plt.ylim(bottom=-0.1,top=0.5)
        plt.xlim(left=12,right=27)
        plt.xlabel('Magnitude',fontsize=12)
        plt.ylabel('Magnitude Error',fontsize=12)
        plt.legend(['N263M; NUVB4'],loc='upper left')

        num_of_f5_sources = len(df5_mag)
        print(f'Number of sources for F5 filter: {num_of_f5_sources}')

    if len(mag_and_error_f6) == 0:
        print('No f6 info. ')
    else:
        df6_new = pd.concat(mag_and_error_f6)

        df6_list_sorted = df6_new.sort_values(by='Magnitude')
        print(df6_list_sorted)

        plt.figure(6)
        df6_mag = df6_list_sorted['Magnitude']
        df6_magerr = df6_list_sorted['Magnitude Error']
        number6 = plt.scatter(df6_mag,df6_magerr,color='olive',s=1)
        plt.ylim(bottom=-0.1,top=0.5)
        plt.xlim(left=12,right=27)
        plt.xlabel('Magnitude',fontsize=12)
        plt.ylabel('Magnitude Error',fontsize=12)
        plt.legend(['N279N; NUVN2'],loc='upper left')

        num_of_f6_sources = len(df6_mag)
        print(f'Number of sources for F6 filter: {num_of_f6_sources}')
    
     
    if len(mag_and_error_f7) == 0:
        print('No f7 info. ')
    else:
        df7_new = pd.concat(mag_and_error_f7)

        df7_list_sorted = df7_new.sort_values(by='Magnitude')
        print(df7_list_sorted)

        plt.figure(7)
        df7_mag = df7_list_sorted['Magnitude']
        df7_magerr = df7_list_sorted['Magnitude Error']
        number7 = plt.scatter(df7_mag,df7_magerr,color='green',s=1)
        plt.ylim(bottom=-0.1,top=0.5)
        plt.xlim(left=12,right=27)
        plt.xlabel('Magnitude',fontsize=12)
        plt.ylabel('Magnitude Error',fontsize=12)
        plt.legend(['F148Wa; CaF2-2'],loc='upper left')

        num_of_f7_sources = len(df7_mag)
        print(f'Number of sources for F7 filter: {num_of_f7_sources}')
    
    
  
    plt.show()


# count_vs_magnitude_distribution()
# source_dist_plot()
# mag_vs_magerr()


def eff_area():

    fuv_path = '/Users/swagat98/UVIT-cat.fits'
    nuv_path = '/Users/swagat98/NUV effective area'

    file = fits.open(fuv_path)
    file.info

    data = file[1].data
    

    filedata = Table(data)
    df = filedata.to_pandas()

    f1_df = df[df['F1_MAG_AUTO'] > 0]

    
    mag = f1_df['F1_MAG_AUTO']
    magerr = f1_df['F1_MAGERR_AUTO']
 
    #to find the faintest magnitude which does not have error of .0198

    f1_errorcut = f1_df[f1_df['F1_MAGERR_AUTO'] < 0.198]
    # print(np.max(f1_errorcut['F1_MAG_AUTO']))

    # plt.scatter(mag,magerr,s=1)
    # plt.axhline(y = 0.1978, color = 'black', linestyle = '-') 
    # plt.ylim(-0.1,0.5)
    # plt.show()



    f2_df = df[df['N2_MAG_AUTO'] > 0]

    # print(f2_df[f2_df['F1_MAG_AUTO'] > 0])

    mag2 = f2_df['N2_MAG_AUTO']
    magerr2 = f2_df['N2_MAGERR_AUTO']

    print(np.min(magerr2))

    # print(mag2,magerr2)
    errorcut = f2_df[f2_df['N2_MAGERR_AUTO']<0.198]
    print(np.max(errorcut['N2_MAG_AUTO']))

    # print(np.max(f2_errorcut['F2_MAG_AUTO']))

    # plt.scatter(mag2,magerr2,s=1)
    # plt.axhline(y = 0.1978, color = 'black', linestyle = '-') 
    # plt.ylim(-0.1,0.5)
    # plt.show()
    

# eff_area()


# print(2.5*np.log10(1+(1/3)))

def positional_offset():

    path = '/Users/swagat98/UVIT Cross match with GAIA.fits'

    file = fits.open(path)
    
    data = file[1].data
    header = file[1].header

    print(header.keys)

    filedf = Table(data)
    df = filedf.to_pandas()

    number,bin,fig = plt.hist(df['angDist'],bins=25,histtype='step')
    plt.xlabel('Angular separation (in arcsec)',fontsize=12)
    plt.ylabel('Number of sources',fontsize=12)
    plt.show()

    print(f'Highest count: {bin[np.argmax(number)]}')

# positional_offset()

#galex : 0.6022544
#gaia : 0.24070276


# exposure time vs number counts of images plot
# take the median value of the exposure images. 
    # or separate the same into three or two different bins, one for low exposure time, one for 1000-2000 exposure - THese will give the bright magnitude limit
    # for faint sources, take images with higher exposure time. 
    # after separating the images into bins, find the limiting magnitude following the same procedure of magnitude vs magnitude error; 
    # i.e. an errorcut of 0.198 magnitude error.

def exposure_numbercounts_2():

    global bin1
    global bin2
    global bin3


    path = '/Users/swagat98/Documents/Combine_cat_2/FUV-cat/'
    path2= '/Users/swagat98/Documents/Combine_cat_2/NUV-cat/'

    files = glob( f'{ path }*.fits' )


    bin1 = []
    bin2 = []
    bin3 = []

    for i in range(len(files)):
        
        img = fits.open(files[i])

        data = img[1].data
        header = img[1].header

        
        datadf = Table(data)
        df = datadf.to_pandas()



        exposure_time = header['EXP_TIME'] 
    
        if exposure_time <= 2000:
            bin1.append(files[i])

        if exposure_time > 2000 and exposure_time < 10000:
            bin2.append(files[i])

        if exposure_time > 10000:
            bin3.append(files[i])


    
    
    print(f'Number of images in bin1: {len(bin1)}')
    print(f'Number of images in bin2: {len(bin2)}')
    print(f'Number of images in bin3: {len(bin3)}')


    exp_time_list_bin1 = []
    exp_time_list_bin2 = []
    exp_time_list_bin3 = []

    for len1 in range(len(bin1)):

        open_bin1 = fits.open(bin1[len1])

        data = open_bin1[1].data
        header = open_bin1[1].header


        opendf = Table(data)
        df_b1 = opendf.to_pandas()

        exp_time = header['EXP_TIME']
        exp_time_list_bin1.append(exp_time)


    #this is for bin1




# 327


    exp_time_list_bin2 = []

    for len2 in range(len(bin2)):

        open_bin2 = fits.open(bin2[len2])

        data = open_bin2[1].data
        header = open_bin2[1].header

        opendf = Table(data)
        df_b2 = opendf.to_pandas()

        exp_time = header['EXP_TIME']
        exp_time_list_bin2.append(exp_time)



    exp_time_list_b3 = []

    for len3 in range(len(bin3)):

        open_bin3 = fits.open(bin3[len3])

        data = open_bin3[1].data
        header = open_bin3[1].header

        opendf = Table(data)
        df_b3 = opendf.to_pandas()

        exp_time = header['EXP_TIME']
        exp_time_list_bin3.append(exp_time)

#
#for bin1

#     plt.figure(1,figsize=(10,6))
#     numb1, binb1, patchb1 = plt.hist(exp_time_list_bin1,bins=15,histtype='step',color='blue')
#     # plt.xlabel('Exposure time (sec)',fontsize=12)
#     # plt.ylabel('# of fields',fontsize=12)
#     # plt.show()
#     #
#     first_bin_b1 = binb1[np.argmax(numb1)]
#     second_bin_b1 = binb1[np.argmax(numb1)+1]
#     final_bin_b1 = (first_bin_b1 + second_bin_b1)/2
#     print(f'Highest number of images has exposure time: {final_bin_b1}, that are in bin1')
#     # print('Bins:',binb1)
#     # print(f'Number of b1 bins: {len(binb1)}')
# #     
# #
# #for bin2
# #
#     # plt.figure(2)
#     num2, bin2, patch2 = plt.hist(exp_time_list_bin2,bins=24,histtype='step',color='blue')
#     # plt.xlabel('Exposure time (sec)',fontsize=12)
#     # plt.ylabel('# of fields / 200 sec bin',fontsize=12)
#     # plt.show()

#     first_bin_b2 = bin2[np.argmax(num2)]
#     second_bin_b2 = bin2[np.argmax(num2)+1]
#     final_bin_b2 = (first_bin_b2 + second_bin_b2)/2 
#     print(f'Highest number of images has exposure time: {final_bin_b2}, that are in bin2')
#     # print('Bins: ',bin2)
#     # print('Number of bins: ',len(bin2))

#     plt.axvline(2000,linestyle='-',color='green',linewidth=2)
#     plt.axvline(10001,linestyle='-',color='red',linewidth=2)

# #
# # #3rd bin
# #
# #     plt.figure(3)
#     num3, bin3, patch3 = plt.hist(exp_time_list_bin3, bins=15, histtype='step',color='blue')
#     plt.xlabel('Exposure time (sec)',fontsize=12)
#     plt.ylabel('# of fields',fontsize=12)
#     plt.xscale('log')
#     plt.show()
# # #
#     first_bin_b3 = bin3[np.argmax(num3)]
#     second_bin_b3 = bin3[np.argmax(num3)+1]
#     final_bin_b3 = (first_bin_b3 + second_bin_b3)/2
#     print(f'Highest number of images has exposure time: {final_bin_b3}, that are in bin3')
    # print(f'Bins: ',bin3)
    # print('Number of b3 bins: ',len(bin3))




    #now separate the images in the respective bins and find the AB magnitude distributions

    nuvfiles = glob(f'{ path2 }*.fits')  


    nuvexp_time_list = []

    for j in nuvfiles:

        nuvopen = fits.open(j)

        nuvdata = nuvopen[1].data
        headernuv = nuvopen[1].header


        nuvexp_time = headernuv['EXP_TIME']
        nuvexp_time_list.append(nuvexp_time)


    # print(len(nuvexp_time_list))

                                                                                                                                                         
    # plt.figure(4)
    # nuvnum, nuvbin, nuvpatch = plt.hist(nuvexp_time_list,histtype='step',bins=45)
    # plt.xlim(0,5000)
    # plt.xlabel('Exposure time (sec)',fontsize=12)
    # plt.ylabel('# of images / 350 sec bin',fontsize=12)
    # plt.show()
    #
    # highest_nuvbin = nuvbin[np.argmax(nuvnum)]
    # second         = nuvbin[np.argmax(nuvnum)+1]
    # final = (highest_nuvbin + second)/2
    # print('final',final)
    # print(nuvbin)



exposure_numbercounts_2()


def bin_ab_mag_dist():


    f1_bin1 = []
    f2_bin1 = []
    f3_bin1 = []
    f5_bin1 = []
    f7_bin1 = []

    for l1 in bin1: 


        file1 = fits.open(l1)

        data1 = file1[1].data
        header1 = file1[1].header

        filter1 = header1['FILTER']

        data1df = Table(data1)
        df = data1df.to_pandas()

        if filter1 == 'F1':
            f1_bin1.append(df)
        if filter1== 'F2':
            f2_bin1.append(df)
        if filter1 == 'F3':
            f3_bin1.append(df)
        if filter1 == 'F5':
            f5_bin1.append(df)
        if filter1 == 'F7':
            f7_bin1.append(df)


    # plt.figure(1)
    # for f1 filter concat:

    if len(f1_bin1) == 0:
        print('No data in f1')
    else:
        f1b1_concat = pd.concat(f1_bin1)
        f1b1df = f1b1_concat[(f1b1_concat['MAG_AUTO']!=99) & (f1b1_concat['MAGERR_AUTO'] < 0.198)]

        # num1bin1, bin1bin1, patch1bin1 = plt.hist(f1b1df['MAG_AUTO'],histtype='step',bins = 50,color='black')


#         # sns.scatterplot(f1b1df,x='MAG_AUTO',y='MAGERR_AUTO',hue='red',size=1)
#         plt.figure(1)
        # plt.scatter(f1b1df['MAG_AUTO'],f1b1df['MAGERR_AUTO'],s=1,color='red')
        # plt.xlabel('Magnitude', fontsize=12)
        # plt.ylabel('Magnitude Error',fontsize=12)
        # plt.legend(['F148W; CaF2-1'],loc='upper left')
        # plt.show()



        print('Max f1 bin1',np.max(f1b1df['MAG_AUTO']))
        
        print('BIN1')
        # print(f'Highest count in magnitude f1: {bin1bin1[np.argmax(num1bin1)]}')

    if len(f2_bin1) == 0:
        print('No data in f2')
    else:
        f2b1_concat = pd.concat(f2_bin1)
        f2b1df = f2b1_concat[(f2b1_concat['MAG_AUTO']!=99) & (f2b1_concat['MAGERR_AUTO'] < 0.198)]
        
        # num2bin1,bin2bin1,patch2bin1 = plt.hist(f2b1df['MAG_AUTO'],histtype='step',bins = 50, color='green')
# #
        # plt.figure(2)
        # plt.scatter(f2b1df['MAG_AUTO'],f2b1df['MAGERR_AUTO'],s=1,color='blue')
        # plt.xlabel('Magnitude', fontsize=12)
        # plt.ylabel('Magnitude Error',fontsize=12)
        # plt.legend(['F154W; BaF2'],loc = 'upper left')
        # plt.show()
        # print(f'Highest count in magnitude f2: {bin2bin1[np.argmax(num2bin1)]}')

        print('Max f2 bin1: ',np.max(f2b1df['MAG_AUTO']))
        
    if len(f3_bin1) == 0:
        print('No data in f3')
    else:
        f3b1_concat = pd.concat(f3_bin1)
        f3b1df      = f3b1_concat[f3b1_concat['MAG_AUTO']!=99]
# 
        # num3bin1, bin3bin1, patch3bin1  = plt.hist(f3b1df['MAG_AUTO'],histtype='step',bins=50,color='purple')

        # print(f"Highest count in magnitude f3: {bin3bin1[np.argmax(num3bin1)]}")
 

    if len(f5_bin1) == 0:
        print('No data in f5')
    else:
        f5b1_concat = pd.concat(f5_bin1)
        f5b1df      = f5b1_concat[f5b1_concat['MAG_AUTO']!=99]

        # num5bin1, bin5bin1, patch5bin1 = plt.hist(f5b1df['MAG_AUTO'],histtype='step',bins=50,color='blue')

#         print(f"Highest count in magnitude F5: {bin5bin1[np.argmax(num5bin1)]}")

        


    if len(f7_bin1) == 0:
        print('No data in f7')
    else:
        f7b1_concat = pd.concat(f7_bin1)
        f7b1df      = f7b1_concat[f7b1_concat['MAG_AUTO']!=99]

        # num7bin1, bin7bin1, patch7bin1  = plt.hist(f7b1df['MAG_AUTO'],histtype='step',bins=50,color='red')

#         print(f"Highest count in magnitude F7: {bin7bin1[np.argmax(num7bin1)]}")


# #         # print(f'Highest source count in: {bin7bin1[np.argmax(num7bin1)]}')

#         fuv_legends =['F148W; CaF2-1','F154W; BaF2','F169M; Sapphire','F172M; Silica','F148Wa; CaF2-2']


#         plt.xlim(10,27)
#         plt.yscale('log')
#         plt.xlabel('AB Magnitude',fontsize=12)
#         plt.ylabel('Source count',fontsize=12)
#         plt.legend(fuv_legends,loc='upper left',fontsize=7)
#         plt.show()


    # print(len(f1_bin1))
    # print(len(f2_bin1))
    # print(len(f3_bin1))
    # print(len(f5_bin1))
    # print(len(f7_bin1))
#     # 

    
    f1_bin2 = []
    f2_bin2 = []
    f3_bin2 = []
    f5_bin2 = []
    f7_bin2 = []


# #     # print(bin2)
# #     # print(len(bin2))


    for l2 in bin2:

        file2 = fits.open(l2)

        data2 = file2[1].data
        header2 = file2[1].header

        filter2 = header2['FILTER']

        data2df = Table(data2)
        df      = data2df.to_pandas()

        if filter2 == 'F1':
            f1_bin2.append(df)
        if filter2 == 'F2':
            f2_bin2.append(df)
        if filter2 == 'F3':
            f3_bin2.append(df)
        if filter2 == 'F5':
            f5_bin2.append(df)
        if filter2 == 'F7':
            f7_bin2.append(df)

#     # plt.figure(2)

    if len(f1_bin2) == 0:
        print('No data in f1')
    else:
        f1b2_concat = pd.concat(f1_bin2)
        f1b2df      = f1b2_concat[(f1b2_concat['MAG_AUTO']!=99) & (f1b2_concat['MAGERR_AUTO']<0.198)]

        # num1bin2, bin1bin2, patch1bin2 = plt.hist(f1b2df['MAG_AUTO'],histtype='step',bins=50, color='black')

#         print('BIN 2')
#         print(f"Highest count in magnitude F1: {bin1bin2[np.argmax(num1bin2)]}")
# #         plt.figure(3)
        # plt.scatter(f1b2df['MAG_AUTO'],f1b2df['MAGERR_AUTO'],s=1,color='red')
        # plt.xlabel('Magnitude', fontsize=12)
        # plt.ylabel('Magnitude Error',fontsize=12)
        # plt.legend(['F148W; CaF2-1'],loc = 'upper left')
        # plt.show()

# #         print('Max f1 bin2: ',np.max(f1b2df['MAG_AUTO']))
    if len(f2_bin2) == 0:
        print('No data in f2')
    else:
        f2b2_concat = pd.concat(f2_bin2)
        f2b2df      = f2b2_concat[(f2b2_concat['MAG_AUTO']!=99) & ( f2b2_concat['MAGERR_AUTO']<0.198 )]

        # num2bin2, bin2bin2, patch2bin2 = plt.hist(f2b2df['MAG_AUTO'],histtype='step',bins=50, color='green')

#         print(f'Highest count in magnitude F2: {bin2bin2[np.argmax(num2bin2)]}')

# #         plt.figure(4)
        # plt.scatter(f2b2df['MAG_AUTO'],f2b2df['MAGERR_AUTO'],s=1,color='blue')
        # plt.xlabel('Magnitude', fontsize=12)
        # plt.ylabel('Magnitude Error',fontsize=12)
        # plt.legend(['F154W; BaF2'],loc = 'upper left')
        # plt.show()


# #         print('Max f2 bin2: ',np.max(f2b2df['MAG_AUTO']))


    if len(f3_bin2) == 0:
        print('No data in f3')
    else:
        f3b2_concat = pd.concat(f3_bin2)
        f3b2df      = f3b2_concat[f3b2_concat['MAG_AUTO']!=99]
# #
        # num3bin2, bin3bin2, patch3bin2 = plt.hist(f3b2df['MAG_AUTO'], histtype='step',bins=50,color='purple')

#         print(f'Highest count in magnitude F3: {bin3bin2[np.argmax(num3bin2)]}')
 

    if len(f5_bin2) == 0:
        print('No data in f5')
    else:
        f5b2_concat = pd.concat(f5_bin2)
        f5b2df      = f5b2_concat[f5b2_concat['MAG_AUTO']!=99]

        # num5bin2, bin5bin2, patch5bin2 = plt.hist(f5b2df['MAG_AUTO'],histtype='step',bins=50,color='blue')


#         print(f'Highest count in magnitude F5: {bin5bin2[np.argmax(num5bin2)]}')


    if len(f7_bin2) == 0:
        print('No data in f7')
    else:
        f7b2_concat = pd.concat(f7_bin2)
        f7b2df      = f7b2_concat[f7b2_concat['MAG_AUTO']!=99]

        # num7bin2, bin7bin2, patch7bin2 = plt.hist(f7b2df['MAG_AUTO'],histtype='step',bins=50, color='red')
#         print(f'Highest count in magnitude F7: {bin7bin2[np.argmax(num7bin2)]}')

# #         # print(f'Max value: {bin7bin2[np.argmax(num7bin2)]}')

#     plt.xlim(10,27)
#     plt.yscale('log')
#     plt.xlabel('AB Magnitude',fontsize=12)
#     plt.ylabel('Source count',fontsize=12)
#     plt.legend(fuv_legends,loc='upper left',fontsize=7)
#     plt.show()

    
#     # print(len(f1_bin2))
#     # print(len(f2_bin2))
#     # print(len(f3_bin2))
#     # print(len(f5_bin2))
#     # print(len(f7_bin2))






    f1_bin3 = []
    f2_bin3 = []
    f3_bin3 = []
    f5_bin3 = []
    f7_bin3 = []
    

    for l3 in bin3:
# # 
        file3 = fits.open(l3)

        data3 = file3[1].data
        header3 = file3[1].header

        filter3 = header3['FILTER']

        data3df = Table(data3)
        df = data3df.to_pandas()

        if filter3 == 'F1':
            f1_bin3.append(df)
        if filter3 == 'F2':
            f2_bin3.append(df)
        if filter3 == 'F3':
            f3_bin3.append(df)
        if filter3 == 'F5':
            f5_bin3.append(df)
        if filter3 == 'F7':
            f7_bin3.append(df)

    
# #     plt.figure(3)

    if len(f1_bin3)==0:
        print('No data in f1')
    else:
        f1b3_concat = pd.concat(f1_bin3)
        f1b3df      = f1b3_concat[( f1b3_concat['MAG_AUTO']!=99) & (f1b3_concat['MAGERR_AUTO']<0.198)]

        # num1bin3, bin1bin3, patch1bin3 = plt.hist(f1b3df['MAG_AUTO'],histtype='step',bins=50,color='black')

#         print('BIN 3')
#         print(f'Highest count in magnitude F1: {bin1bin3[np.argmax(num1bin3)]}')
# #         plt.figure(5)
        # plt.scatter(f1b3df['MAG_AUTO'],f1b3df['MAGERR_AUTO'],color='red',s=1)
        # plt.xlabel('Magnitude',fontsize=12)
        # plt.ylabel('Magnitude Error', fontsize=12)
        # plt.legend(['F148W; CaF2-1'],loc='upper left')
#         plt.show()

# # #         print('Max f1 bin3: ',np.max(f1b3df['MAG_AUTO']))


    if len(f2_bin3)==0:
        print('No data in f2')
    else:
        f2b3_concat = pd.concat(f2_bin3)
        f2b3df      = f2b3_concat[(f2b3_concat['MAG_AUTO']!=99) & (f2b3_concat['MAGERR_AUTO']<0.198)]

        # num2bin3, bin2bin3, patch2bin3 = plt.hist(f2b3df['MAG_AUTO'],histtype='step',bins=50,color='green')

#         print(f'Highest count in magnitude F2: {bin2bin3[np.argmax(num2bin3)]}')
# #         plt.figure(6)
        # plt.scatter(f2b3df['MAG_AUTO'],f2b3df['MAGERR_AUTO'],color='blue',s=1)
        # plt.xlabel('Magnitude',fontsize=12)
        # plt.ylabel('Magnitude Error', fontsize=12)
        # plt.legend(['F154W, BaF2'],loc='upper left')
        # plt.show()

# #         print('Max f2 bin3: ',np.max(f2b3df['MAG_AUTO']))
    if len(f3_bin3) == 0:
        print('No data in f3')
    else:
        f3b3_concat = pd.concat(f3_bin3)
        f3b3df      = f3b3_concat[f3b3_concat['MAG_AUTO']!=99]
# 
        # num3bin3, bin3bin3, patch3bin3 = plt.hist(f3b3df['MAG_AUTO'],histtype='step',bins=50,color='purple')

#         print(f'Highest count in magnitude F3: {bin3bin3[np.argmax(num3bin3)]}')


    if len(f5_bin3) == 0:
        print('No data in f5')
    else:
        f5b3_concat = pd.concat(f5_bin3)
        f5b3df      = f5b3_concat[f5b3_concat['MAG_AUTO']!=99]

        num5bin3, bin5bin3, patch5bin3 = plt.hist(f5b3df['MAG_AUTO'],histtype='step',bins=50,color='blue')
        # first_falloff = bin5bin3[np.argmax(num5bin3)+9]
        # second_falloff = bin5bin3[np.argmax(num5bin3)+10]
        # final_falloff = (first_falloff + second_falloff)/2
        # print(f'Magnitude limit: {final_falloff}')
        # plt.axvline(final_falloff,linestyle='-')
        # plt.show()

#         print(f'Highest count in magnitude F5: {bin5bin3[np.argmax(num5bin3)]}')

    if len(f7_bin3) == 0:
        print('No data in f7')
#     else:
#         f7b3_concat = pd.concat(f7_bin3)
#         f7b3df      = f7b3_concat[f7b3_concat['MAG_AUTO']!=99]

#         num7bin3, bin7bin3, patch7bin3 = plt.hist(f7b3df['MAG_AUTO'],histtype='step',bins=50,color='red')

#         print(f'Highest count in magnitude F7: {bin7bin3[np.argmax(num7bin3)]}')


#     plt.xlim(10,29)
#     plt.yscale('log')
#     plt.xlabel('AB Magnitude',fontsize=12)
#     plt.ylabel('Source count',fontsize=12)
#     plt.legend(fuv_legends,loc='upper left',fontsize=7)
#     plt.show()


bin_ab_mag_dist()

