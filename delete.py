
# #                     # print(df.columns)
#                     mean_num_f1 = 0
#                     mean_num_f2 = 0
#                     mean_num_f3 = 0
#                     mean_num_f5 = 0
#                     mean_num_f7 = 0

                    
#                     count_f1 = 0
#                     count_f2 = 0
#                     count_f3 = 0
#                     count_f5 = 0
#                     count_f7 = 0


#                     # print(matched_rows)
#                     for j in range(len(index_list)):
#                         # print(j)
                        
#                         try:
#                             if drop_rows['F1_M_A'].iloc[j] > 0:
#                                 mean_num_f1 += drop_rows['F1_M_A'].iloc[j]
#                                 count_f1 += 1
#                         except:
#                             print('F1 not present')
#                         try:
#                             if drop_rows['F2_M_A'].iloc[j]> 0:
#                                 mean_num_f2 += drop_rows['F2_M_A'].iloc[j]
#                                 count_f2 += 1
#                         except:
#                             print('F2 not present')
                        
#                         try:
#                             if drop_rows['F3_M_A'].iloc[j] > 0:
#                                 mean_num_f3 += drop_rows['F3_M_A'].iloc[j]

#                                 count_f3 += 1
#                         except:
#                             print('F3 data not present')

#                         try:
#                             if drop_rows['F5_M_A'].iloc[j] > 0:
#                                 mean_num_f5 += drop_rows['F5_M_A'].iloc[j]
#                                 count_f5 += 1
#                         except:
#                             print('F5 data not present')

#                         try:
#                             if drop_rows['F7_M_A'].iloc[j]>0:
#                                 mean_num_f7 += drop_rows['F7_M_A'].iloc[j]
#                                 count_f7 += 1
#                         except:
#                             print('F7 not present')

#                     print(mean_num_f1)
#                     print(mean_num_f2)
#                     print(mean_num_f3)
#                     print(mean_num_f5)
#                     print(mean_num_f7)

#                     # print(mean_num)
#                     print('Count')

#                     try:

#                         print(mean_num_f1/count_f1)

#                         f1_mag_auto = mean_num_f1/count_f1
#                         drop_rows.loc[0,'F1_M_A'] = f1_mag_auto
                        
                        
#                     except:
#                         print('No f1 value')

#                     try:
#                         print(mean_num_f2/count_f2)

#                         f2_mag_auto = mean_num_f2/count_f2
#                         drop_rows.loc[0,'F2_M_A'] = f2_mag_auto
#                     except:
#                         print('No f2 value')
#                     try:
#                         print(mean_num_f3/count_f3)

#                         f3_mag_auto = mean_num_f3/count_f3
#                         drop_rows.loc[0,'F3_M_A'] = f3_mag_auto
#                     except:
#                         print('No f3 value')
#                     try:
#                         print(mean_num_f5/count_f5)

#                         f5_mag_auto = mean_num_f5/count_f5
#                         drop_rows.loc[0,'F5_M_A'] = f5_mag_auto
#                     except:
#                         print('No f5 values')
#                     try:
#                         print(mean_num_f7/count_f7)

#                         f7_mag_auto = mean_num_f7/count_f7
#                         drop_rows.loc[0,'F7_M_A'] = f7_mag_auto
#                     except:
#                         print('No f7 values') 





# # # ### match rows
# # # ### Doing iteratively for every matched rows
#                 if len(index_list) == 1:
#                     pass
#                 elif len(index_list) > 1:

#                     df.drop(index_list,axis='rows',inplace=True)
#                     df.reset_index(drop=True,inplace=True)
#                     print(df.shape)
#                     print(df)
# # match found except the first row
#                     replace_rows = drop_rows.iloc[0].to_frame().T
# #     # #                     # print(replace_rows['F1_M_A'])


#                     print(df.shape)
#                     df = pd.concat([df,replace_rows],ignore_index=True)
#                     print(df)
#                     # print(df)
#                     print(df.shape)

#                     print('\n')

import pandas as pd

# Sample DataFrame
data = {'A': [1, 2, 3, 4], 'B': [5, 6, 7, 8]}
df = pd.DataFrame(data)
print(df)

# Using iterrows to iterate through the DataFrame
for index, row in df.iterrows():
    try:
        value_a = row['A']
        value_b = row['B']
        print(f"Row {index}: A = {value_a}, B = {value_b}")
    except:
        print('Row not present')

    df.drop(3,axis='rows',inplace=True)
    df.reset_index(drop=True,inplace=True)