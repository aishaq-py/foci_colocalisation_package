import pandas as pd
import numpy as np
import statistics as stats
import math

root = 'H:\\Abbas\\2018\\Deltapix backup 26062018\\15112018 3t3l1 droplet size r3\\'
file_name = '20181128 LD sizes R3'
_format = '.xlsx'
input_file = file_name + _format
output_file = file_name + ' output' + '.csv'

df_green = pd.read_excel(input_file, header=None)
df_all = np.array(df_green)
df_dataset = np.array(df_green[0])
df_interior = np.array(df_green[4])
df_round = np.array(df_green[6])
#useful values: 4=interior, 6=roundness

# output = '1' for descriptive stats (Mean, SD, SE, n),
# '2' for descriptive stats of ALL DATA on sheet:
# (Mean, SD, SE - ALL, nALL, SE - cells, nCells by header)
# Do not modify this function
def summarise(treatments,df_interior,df_round,output):
    all_obj, desc_stats = ({} for i in range(2))
    obj_list, all_desc = ([] for i in range(2))
    num = 0
    for index, obj in enumerate(df_interior):
        if index == 0:
            pass
        elif obj == df_interior[0]:
            all_obj[treatments[num]] = obj_list
            desc_stats[treatments[num]] = (stats.mean(obj_list), 
                      stats.stdev(obj_list),
                      (stats.stdev(obj_list)/math.sqrt(len(obj_list))),
                      len(obj_list))
            num += 1
            obj_list = []
        elif index+1 == len(df_interior):
            obj_list.append(obj)
            all_desc.append(obj)
            all_obj[treatments[num]] = obj_list
            all_desc = [stats.mean(all_desc), 
                              stats.stdev(all_desc), 
                              stats.stdev(all_desc)/math.sqrt(len(all_desc)),
                              len(all_desc),
                              (stats.stdev(obj_list)/math.sqrt(len(all_obj))),
                              len(all_obj)]
        elif df_round[index] >= float(50.0):
            obj_list.append(obj)
            all_desc.append(obj)
        else:
            pass
    if output == '1':
        return desc_stats
    elif output == '2':
        return all_desc       

def sortby_treatment(dataset):
    obj_list, index_list = [],[]
    for index, obj in enumerate(dataset):
        if index == 0:
            pass
        elif obj == dataset[0]:
            pass
        elif obj not in obj_list:
            obj_list.append(obj)
            index_list.append(max(0,index-2))
        elif index+1 == len(dataset):
            index_list.append(index)
            break
    return obj_list

treatments = sortby_treatment(df_dataset)
output_1 = summarise(treatments,df_interior,df_round,'1')
output_2 = summarise(treatments,df_interior,df_round,'2')

# =============================================================================
# def summarise_all(df):
#     df_sum = {}
#     for array in df.T:
#         try:
#             x = summarise(array[len(array)],'2')
#             df_sum[array[0]] = x
#         except Exception:
#             print("ASD")
# =============================================================================

desc = summarise(df_round,'1') # change df_round to df_whatever you like - change value for df_whatever above
all_desc = summarise(df_round,'2') 
# run multiple summarise functions to get everything right now, 
# working on a solution to iterate over entire excel sheet
#ass = summarise_all(df_green)
#fsad = summarise_all(df_all)

dfdesc = pd.DataFrame.from_dict(desc, orient='index', 
                                columns=['Mean','SD','SE','n'])
dfall_desc = pd.DataFrame(all_desc, index=['Mean','SD all','SE - all',
                                           'n - all','SE by cell','n by cell'])
cleanup = pd.concat([dfall_desc,dfdesc])
cleanup.to_csv(output_file, index=True)