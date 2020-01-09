import pandas as pd
import numpy as np
import sys
from numba import jit
import time
from datetime import date
import statistics as st
from math import sqrt
import os

print('Python Version ' + sys.version)
print('Pandas Version ' + pd.__version__)
print('Np Version ' + np.__version__)

root = 'J:\\TvZ projects\\Muscle\\DMi8\\Reanalysis\\20200109_C4_remaining\\' #remember to use double backslash or single forward slash (is a requirement for python convention)
input_file_H2AX = root + '20200109_H2AX_C4_remaining.xlsx'
input_file_TELO = root + '20200109_TELO_C4_remaining.xlsx'
input_file_DAPI = root + '20200109_DAPI_C4_remaining.xlsx'
output_details = '_threshold_1+.xlsx'
output_file = root + str(date.today()) + output_details
diagnostics_file = root + str(date.today()) + '_diagnostics' + output_details

#PARAMETERS FOR ANALYSIS
print("\n *** Reminder to change input and output values if needed. Check lines 23 - 31, and 55. Check readme for help. ***")
px = float(0.1) #float(input("Define pixel length or width in micron: \n"))n micron: \n")) 
px_area = float(px*px)
z_size = float(0.10) #float(input("Define size of z-step in mixron: \n"))
top_overlap_ratio = float(10) #input("Upper threshold for H2AX:TAF overlap ratio (> 1 means H2AX encompasses Telo, max 5): \n")
bottom_overlap_ratio = float(0.3) #input("Lower threshold for H2AX:TAF overlap ratio (0 means no overlap): \n")
TAF_size_threshold = float(0.4) #input("TAF size threshold in micron (0 will return no threshold): \n")
TAF_positive_threshold = float(1) #input("Number of TAF to qualify as senescence-positive: \n")
upper_TAF_positive_threshold = float(10000) #input("Maximum number of TAF per nucleus (to filter aberrant nuclei: \n")
H2AX_size_threshold = 0 #input("H2AX foci size threshold in micron (0 will return no threshold): \n")
TELO_size_threshold = 0 #input("Telo foci size threshold in micron (0 will return no threshold): \n")
start = time.time()
df_H2AX = pd.read_excel(input_file_H2AX, header=None)
df_TELO = pd.read_excel(input_file_TELO, header=None)
df_DAPI = pd.read_excel(input_file_DAPI, header=None)


dataset_H2AX = np.array(df_H2AX[0])
x_H2AX, y_H2AX, z_H2AX  = np.array(df_H2AX[3]), np.array(df_H2AX[4]), np.array(df_H2AX[5])
width_H2AX, height_H2AX, depth_H2AX = np.array(df_H2AX[6]), np.array(df_H2AX[7]), np.array(df_H2AX[8])
dataset_TELO = np.array(df_TELO[0])
x_TELO, y_TELO, z_TELO = np.array(df_TELO[3]), np.array(df_TELO[4]), np.array(df_TELO[5])
width_TELO, height_TELO, depth_TELO = np.array(df_TELO[6]), np.array(df_TELO[7]), np.array(df_TELO[8])
maxint_TELO = np.array(df_TELO[12])
avint_TELO = np.array(df_TELO[11])
dataset_DAPI = np.array(df_DAPI[0])
x_DAPI,y_DAPI, z_DAPI = np.array(df_DAPI[3]), np.array(df_DAPI[4]), np.array(df_DAPI[5])
width_DAPI, height_DAPI, depth_DAPI = np.array(df_DAPI[6]), np.array(df_DAPI[7]), np.array(df_DAPI[8])
interior_DAPI = np.array(df_DAPI[9])
#os.chdir("C:\\Users\\Ishaq\\Anaconda3\\") # change/activate this to your Anaconda3 directory - Anaconda prompt > "where Anaconda"

def floatify(val):
    if type(val) == int:
        return np.float32(val)
    else:
        return val
    
def positive(val):
    if val < 0:
        return val * -1
    else:
        return val
    
def nuclear_filter(p,DAPI1,DAPI2,DAPI3,DAPI4):
    if (floatify(p[0]) > floatify(DAPI1) #filters for within nuclear regions only
        and floatify(p[0]) < floatify(DAPI2) #dapi_vector0 < point0 < dapi_vector3
        and floatify(p[1]) > floatify(DAPI3) 
        and floatify(p[1]) < floatify(DAPI4)):
        return True
    else:
        return False

def convert_micron(pixel_size,V1):
    return floatify(V1*pixel_size)

def convert_size_micron(pixel_size,V1,V2):
    return floatify(V1*pixel_size)+floatify(V2*pixel_size)

def convert_start(pixel_size,centre,vector):
    return floatify((centre - (vector/2)) * pixel_size)
    
def convert_end(pixel_size,centre,vector):
    return floatify((centre + (vector/2)) * pixel_size)

@jit
def colocalisation(x1,y1,w1,h1,x2,y2,w2,h2):
    left = max(x1,x2)
    right = min(w1,w2)
    bottom = max(y1,y2)
    top = min(h1,h2)
    Aoverlap = (bottom - top)*(right - left)
    ratio_to_TELO = Aoverlap/((w1-x1)*(h1-y1))
    if not ((w1-x1)*(h1-y1)) == 0: #this bit was added to skirt an unknown bug
        ratio_to_TELO = Aoverlap/((w1-x1)*(h1-y1))
    #ratio_to_H2AX = Aoverlap/((w2-x2)*(h2-y2))
    if ratio_to_TELO > 0 and ratio_to_TELO < 5:
        return ratio_to_TELO
    else:
        return False

def full_analysis(index_1,index_2,index_3,index_4,index_5,index_6):
    all_H2AX = list(zip(x_H2AX, y_H2AX, z_H2AX, width_H2AX, height_H2AX, depth_H2AX))
    all_TELO = list(zip(x_TELO, y_TELO, z_TELO, width_TELO, height_TELO, depth_TELO, avint_TELO))
    all_DAPI = list(zip(x_DAPI, y_DAPI, z_DAPI, width_DAPI, height_DAPI, depth_DAPI))
    x_DAPI_start, y_DAPI_start, z_dim_DAPI = ([] for i in range(3))
    x_DAPI_end, y_DAPI_end,nuclear_count = ([] for i in range(3))
    xmicron_H2AX, ymicron_H2AX, zmicron_H2AX = ([] for i in range(3))
    filt_H2AX,filt_TELO,rellen_TELO = ([] for i in range(3))
    xmicron_H2AX_end, ymicron_H2AX_end, zmicron_H2AX_end = ([] for i in range(3))
    totaled_H2AX_count,totaled_TELO_count = [],[]
    xmicron_TELO, ymicron_TELO, zmicron_TELO = ([] for i in range(3))
    xmicron_TELO_end, ymicron_TELO_end, zmicron_TELO_end = ([] for i in range(3))
    xmicron_H2AX_start, ymicron_H2AX_start, zmicron_H2AX_start = [],[],[]
    xmicron_TELO_start, ymicron_TELO_start, zmicron_TELO_start = [],[],[]
    
    for point in all_DAPI[index_1:index_2]:
        x_DAPI_start.append(max(point[0] - (point[3]/2),0)) #x_dim_calc
        x_DAPI_end.append(point[0] + (point[3]/2))
        y_DAPI_start.append(max(point[1] - (point[4]/2),0)) #y_dim_calc
        y_DAPI_end.append(point[1] + (point[4]/2))
        z_dim_DAPI.append(float(1)) #z_dim_calc
    nuclear_count.append(len(x_DAPI_start))
    nuclear_count.remove(nuclear_count[0])
    all_DAPI = list(zip(x_DAPI[index_1:index_2], y_DAPI[index_1:index_2], 
                        z_DAPI[index_1:index_2], interior_DAPI[index_1:index_2], 
                        x_DAPI_start, x_DAPI_end, y_DAPI_start,y_DAPI_end,
                        z_dim_DAPI)) #one less list level to iterate
    
    #converts all pixels into microns, point comparisons in microns, point vs DAPI comparisons in vectors
    #px = pixel size
    for vector in all_H2AX[index_3:index_4]:
        xmicron_H2AX.append(convert_micron(px,vector[0]))
        ymicron_H2AX.append(convert_micron(px,vector[1]))
        zmicron_H2AX.append(convert_micron(px,vector[2]))
        xmicron_H2AX_start.append(convert_start(px,vector[0],vector[3]))
        ymicron_H2AX_start.append(convert_start(px,vector[1],vector[4]))
        zmicron_H2AX_start.append(convert_start(px,vector[2],vector[5]))            
        xmicron_H2AX_end.append(convert_end(px,vector[0],vector[3]))
        ymicron_H2AX_end.append(convert_end(px,vector[1],vector[4]))
        zmicron_H2AX_end.append(convert_end(px,vector[2],vector[5]))
    all_H2AX = list(zip(x_H2AX[index_3:index_4],y_H2AX[index_3:index_4],
                        z_H2AX[index_3:index_4],width_H2AX[index_3:index_4], 
                        height_H2AX[index_3:index_4],depth_H2AX[index_3:index_4], 
                        xmicron_H2AX, ymicron_H2AX, zmicron_H2AX,
                        xmicron_H2AX_start,ymicron_H2AX_start,zmicron_H2AX_end,
                        xmicron_H2AX_end,ymicron_H2AX_end,zmicron_H2AX_end))

    #converts all pixels into microns, point comparisons in microns, point vs DAPI comparisons in vectors
    for vector in all_TELO[index_5:index_6]:
        xmicron_TELO.append(convert_micron(px,vector[0]))
        ymicron_TELO.append(convert_micron(px,vector[1]))
        zmicron_TELO.append(convert_micron(px,vector[2]))
        xmicron_TELO_start.append(convert_start(px,vector[0],vector[3]))
        ymicron_TELO_start.append(convert_start(px,vector[1],vector[4]))
        zmicron_TELO_start.append(convert_start(px,vector[2],vector[5]))            
        xmicron_TELO_end.append(convert_end(px,vector[0],vector[3]))
        ymicron_TELO_end.append(convert_end(px,vector[1],vector[4]))
        zmicron_TELO_end.append(convert_end(px,vector[2],vector[5]))
    all_TELO = list(zip(x_TELO[index_5:index_6], y_TELO[index_5:index_6], 
                        z_TELO[index_5:index_6], avint_TELO[index_5:index_6]))
    
    maxsize_TELO = max(all_TELO[0])
    def telo_rellen(avint_TELO): #relative telo len in spreadsheet
        return avint_TELO/maxsize_TELO
    
    for vector in all_TELO:
        rellen_TELO.append(telo_rellen(vector[3]))
    all_TELO = list(zip(x_TELO[index_5:index_6],y_TELO[index_5:index_6],
                        z_TELO[index_5:index_6],width_TELO[index_5:index_6],
                        height_TELO[index_5:index_6],depth_TELO[index_5:index_6],
                        xmicron_TELO,ymicron_TELO,zmicron_TELO,xmicron_TELO_start,
                        ymicron_TELO_start,zmicron_TELO_start,
                        xmicron_TELO_end,ymicron_TELO_end,zmicron_TELO_end,
                        rellen_TELO,maxint_TELO[index_5:index_6],
                        avint_TELO[index_5:index_6]))
    
    #converts all pixels into microns, point comparisons in microns, point vs DAPI comparisons in vectors
    for DAPI_vectors in all_DAPI:
        for point in all_H2AX:
            if (nuclear_filter(point,DAPI_vectors[4],
                    DAPI_vectors[5],DAPI_vectors[6],DAPI_vectors[7]) == True):
                if (not H2AX_size_threshold == 0) and (point[9]       #filters by H2AX foci size
                     or point[10] or point[11] >= H2AX_size_threshold):
                    pass
                    filt_H2AX.append(point)
                elif H2AX_size_threshold == 0:
                    filt_H2AX.append(point)
            else:
                pass
        totaled_H2AX_count.append(len(filt_H2AX))
        for point in all_TELO:
            if (nuclear_filter(point,DAPI_vectors[4],
                    DAPI_vectors[5],DAPI_vectors[6],DAPI_vectors[7]) == True):
                if (not TELO_size_threshold == 0) and (point[9]       #filters by H2AX foci size
                     or point[10] or point[11] >= TELO_size_threshold):
                    pass
                    filt_H2AX.append(point)
                elif TELO_size_threshold == 0:
                    filt_TELO.append(point)
            else:
                pass
        totaled_TELO_count.append(len(filt_TELO))
    average_H2AX_count = len(filt_H2AX)/len(all_DAPI)

    dict_nuclei_H2AX, dict_nuclei_TELO, dict_interior = {},{},{}
    dict_H2AX_count, dict_TELO_count, dict_nuclei = {},{},{}
    for i in range(len(all_DAPI)):
        dict_nuclei["Nucleus no. " + str(i)] = all_DAPI[i][0:2]
        dict_nuclei_H2AX["Nucleus no. " + str(i)] = (filt_H2AX[totaled_H2AX_count[max(i-1,0)]:totaled_H2AX_count[i]])
        dict_nuclei_TELO["Nucleus no. " + str(i)] = (filt_TELO[totaled_TELO_count[max(i-1,0)]:totaled_TELO_count[i]])
        dict_H2AX_count["Nucleus no. " + str(i)] = len(filt_H2AX[totaled_H2AX_count[max(i-1,0)]:totaled_H2AX_count[i]])
        dict_TELO_count["Nucleus no. " + str(i)] = len(filt_TELO[totaled_TELO_count[max(i-1,0)]:totaled_TELO_count[i]])
        dict_interior["Nucleus no. " + str(i)] = all_DAPI[i][3]*px_area

    TTAF, TELO_len, HTAF, n_TAF = {},{},{},{} # TTAF = TAF with telomeres as comparisons
    percent_TELO_vs_TAF, percent_H2AX_vs_TAF = {}, {} # for 53BP1 diagnostic
    TAF_positive_nuclei, TAF_percent_positive, n_TAF_TELO, histo_nTAF = [],[],[],[]
    #interior_positive_nuclei, interior_negative_nuclei = {},{}
    _TELO_total, H2AX_total = [],[] # for percent TAF vs TELO or H2AX, 53BP1 diagnostic
    #colocalisation of H2AX and TELO
    z_stacks_per_TAF = float(TAF_size_threshold) / float(z_size)
    for (Tkey, Tval), (Hkey, Hval) in zip(dict_nuclei_TELO.items(),dict_nuclei_H2AX.items()):
        TAF_TELO,TAF_H2AX,TELO_length,H2AX_total = [],[],[],[]
        if Hkey == Tkey:
            for Tval2 in Tval:
                _TELO_total.append(Tval2[0:3])
                for Hval2 in Hval:
                    H2AX_total.append(Hval2[0:3])
                    coloc = colocalisation(Tval2[9],Tval2[10],Tval2[12],
                            Tval2[13],Hval2[9],Hval2[10],Hval2[12],Hval2[13])
                    if positive(Hval2[2] - Tval2[2]) > float(z_stacks_per_TAF):
                        pass
                    elif (coloc > float(bottom_overlap_ratio) and 
                          coloc < float(top_overlap_ratio) and
                          not Tval2[0:3] in TAF_TELO): # anti-ghosting included here
                        TAF_TELO.append(Tval2[0:3])
                        TELO_length.append(Tval2[15]) #15 is relative telomere length
                        TAF_H2AX.append(Hval2[0:3])
        n_TAF_TELO.append(len(TAF_TELO))
        if len(TAF_TELO) or len(_TELO_total) > 0:
            percent_TELO_vs_TAF[Tkey] = (len(TAF_TELO) / len(_TELO_total)) * 100
        if len(TAF_TELO) or len(H2AX_total) > 0:
            percent_H2AX_vs_TAF[Tkey] = (len(TAF_TELO) / len(H2AX_total)) * 100
        TTAF[Tkey] = TAF_TELO[:]
        TELO_len[Tkey] = TELO_length[:]
        HTAF[Tkey] = TAF_H2AX[:]
        n_TAF[Tkey] = len(TAF_TELO[:])
        histo_nTAF.append(len(TAF_TELO))
        if ((len(TTAF.get(Tkey)) < int(TAF_positive_threshold)) or 
            (len(TTAF.get(Tkey)) > int(upper_TAF_positive_threshold))):
            pass
        else:
            TAF_positive_nuclei.append("1")
        
    interior_negative_nuclei, interior_positive_nuclei = {},{}
    #all_nuclear_sizes = {}
    average_negative_nuclei, average_positive_nuclei = np.nan,np.nan
    average_interior = np.nan
    temp_neg, temp_pos, temp_all,av_int_all = [],[],[],[]
    for (Tkey, Tval), (Nkey, Nval) in zip(n_TAF.items(),dict_interior.items()):
        if Tkey == Nkey:
            #for (Tkey2,Tval2), (Nkey2,Nval2) in Tval,Nval:
            if Tval < int(TAF_positive_threshold) or Tval > int(upper_TAF_positive_threshold):
                temp_neg.append(Nval)
                temp_all.append(Nval)
                interior_negative_nuclei[Tkey] = Nval
            else:
                temp_pos.append(Nval)
                temp_all.append(Nval)
                interior_positive_nuclei[Tkey] = Nval
    if len(temp_neg) > 0:
        average_negative_nuclei = sum(temp_neg)/len(temp_neg)
    if len(temp_pos) > 0:
        average_positive_nuclei = sum(temp_pos)/len(temp_pos)
    average_interior = sum(temp_all)/len(temp_all)
    av_int_all.append(average_interior)
    av_int_all.append(average_positive_nuclei)
    av_int_all.append(average_negative_nuclei)

    nuclei_multidict = {}       #for checking overlapping nuclei and TAF
    for k,v in dict_nuclei.items():
        nuclei_multidict.setdefault(v,set().add(k))
        
    for k,v in nuclei_multidict.items():
        if not v is None and len(v) > 1:
            for k2,v2 in dict_nuclei.items():
                if v[0] == k2 and len(v) > 1:
                    dict_nuclei.pop('k2')
        else:
            pass

    if len(TTAF) > 0:
        TAF_percent_positive.append(len(TAF_positive_nuclei)/len(TTAF)*100) #percentage, count positive, count total
        TAF_percent_positive.append(len(TAF_positive_nuclei))
        TAF_percent_positive.append(len(TTAF))
    else:
        TAF_percent_positive.append(float(0))
        TAF_percent_positive.append(float(0))
        TAF_percent_positive.append(len(TTAF))
        
    return [TTAF, TELO_len, n_TAF, TAF_percent_positive, dict_nuclei, #0-4
            dict_H2AX_count, dict_TELO_count,average_H2AX_count,percent_TELO_vs_TAF,percent_H2AX_vs_TAF, #5-9
            dict_interior,histo_nTAF,av_int_all] #14-16
    
def sortby_treatment(dataset):
    obj_list = []
    for index, obj in enumerate(dataset):
        if index == 0:
            pass
        elif obj == dataset_H2AX[0]:
            pass
        elif obj not in obj_list:
            obj_list.append(obj)
        elif index+1 == len(dataset):
            break
    return obj_list

def treatment_index(dataset):
    obj_list, index_list = [],[]
    for index, obj in enumerate(dataset):
        if index == 0:
            pass
        elif obj == dataset_H2AX[0]:
            pass
        elif obj not in obj_list:
            obj_list.append(obj)
            index_list.append(max(0,index-2))
        elif index+1 == len(dataset):
            index_list.append(index)
            break
    return index_list

def retrieve_index(df):
    index_list = []
    for index, obj in enumerate(df):
        if index == 0:
            index_list.append(0)
        elif obj == x_H2AX[0]:
            index_list.append(index)
        elif index+1 == len(df):
            index_list.append(index)
            break
        else:
            pass
    return index_list

dataset_obj = sortby_treatment(dataset_DAPI)
dataset_indices = list(zip(treatment_index(dataset_DAPI),treatment_index(dataset_H2AX),treatment_index(dataset_TELO)))
image_indices = list(zip(retrieve_index(x_DAPI),retrieve_index(x_H2AX),retrieve_index(x_TELO)))
        
treatments_TTAF,treatments_pos,treatments_nTAF,treatments_Tlen = {},{},{},{}
treatments_nuclei,treatments_H2AX,treatments_TELO,treatments_nH2AX = {},{},{},{}
treatments_interior,treatments_avg_int = {},{}
percent_TELO_vs_TAF,percent_H2AX_vs_TAF = {},{}
histogram_nTAF = []
treatments_interior_neg,treatments_interior_pos = {},{}
treatments_avg_int_neg,treatments_avg_int_pos,treatments_avg_int = {},{},{}
for n, obj in enumerate(dataset_indices):
    images_TTAF, images_pos, images_Tlen, images_nTAF = {},{},{},{}
    images_nuclei, images_H2AX, images_TELO, images_nH2AX = {},{},{},{}
    images_TELO_vs_TAF, images_H2AX_vs_TAF, images_interior = {},{},{}
    images_avg_int = {}
    for m, obj_2 in enumerate(image_indices):
        if m > 0 and n <= len(dataset_indices):
            if ((image_indices[m-1][0] >= (dataset_indices[n-1][0])) and
                (image_indices[m][0] <= (dataset_indices[n][0])+1)):
                Image_num = "Image_" + str(m)
                analysis = full_analysis((image_indices[m-1][0]+1),image_indices[m][0],
                                     (image_indices[m-1][1]+1),image_indices[m][1],
                                     (image_indices[m-1][2]+1),image_indices[m][2])
                # output for analysis is a list of dicts
                images_TTAF[Image_num] = analysis[0] # sorts data by image
                images_Tlen[Image_num] = analysis[1]
                images_nTAF[Image_num] = analysis[2]
                images_pos[Image_num] = analysis[3]
                images_nuclei[Image_num] = analysis[4]
                images_H2AX[Image_num] = analysis[5]
                images_TELO[Image_num] = analysis[6]
                images_nH2AX[Image_num] = analysis[7]
                images_TELO_vs_TAF[Image_num] = analysis[8]
                images_H2AX_vs_TAF[Image_num] = analysis[9]
                #all interior values converted to micron
                images_interior[Image_num] = analysis[10]
                images_avg_int[Image_num] = analysis[12]
                #histogram_nTAF = analysis[11]
                treatments_TTAF.update({dataset_obj[n-1] : images_TTAF}) # sorts images by treatments
                treatments_pos.update({dataset_obj[n-1] : images_pos})
                treatments_nTAF.update({dataset_obj[n-1] : images_nTAF})
                treatments_Tlen.update({dataset_obj[n-1] : images_Tlen})
                treatments_nuclei.update({dataset_obj[n-1] : images_nuclei})
                treatments_H2AX.update({dataset_obj[n-1] : images_H2AX})
                treatments_TELO.update({dataset_obj[n-1] : images_TELO})
                treatments_nH2AX.update({dataset_obj[n-1] : images_nH2AX})
                treatments_interior.update({dataset_obj[n-1] : images_interior})
                percent_TELO_vs_TAF.update({dataset_obj[n-1] : images_TELO_vs_TAF})
                percent_H2AX_vs_TAF.update({dataset_obj[n-1] : images_H2AX_vs_TAF})
                #all interior values converted to micron
                treatments_avg_int.update({dataset_obj[n-1] : images_avg_int})
            else:
                pass

# =============================================================================
# for (treatments,images), (treatments2,images2) in zip(treatments_interior.items(),treatments_nTAF.items())):
#     for k,v in images.items():
#         print(k,v)
#         if treatments == dataset_obj[n]:
#             temp_percent.append(v[0])
# =============================================================================

def mean_percentage(treatments_list,list_position):
    mean_dict = {}
    temp_average = []
    for n,(treatments,images) in enumerate(treatments_list.items()):
        for k,v in images.items():
            if treatments == dataset_obj[n]:
                try:
                    temp_average.append(v[list_position])
                except:
                    temp_average.append(v)
        if len(temp_average) == len(images) and len(temp_average) > 1:
            temp_average = [st.mean(temp_average),st.stdev(temp_average),
                (st.stdev(temp_average)/sqrt(len(images))),len(images)]
            mean_dict.update({dataset_obj[n] : temp_average})
            temp_average = []
        else:
            mean_dict.update({dataset_obj[n] : float(0)})
            temp_average = []
    return mean_dict

mean_pos_treatment = mean_percentage(treatments_pos,0)
mean_interior_all = mean_percentage(treatments_avg_int,0)
mean_interior_pos = mean_percentage(treatments_avg_int,1)
mean_interior_neg = mean_percentage(treatments_avg_int,2)
mean_H2AX_pos = mean_percentage(treatments_nH2AX,0)

# # for use with multiple images per treatment. Comment out and use below if doing single image per treatment
# # percentage nuclei positive for senescence by treatment
dfmeanpos = pd.DataFrame.from_dict(mean_pos_treatment,orient='index',
        columns=['Percent TAF positive','SD','SE', 'n'])
# # average H2AX per nucleus
dfmeanH2AX = pd.DataFrame.from_dict(mean_H2AX_pos,orient='index',
        columns=['Average H2AX per nucleus', 'SD','SE','n'])
# # percentage nuclei positive for H2AX
dfmean_int_all = pd.DataFrame.from_dict(mean_interior_all,orient='index',
        columns=['Average all nuclear size (um)', 'SD','SE','n'])
dfmean_int_pos = pd.DataFrame.from_dict(mean_interior_pos,orient='index',
        columns=['Average TAF+ nuclear size (um)', 'SD','SE','n'])
dfmean_int_neg = pd.DataFrame.from_dict(mean_interior_neg,orient='index',
        columns=['Average TAF- nuclear size (um)', 'SD','SE','n'])
# =============================================================================
# dfH2AXpos = pd.DataFrame.from_dict(mean_pos_treatment,orient='index',
#         columns=['Percent H2AX positive','SD','SE', 'n'])
# =============================================================================

def df_from_dict(input_dict):
    return pd.DataFrame.from_dict({(i,j): input_dict[i][j]
                    for i in input_dict.keys()
                    for j in input_dict[i].keys()},
                    orient='index')

def dfs_tabs(df_list, sheet_list, file_name):
    writer = pd.ExcelWriter(file_name,engine='xlsxwriter')   
    for dataframe, sheet in zip(df_list, sheet_list):
        dataframe.to_excel(writer, sheet_name=sheet, startrow=0 , startcol=0)   
    writer.save()

# percentage of nuclei positive for senescence by images
dfpos = pd.DataFrame.from_dict({(i,j): treatments_pos[i][j]
        for i in treatments_pos.keys()
        for j in treatments_pos[i].keys()},
        orient='index', columns=['Percent TAF positive',
                                 'TAF positive count','Total nuclei count'])
dfnuclear_sizes = pd.DataFrame.from_dict({(i,j): treatments_avg_int[i][j]
        for i in treatments_avg_int.keys()
        for j in treatments_avg_int[i].keys()},        
        orient='index', columns=['Average all Nuclear Size (um)',
                                 'Average TAF+ Nuclear Size (um)','Average TAF- Nuclear Size (um)'])
dfTTAF = df_from_dict(treatments_TTAF)  # coordinate of each TAF positive telomere
dfTlen = df_from_dict(treatments_Tlen)  # relative length of each TAF positive telomere
dfnTAF = df_from_dict(treatments_nTAF)  # summarised number of TAF per nucleus

df_all_TAF = pd.concat([dfmeanpos,dfpos,dfmeanH2AX], axis=1, sort=False) #removed dfH2AXpos, needs rewriting
df_all_sizes = pd.concat([dfmean_int_all,dfmean_int_pos,dfmean_int_neg], axis=1, sort=True)
dfnuclei = df_from_dict(treatments_nuclei)  # position of nuclei
dfinterior = df_from_dict(treatments_interior)
dfnH2AX = df_from_dict(treatments_H2AX) # summarised number of H2AX per nucleus
dfnTELO = df_from_dict(treatments_TELO) # summarised number of TELO per nucleus
dfpposTELOvsTAF = df_from_dict(percent_TELO_vs_TAF)# % 53BP1 positive vs correlate
dfpposH2AXvTAF = df_from_dict(percent_H2AX_vs_TAF)  # % H2AX positive vs total correlate        
dfnTAF_hist = pd.DataFrame(histogram_nTAF)                           


df_main,sheets_main = ([df_all_TAF,df_all_sizes,dfnTAF],
["Percent TAF positive","Nuclear size summary","nTAF per nucleus"])
df_diagnostics, sheets_diagnostics = ([dfTTAF,dfTlen,dfnuclear_sizes,dfnuclei,dfinterior,dfnH2AX,dfnTELO],
["TAF coordinates","Relative Telo length","Nuclear size (img)","Nuclear coordinates","Raw Nuclear Size","H2A.X count","Telo count"])
dfs_tabs(df_main,sheets_main,output_file)
dfs_tabs(df_diagnostics,sheets_diagnostics,diagnostics_file)
end = time.time()
print("Runtime = %s" % (end - start) + " s")