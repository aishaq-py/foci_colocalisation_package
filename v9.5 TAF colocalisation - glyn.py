import pandas as pd
import numpy as np
import sys
from numba import jit
import time
from datetime import date
import statistics as st
from math import sqrt

print('Python Version ' + sys.version)
print('Pandas Version ' + pd.__version__)
print('Np Version ' + np.__version__)

root = 'C:\\Users\\Ishaq\\Documents\\All spreadsheets\\Glyn\\' #remember to use double backslash or single forward slash (is a requirement for python convention)
input_file_H2AX = root + 'HDFn434 H2AX.xlsx'
input_file_TELO = root + 'HDFn434 TELO.xlsx'
input_file_DAPI = root + 'HDFn434 DAPI.xlsx'
output_file = root + str(date.today()) + '_HDFn434_after_script.xlsx'

#parameters for analysis - change to absolute values if needed as such
print("\n Input all requested values as decimal numbers (floats).")
px = float(0.065) #float(input("Define pixel length or width in micron: \n")) 
z_size = float(0.18) #float(input("Define size of z-step in mixron: \n"))
top_overlap_ratio = float(1.5) #input("Upper threshold for H2AX:TAF overlap ratio (> 1 means H2AX encompasses Telo, max 5): \n")
bottom_overlap_ratio = float(0.3) #input("Lower threshold for H2AX:TAF overlap ratio (0 means no overlap): \n")
TAF_size_threshold = float(0.5) #input("TAF size threshold in micron (0 will return no threshold): \n")
TAF_positive_threshold = float(2) #input("Number of TAF to qualify as senescence-positive: \n")
upper_TAF_positive_threshold = float(10) #input("Maximum number of TAF per nucleus (to filter aberrant nuclei: \n")
H2AX_size_threshold = 0 #input("H2AX foci size threshold in micron (0 will return no threshold): \n")
TELO_size_threshold = 0 #input("Telo foci size threshold in micron (0 will return no threshold): \n")

start = time.time()

df_H2AX = pd.read_excel(input_file_H2AX, header=None)
df_TELO = pd.read_excel(input_file_TELO, header=None)
df_DAPI = pd.read_excel(input_file_DAPI, header=None)

#x_H2AX = pd.to_numeric(df_H2AX[3], errors='coerce')
#based on my mod of Antho's Icy protocol, df 0-Full path, 1-Parent Folder, 2- Dataset, 3- ROI
#4- Colour, 5- X, 6- Y, 7- Z, 8- Width, 9- Height, 10- Depth, 11-Contour, 12- Interior,
#13- Sphericity, 14- Roundness, 15- Convexity, 16- Min intensity, 17- Avg. intensity,
#18- Max intensity, 19- Sum Intensity, 20- Std dev Intensity
#channel 2 for red
dataset_H2AX = np.array(df_H2AX[0])
x_H2AX, y_H2AX, z_H2AX  = np.array(df_H2AX[3]), np.array(df_H2AX[4]), np.array(df_H2AX[5])
width_H2AX, height_H2AX, depth_H2AX = np.array(df_H2AX[6]), np.array(df_H2AX[7]), np.array(df_H2AX[8])
dataset_TELO = np.array(df_TELO[0])
x_TELO, y_TELO, z_TELO = np.array(df_TELO[3]), np.array(df_TELO[4]), np.array(df_TELO[5])
width_TELO, height_TELO, depth_TELO = np.array(df_TELO[6]), np.array(df_TELO[7]), np.array(df_TELO[8])
maxint_TELO = np.array(df_TELO[21])
avint_TELO = np.array(df_TELO[20])
dataset_DAPI = np.array(df_DAPI[0])
x_DAPI,y_DAPI, z_DAPI = np.array(df_DAPI[3]), np.array(df_DAPI[4]), np.array(df_DAPI[5])
width_DAPI, height_DAPI, depth_DAPI = np.array(df_DAPI[6]), np.array(df_DAPI[7]), np.array(df_DAPI[8])

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
    x_DAPI_start, y_DAPI_start, z_dim_DAPI, ROI_end_DAPI = ([] for i in range(4))
    x_DAPI_end, y_DAPI_end,nuclear_count = ([] for i in range(3))
    x_dim_H2AX, y_dim_H2AX, z_dim_H2AX = ([] for i in range(3))
    xmicron_H2AX, ymicron_H2AX, zmicron_H2AX = ([] for i in range(3))
    filt_H2AX,filt_TELO,rellen_TELO = ([] for i in range(3))
    xmicron_H2AX_end, ymicron_H2AX_end, zmicron_H2AX_end = ([] for i in range(3))
    totaled_H2AX_count, H2AX_count, H2AX_volume = ([] for i in range(3))
    x_dim_TELO, y_dim_TELO, z_dim_TELO = ([] for i in range(3))
    xmicron_TELO, ymicron_TELO, zmicron_TELO = ([] for i in range(3))
    sxmicron_TELO, symicron_TELO, szmicron_TELO = ([] for i in range(3))
    xmicron_TELO_end, ymicron_TELO_end, zmicron_TELO_end = ([] for i in range(3))
    totaled_TELO_count, TELO_count, TELO_volume = ([] for i in range(3))
    ROI_end_TELO, IMG_no_TELO, values = ([] for i in range(3))
    xmicron_H2AX_start, ymicron_H2AX_start, zmicron_H2AX_start = [],[],[]
    xmicron_TELO_start, ymicron_TELO_start, zmicron_TELO_start = [],[],[]

    for point in all_DAPI[index_1:index_2]:
        if point[0] == all_DAPI[0][0]: #numbers the images
            x_DAPI_start.append('Position X start')
            x_DAPI_end.append('Position X end')
            y_DAPI_start.append('Position Y start')
            y_DAPI_end.append('Position Y end')
            z_dim_DAPI.append('Position Z')
            nuclear_count.append(len(x_DAPI_start)-1) #first iter will show 0, last iter outside loop, nuclei = total so far. Do subtractions
        else:
            x_DAPI_start.append(max(point[0] - (point[3]/2),0)) #x_dim_calc
            x_DAPI_end.append(point[0] + (point[3]/2))
            y_DAPI_start.append(max(point[1] - (point[4]/2),0)) #y_dim_calc
            y_DAPI_end.append(point[1] + (point[4]/2))
            z_dim_DAPI.append(float(1)) #z_dim_calc
    nuclear_count.append(len(x_DAPI_start))
    nuclear_count.remove(nuclear_count[0])
    all_DAPI = list(zip(x_DAPI[index_1:index_2], y_DAPI[index_1:index_2], 
                        z_DAPI[index_1:index_2], x_DAPI_start,
                        x_DAPI_end, y_DAPI_start,y_DAPI_end,
                        z_dim_DAPI)) #one less list level to iterate
    all_DAPI.remove(all_DAPI[0])

    #converts all pixels into microns, point comparisons in microns, point vs DAPI comparisons in vectors
    #px = pixel size
    for vector in all_H2AX[index_3:index_4]:
        if vector[0] == all_H2AX[0][0]:
            xmicron_H2AX.append("X") #converts x-val into microns (max is 2080x0.16)
            ymicron_H2AX.append("Y")
            zmicron_H2AX.append("Z")
            xmicron_H2AX_start.append("sX") #converts box width values into microns
            ymicron_H2AX_start.append("sY")
            zmicron_H2AX_start.append("sZ")
            xmicron_H2AX_end.append("X end") #adds box width in microns to x-val - calculates len of box side
            ymicron_H2AX_end.append("Y end")
            zmicron_H2AX_end.append("Z end")
        else:
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
    all_H2AX.remove(all_H2AX[0])
    
    #converts all pixels into microns, point comparisons in microns, point vs DAPI comparisons in vectors
    for vector in all_TELO[index_5:index_6]:
        if vector[0] == all_TELO[0][0]:
            xmicron_TELO.append("X") #converts x-val into microns (max is 2080x0.16)
            ymicron_TELO.append("Y")
            zmicron_TELO.append("Z")
            xmicron_TELO_start.append("sX") #converts box width values into microns
            ymicron_TELO_start.append("sY")
            xmicron_TELO_start.append("sZ")
            xmicron_TELO_end.append("X end") #defines rightmost edges of the foci
            ymicron_TELO_end.append("Y end")
            zmicron_TELO_end.append("Z end")
        else:
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
    
    maxsize_TELO = max(all_TELO[17])
    def telo_rellen(avint_TELO): #relative telo len in spreadsheet
        return avint_TELO/maxsize_TELO
    
    for vector in all_TELO:
        if vector[0] == all_TELO[0][0]:
            rellen_TELO.append("Relative Telo Length")
        else:
            rellen_TELO.append(telo_rellen(vector[3]))
    all_TELO = list(zip(x_TELO[index_5:index_6],y_TELO[index_5:index_6],
                        z_TELO[index_5:index_6],width_TELO[index_5:index_6],
                        height_TELO[index_5:index_6],depth_TELO[index_5:index_6],
                        xmicron_TELO,ymicron_TELO,zmicron_TELO,xmicron_TELO_start,
                        ymicron_TELO_start,zmicron_TELO_start,
                        xmicron_TELO_end,ymicron_TELO_end,zmicron_TELO_end,
                        rellen_TELO,maxint_TELO[index_5:index_6],
                        avint_TELO[index_5:index_6]))
    all_TELO.remove(all_TELO[0])
    
    #converts all pixels into microns, point comparisons in microns, point vs DAPI comparisons in vectors
    num = 0
    for DAPI_vectors in all_DAPI:
        for point in all_H2AX:
            if point[0] == all_H2AX[0][0]:
                totaled_H2AX_count.append(len(filt_H2AX))
            elif DAPI_vectors[0] == all_DAPI[0][0]:
                pass #dapi_vector0 < point0 < dapi_vector3
            elif (nuclear_filter(point,DAPI_vectors[3],
                    DAPI_vectors[4],DAPI_vectors[5],DAPI_vectors[6]) == True):
                if (not H2AX_size_threshold == 0) and (point[9]       #filters by H2AX foci size
                     or point[10] or point[11] >= H2AX_size_threshold):
                    pass
                    filt_H2AX.append(point)
                elif H2AX_size_threshold == 0:
                    filt_H2AX.append(point)
            else:
                pass
        
    for DAPI_vectors in all_DAPI:
        for point in all_TELO:
            if point[0] == all_TELO[0][0]:
                totaled_TELO_count.append(len(filt_TELO))
            elif DAPI_vectors[0] == all_DAPI[0][0]:
                pass 
            elif (nuclear_filter(point,DAPI_vectors[3],
                    DAPI_vectors[4],DAPI_vectors[5],DAPI_vectors[6]) == True):
                if (not TELO_size_threshold == 0) and (point[9]       #filters by H2AX foci size
                     or point[10] or point[11] >= TELO_size_threshold):
                    pass
                    filt_H2AX.append(point)
                elif TELO_size_threshold == 0:
                    filt_TELO.append(point)
            else:
                pass

    dict_nuclei_H2AX, dict_nuclei_TELO = {},{}
    dict_H2AX_count, dict_TELO_count, dict_nuclei = {},{},{}
    for i in range(len(all_DAPI)):
        dict_nuclei["Nucleus no. " + str(i)] = all_DAPI[i][0:2]
        dict_nuclei_H2AX["Nucleus no. " + str(i)] = (filt_H2AX[totaled_H2AX_count[max(i-1,0)]:totaled_H2AX_count[i]])
        dict_nuclei_TELO["Nucleus no. " + str(i)] = (filt_TELO[totaled_TELO_count[max(i-1,0)]:totaled_TELO_count[i]])
        dict_H2AX_count["Nucleus no. " + str(i)] = len((filt_H2AX[totaled_H2AX_count[max(i-1,0)]:totaled_H2AX_count[i]]))
        dict_TELO_count["Nucleus no. " + str(i)] = len((filt_TELO[totaled_TELO_count[max(i-1,0)]:totaled_TELO_count[i]]))
    TTAF, TELO_len, HTAF, n_TAF = {},{},{},{}
    TAF_TELO, TAF_H2AX, TELO_length, n_TAF_TELO = [],[],[],[]
    TAF_positive_nuclei, TAF_percent_positive = [],[]
    num = 0
    #colocalisation of H2AX and TELO
    z_stacks_per_TAF = float(TAF_size_threshold) / float(z_size)
    for (Tkey, Tval), (Hkey, Hval) in zip(dict_nuclei_TELO.items(),dict_nuclei_H2AX.items()):
        if Hkey == Tkey:
            for Tval2 in Tval:
                for Hval2 in Hval:
                    coloc = colocalisation(Tval2[9],Tval2[10],Tval2[12],
                            Tval2[13],Hval2[9],Hval2[10],Hval2[12],Hval2[13])
                    if positive(Hval2[2] - Tval2[2]) > float(z_stacks_per_TAF):
                        pass
                    elif (coloc > float(bottom_overlap_ratio) and 
                          coloc < float(top_overlap_ratio) and
                          not Tval2[0:3] in TAF_TELO): #anti-ghosting included here
                        TAF_TELO.append(Tval2[0:3])
                        TELO_length.append(Tval2[15]) #15 is relative telomere length
                        TAF_H2AX.append(Hval2[0:3])
        n_TAF_TELO.append(len(TAF_TELO))
        TTAF[Tkey] = TAF_TELO[:]
        TELO_len[Tkey] = TELO_length[:]
        HTAF[Tkey] = TAF_H2AX[:]
        n_TAF[Tkey] = len(TAF_TELO[:])
        TAF_TELO.clear()
        TAF_H2AX.clear()
        TELO_length.clear()
        TAF_TELO.clear()
        num += 1
        if ((len(TTAF.get(Tkey)) < int(TAF_positive_threshold)) or 
            (len(TTAF.get(Tkey)) > int(upper_TAF_positive_threshold))):
            pass
        else:
            TAF_positive_nuclei.append("1")
    
    nuclei_multidict = {}
    for k,v in dict_nuclei.items():
        nuclei_multidict.setdefault(v,set().add(k))
        
    for k,v in nuclei_multidict.items():
        if not v is None and len(v) > 1:
            for k2,v2 in dict_nuclei.items:
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
    print(dict_nuclei)
    return [TTAF, TELO_len, n_TAF, TAF_percent_positive, dict_nuclei]

    
def sortby_treatment(dataset):
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
#have to set separate indices for DAPI, H2AX and TELO since they all have diff lengths

dataset_obj = sortby_treatment(dataset_DAPI)
dataset_indices = list(zip(treatment_index(dataset_DAPI),
                  treatment_index(dataset_H2AX),treatment_index(dataset_TELO)))
image_indices = list(zip(retrieve_index(x_DAPI),retrieve_index(x_H2AX),
                         retrieve_index(x_TELO)))
        
treatments_TTAF,treatments_pos,treatments_nTAF,treatments_Tlen = {},{},{},{}
treatments_nuclei = {}
images_TTAF, images_pos, images_Tlen, images_nTAF = {},{},{},{}
images_nuclei = {}
# =============================================================================
# for n, obj in enumerate(dataset_indices):
#     images_TTAF, images_pos, images_Tlen, images_nTAF = {},{},{},{}
#     images_nuclei = {}
#     for m, obj_2 in enumerate(image_indices):
#         if m > 0 and n <= len(dataset_indices):
#             if ((image_indices[m-1][0] >= (dataset_indices[n-1][0])) and
#                 (image_indices[m][0] <= (dataset_indices[n][0])+1)):
#                 Image_num = "Image_" + str(m)
#                 analysis = full_analysis(image_indices[m-1][0],image_indices[m][0],
#                                      image_indices[m-1][1],image_indices[m][1],
#                                      image_indices[m-1][2],image_indices[m][2])
#                 print(m)
#                 print(n)
#                 # output for analysis is a list of dicts
#                 images_TTAF[Image_num] = analysis[0]
#                 images_Tlen[Image_num] = analysis[1]
#                 images_nTAF[Image_num] = analysis[2]
#                 images_pos[Image_num] = analysis[3]
#                 images_nuclei[Image_num] = analysis[4]
#                 treatments_TTAF.update({dataset_obj[n-1] : images_TTAF})
#                 treatments_pos.update({dataset_obj[n-1] : images_pos})
#                 treatments_nTAF.update({dataset_obj[n-1] : images_nTAF})
#                 treatments_Tlen.update({dataset_obj[n-1] : images_Tlen})
#                 treatments_nuclei.update({dataset_obj[n-1] : images_nuclei})
#             else:
#                 pass
# 
# =============================================================================
analysis_1 = full_analysis(0,2,0,70,0,65)
analysis_2 = full_analysis(2,4,70,123,65,133)
analysis_3 = full_analysis(4,7,123,238,133,284)
num = 0
# =============================================================================
# while num < 3:
#     Image_num = "Image_" + str(num)
#     images_TTAF[Image_num] = "analysis_" + str(num) + "[0]") #this doesnt work > excel output incorrect
#     images_Tlen[Image_num] = "analysis_" + str(num) + "[1]") #check hack method below
#     images_nTAF[Image_num] = "analysis_" + str(num) + "[2]"
#     images_pos[Image_num] = "analysis_" + str(num) + "[3]"
#     images_nuclei[Image_num] = "analysis_" + str(num) + "[4]"
#     treatments_TTAF.update({dataset_obj[num-1] : images_TTAF})
#     treatments_pos.update({dataset_obj[num-1] : images_pos})
#     treatments_nTAF.update({dataset_obj[num-1] : images_nTAF})
#     treatments_Tlen.update({dataset_obj[num-1] : images_Tlen})
#     treatments_nuclei.update({dataset_obj[num-1] : images_nuclei})
#     num += 1
#     print("X")
# =============================================================================
    
# =============================================================================
# images_TTAF["Image_1"] = analysis_1[0] #this doesnt work > excel output incorrect
# images_Tlen["Image_1"] = analysis_1[1]
# images_nTAF["Image_1"] = analysis_1[2]
# images_pos["Image_1"] = analysis_1[3]
# images_nuclei["Image_1"] = analysis_1[4]
# =============================================================================
treatments_TTAF.update({dataset_obj[0] : analysis_1[0]})
treatments_pos.update({dataset_obj[0] : analysis_1[1]})
treatments_nTAF.update({dataset_obj[0] : analysis_1[2]})
treatments_Tlen.update({dataset_obj[0] : analysis_1[3]})
treatments_nuclei.update({dataset_obj[0] : analysis_1[4]})

# =============================================================================
# images_TTAF["Image_2"] = analysis_2[0] #this doesnt work > excel output incorrect
# images_Tlen["Image_2"] = analysis_2[1]
# images_nTAF["Image_2"] = analysis_2[2]
# images_pos["Image_2"] = analysis_2[3]
# images_nuclei["Image_2"] = analysis_2[4]
# =============================================================================
treatments_TTAF.update({dataset_obj[1] : analysis_2[0]})
treatments_pos.update({dataset_obj[1] : analysis_2[1]})
treatments_nTAF.update({dataset_obj[1] : analysis_2[2]})
treatments_Tlen.update({dataset_obj[1] : analysis_2[3]})
treatments_nuclei.update({dataset_obj[1] : analysis_2[4]})

# =============================================================================
# images_TTAF["Image_3"] = analysis_2[0] #this doesnt work > excel output incorrect
# images_Tlen["Image_3"] = analysis_2[1]
# images_nTAF["Image_3"] = analysis_2[2]
# images_pos["Image_3"] = analysis_2[3]
# images_nuclei["Image_3"] = analysis_2[4]
# =============================================================================
treatments_TTAF.update({dataset_obj[2] : analysis_3[0]})
treatments_pos.update({dataset_obj[2] : analysis_3[1]})
treatments_nTAF.update({dataset_obj[2] : analysis_3[2]})
treatments_Tlen.update({dataset_obj[2] : analysis_3[3]})
treatments_nuclei.update({dataset_obj[2] : analysis_3[4]})

mean_pos_treatment = {}
temp_percent = []
# =============================================================================
# for n,(treatments,images) in enumerate(treatments_pos.items()):
#     for k,v in images.items():
#         if treatments == dataset_obj[n]:
#             temp_percent.append(v[0])
#     if len(temp_percent) == len(images) and len(temp_percent) > 1:
#         temp_percent = [st.mean(temp_percent),st.stdev(temp_percent),
#             (st.stdev(temp_percent)/sqrt(len(images))),len(images)]
#         mean_pos_treatment.update({dataset_obj[n] : temp_percent})
#         temp_percent = []
#     else:
#         temp_percent = [0,0,0,1]
#         mean_pos_treatment.update({dataset_obj[n] : temp_percent})
#         temp_percent = []
# =============================================================================
        
# percentage nuclei positive for senescence by treatment
dfmeanpos = pd.DataFrame.from_dict(mean_pos_treatment,orient='index',
        columns=['Percent TAF positive','SD','SE', 'n'])
# percentage of nuclei positive for senescence by images
# =============================================================================
# dfpos = pd.DataFrame.from_dict({(i,j): treatments_pos[i][j]
#         for i in treatments_pos.keys()
#         for j in treatments_pos[i].keys()},
#         orient='index', columns=['Percent TAF positive',
#                                  'TAF positive count','Total nuclei count'])
# =============================================================================
# coordinate of each TAF positive telomere
dfTTAF = pd.DataFrame.from_dict({(i,j): treatments_TTAF[i][j]
                                    for i in treatments_TTAF.keys()
                                    for j in treatments_TTAF[i].keys()},
                                    orient='index')
# relative length of each TAF positive telomere
# =============================================================================
# dfTlen = pd.DataFrame.from_dict({(i,j): treatments_Tlen[i][j]
#                                     for i in treatments_Tlen.keys()
#                                     for j in treatments_Tlen[i].keys()},
#                                     orient='index')
# =============================================================================
# summarised number of TAF per nucleus
dfnTAF = pd.DataFrame.from_dict({(i,j): treatments_nTAF[i][j]
                                    for i in treatments_nTAF.keys()
                                    for j in treatments_nTAF[i].keys()},
                                    orient='index')
#df_allpos = pd.concat([dfmeanpos,dfpos], axis=1, sort=False)
dfnuclei = pd.DataFrame.from_dict({(i,j): treatments_nuclei[i][j]
                                    for i in treatments_nuclei.keys()
                                    for j in treatments_TTAF[i].keys()},
                                    orient='index')

#dfTTAF = pd.DataFrame.from_dict(treatments_TTAF, orient='index')
dfTlen = pd.DataFrame.from_dict(treatments_Tlen, orient='index')
#dfnTAF = pd.DataFrame.from_dict(treatments_nTAF, orient='index')
dfnuclei = pd.DataFrame.from_dict(treatments_nuclei, orient='index')

def dfs_tabs(df_list, sheet_list, file_name):
    writer = pd.ExcelWriter(file_name,engine='xlsxwriter')   
    for dataframe, sheet in zip(df_list, sheet_list):
        dataframe.to_excel(writer, sheet_name=sheet, startrow=0 , startcol=0)   
    writer.save()

#df = [df_allpos,dfTTAF,dfTlen,dfnTAF,dfnuclei]
df = [dfTTAF,dfTlen,dfnTAF,dfnuclei]
#sheets = ["Percent positive","TAF coordinates","Relative Telo length","n TAF per nucleus","Nuclear coordinates"]
sheets = ["TAF coordinates","Relative Telo length","n TAF per nucleus","Nuclear coordinates"]

dfs_tabs(df,sheets,output_file)

end = time.time()
print("Runtime = %s" % (end - start))