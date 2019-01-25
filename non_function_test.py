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

root = 'C:\\Users\\Ishaq\\Documents\\All spreadsheets\\Ed\\' #remember to use double backslash or single forward slash (is a requirement for python convention)
input_file_H2AX = root + 'All_H2AX.xlsx'
input_file_TELO = root + 'All_TELO.xlsx'
input_file_DAPI = root + 'All_DAPI.xlsx'
output_file = root + str(date.today()) + '_project_all_after_script.xlsx'

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

for point in all_DAPI:
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
all_DAPI = list(zip(x_DAPI, y_DAPI, 
                    z_DAPI, x_DAPI_start,
                    x_DAPI_end, y_DAPI_start,y_DAPI_end,
                    z_dim_DAPI)) #one less list level to iterate
all_DAPI.remove(all_DAPI[0])

#converts all pixels into microns, point comparisons in microns, point vs DAPI comparisons in vectors
#px = pixel size
for vector in all_H2AX:
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
all_H2AX = list(zip(x_H2AX,y_H2AX,
                    z_H2AX,width_H2AX, 
                    height_H2AX,depth_H2AX, 
                    xmicron_H2AX, ymicron_H2AX, zmicron_H2AX,
                    xmicron_H2AX_start,ymicron_H2AX_start,zmicron_H2AX_end,
                    xmicron_H2AX_end,ymicron_H2AX_end,zmicron_H2AX_end))
all_H2AX.remove(all_H2AX[0])

#converts all pixels into microns, point comparisons in microns, point vs DAPI comparisons in vectors
for vector in all_TELO:
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
all_TELO = list(zip(x_TELO, y_TELO, 
                    z_TELO, avint_TELO))

maxsize_TELO = max(all_TELO[17])
def telo_rellen(avint_TELO): #relative telo len in spreadsheet
    return avint_TELO/maxsize_TELO

for vector in all_TELO:
    if vector[0] == all_TELO[0][0]:
        rellen_TELO.append("Relative Telo Length")
    else:
        rellen_TELO.append(telo_rellen(vector[3]))
all_TELO = list(zip(x_TELO,y_TELO,
                    z_TELO,width_TELO,
                    height_TELO,depth_TELO,
                    xmicron_TELO,ymicron_TELO,zmicron_TELO,xmicron_TELO_start,
                    ymicron_TELO_start,zmicron_TELO_start,
                    xmicron_TELO_end,ymicron_TELO_end,zmicron_TELO_end,
                    rellen_TELO,maxint_TELO,
                    avint_TELO))
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
for i in range(len(all_DAPI)-1):
    dict_nuclei["Nucleus no. " + str(i)] = all_DAPI[i]
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
                      coloc < float(top_overlap_ratio)):
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

TAF_multidict = {}
for k,v in TTAF.items():
    TAF_multidict.setdefault(v,set().add(k))

TAF_percent_positive.append(len(TAF_positive_nuclei)/len(TTAF)*100) #percentage, count positive, count total
TAF_percent_positive.append(len(TAF_positive_nuclei))
TAF_percent_positive.append(len(TTAF))
#return [TTAF, TELO_len, n_TAF, TAF_percent_positive]
