import pandas as pd
import numpy as np
import sys
from numba import jit
import time

print('Python Version ' + sys.version)
print('Pandas Version ' + pd.__version__)
print('Np Version ' + np.__version__)

root = 'H:\\Ed TAF\\20181120 test 2\\' #remember to use double backslash or single forward slash (is a requirement for python convention)
input_file_H2AX = root + 'All_H2AX.xlsx'
input_file_TELO = root + 'All_TELO.xlsx'
input_file_DAPI = root + 'All_DAPI.xlsx'
output_file = root + 'All_TTAF_after_script.csv'
output_file_2 = root + 'All_HTAF_after_script.csv'
output_file_3 = root + 'All_Telo_len_after_script.csv'

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
x_H2AX = np.array(df_H2AX[3])
y_H2AX = np.array(df_H2AX[4])
z_H2AX = np.array(df_H2AX[5])
width_H2AX = np.array(df_H2AX[6])
height_H2AX = np.array(df_H2AX[7])
depth_H2AX = np.array(df_H2AX[8])
x_TELO = np.array(df_TELO[3])
y_TELO = np.array(df_TELO[4])
z_TELO = np.array(df_TELO[5])
width_TELO = np.array(df_TELO[6])  
height_TELO = np.array(df_TELO[7]) 
depth_TELO = np.array(df_TELO[8])
maxint_TELO = np.array(df_TELO[21])
avint_TELO = np.array(df_TELO[20])
x_DAPI = np.array(df_DAPI[3])
y_DAPI = np.array(df_DAPI[4])
z_DAPI = np.array(df_DAPI[5])
width_DAPI = np.array(df_DAPI[6])
height_DAPI = np.array(df_DAPI[7])
depth_DAPI = np.array(df_DAPI[8])
x_dim_DAPI = []
y_dim_DAPI = []
z_dim_DAPI = []
ROI_end_DAPI = []
nuclear_count = []
x_dim_H2AX = []
y_dim_H2AX = []
z_dim_H2AX = []
xmicron_H2AX = []
ymicron_H2AX = []
zmicron_H2AX = []
area_H2AX = []
vol_H2AX = []
sxmicron_H2AX = [] #size
symicron_H2AX = []
szmicron_H2AX = []
filt_H2AX = []
xmicron_H2AX_end = []
ymicron_H2AX_end = []
zmicron_H2AX_end = []
totaled_H2AX_count = []
H2AX_count = []
H2AX_volume = []
ROI_end_H2AX = []
IMG_no_H2AX = []
x_dim_TELO = []
y_dim_TELO = []
z_dim_TELO = []
xmicron_TELO = []
ymicron_TELO = []
zmicron_TELO = []
area_TELO = []
vol_TELO = []
sxmicron_TELO = []
symicron_TELO = []
szmicron_TELO = []
filt_TELO = []
rellen_TELO = []
xmicron_TELO_end = []
ymicron_TELO_end = []
zmicron_TELO_end = []
totaled_TELO_count = []
TELO_count = []
TELO_volume = []
ROI_end_TELO = []
IMG_no_TELO = []
TAF_count = []
totaled_TAF_count = []
values = []

all_H2AX = list(zip(x_H2AX, y_H2AX, z_H2AX, width_H2AX, height_H2AX, depth_H2AX))
all_TELO = list(zip(x_TELO, y_TELO, z_TELO, width_TELO, height_TELO, depth_TELO, avint_TELO))
all_DAPI = list(zip(x_DAPI, y_DAPI, z_DAPI, width_DAPI, height_DAPI, depth_DAPI))

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
    
def nuclear_filter(p,p0,p1,DAPI1,DAPI2,DAPI3,DAPI4):
    if (      floatify(point[0]) > floatify(DAPI1) #filters for within nuclear regions only
                and floatify(point[0]) < floatify(DAPI2) #dapi_vector0 < point0 < dapi_vector3
                and floatify(point[1]) > floatify(DAPI3) 
                and floatify(point[1]) < floatify(DAPI4)):
        return True
    else:
        return False

def convert_micron(pixel_size,V1):
    return floatify(V1*pixel_size)

def convert_size_micron(pixel_size,V1,V2):
    return floatify(V1*pixel_size)+floatify(V2*pixel_size)

def full_analysis(index_1,index_2):
    for point in all_DAPI[index_1,index_2]:
        if point[0] == all_DAPI[0][0]: #numbers the images
            x_dim_DAPI.append('Position X')
            y_dim_DAPI.append('Position Y')
            z_dim_DAPI.append('Position Z')
            nuclear_count.append(len(x_dim_DAPI)-1) #first iter will show 0, last iter outside loop, nuclei = total so far. Do subtractions
        else:
            x_dim_DAPI.append(point[0] + point[3]) #x_dim_calc
            y_dim_DAPI.append(point[1] + point[4]) #y_dim_calc
            z_dim_DAPI.append(float(1)) #z_dim_calc
    nuclear_count.append(len(x_dim_DAPI))
    #nuclear_count.remove(0)
    start_end_vectors_DAPI_merged = list(zip(x_DAPI, y_DAPI, z_DAPI, x_dim_DAPI, y_dim_DAPI, z_dim_DAPI)) #one less list level to iterate
    
    #converts all pixels into microns, point comparisons in microns, point vs DAPI comparisons in vectors
    #px = pixel size
    for vector in all_H2AX[index_1,index_2]:
        if vector[0] == all_H2AX[0][0]:
            xmicron_H2AX.append("X") #converts x-val into microns (max is 2080x0.16)
            ymicron_H2AX.append("Y")
            zmicron_H2AX.append("Z")
            area_H2AX.append("Area")
            vol_H2AX.append("Vol.")
            sxmicron_H2AX.append("sX") #converts box width values into microns
            symicron_H2AX.append("sY")
            szmicron_H2AX.append("sZ")
            xmicron_H2AX_end.append("X end") #adds box width in microns to x-val - calculates len of box side
            ymicron_H2AX_end.append("Y end")
            zmicron_H2AX_end.append("Z end")
        else:
            xmicron_H2AX.append(convert_micron(px,vector[0]))
            ymicron_H2AX.append(convert_micron(px,vector[1]))
            zmicron_H2AX.append(convert_micron(px,vector[2]))
            area_H2AX.append(convert_size_micron(px,vector[0],vector[1]))
            vol_H2AX.append((floatify(vector[0])*px)*(floatify(vector[1])*px)*(floatify(vector[2])*z_size))
            sxmicron_H2AX.append(convert_micron(px,vector[3]))
            symicron_H2AX.append(convert_micron(px,vector[4]))
            szmicron_H2AX.append(convert_micron(px,vector[5]))
            xmicron_H2AX_end.append(convert_size_micron(px,vector[0],vector[3]))
            ymicron_H2AX_end.append(convert_size_micron(px,vector[1],vector[4]))
            zmicron_H2AX_end.append(convert_size_micron(px,vector[2],vector[5]))
    all_H2AX = list(zip(x_H2AX, y_H2AX, z_H2AX, width_H2AX, height_H2AX, 
                        depth_H2AX, xmicron_H2AX, ymicron_H2AX, zmicron_H2AX, 
                        sxmicron_H2AX, symicron_H2AX, szmicron_H2AX, 
                        xmicron_H2AX_end, ymicron_H2AX_end, zmicron_H2AX_end))
    
    #converts all pixels into microns, point comparisons in microns, point vs DAPI comparisons in vectors
    for vector in all_TELO[index_1,index_2]:
        if vector[0] == all_TELO[0][0]:
            xmicron_TELO.append("X") #converts x-val into microns (max is 2080x0.16)
            ymicron_TELO.append("Y")
            zmicron_TELO.append("Z")
            area_TELO.append("Area")
            vol_TELO.append("Vol.")
            sxmicron_TELO.append("sX") #converts box width values into microns
            symicron_TELO.append("sY")
            szmicron_TELO.append("sZ")
            xmicron_TELO_end.append("X end") #adds box width in microns to x-val - calculates len of box side
            ymicron_TELO_end.append("Y end")
            zmicron_TELO_end.append("Z end")
        else:
            xmicron_TELO.append(convert_micron(px,vector[0]))
            ymicron_TELO.append(convert_micron(px,vector[1]))
            zmicron_TELO.append(convert_micron(px,vector[2]))
            area_TELO.append(convert_size_micron(px,vector[0],vector[1]))
            vol_TELO.append((floatify(vector[0])*px)*(floatify(vector[1])*px)*(floatify(vector[2])*z_size))
            sxmicron_TELO.append(convert_micron(px,vector[3]))
            symicron_TELO.append(convert_micron(px,vector[4]))
            szmicron_TELO.append(floatify(vector[5])*z_size)
            xmicron_TELO_end.append(convert_size_micron(px,vector[0],vector[3]))
            ymicron_TELO_end.append(convert_size_micron(px,vector[1],vector[4]))
            zmicron_TELO_end.append(convert_size_micron(px,vector[2],vector[5]))
    all_TELO = list(zip(x_TELO, y_TELO, z_TELO, avint_TELO))
    
    maxsize_TELO = max(all_TELO[17])
    def telo_rellen(avint_TELO): #relative telo len in spreadsheet
        return avint_TELO/maxsize_TELO
    
    for vector in all_TELO:
        if vector[0] == all_TELO[0][0]:
            rellen_TELO.append("Relative Telo Length")
        else:
            rellen_TELO.append(telo_rellen(vector[3]))
    all_TELO = list(zip(x_TELO, y_TELO, z_TELO, width_TELO, height_TELO, 
                        depth_TELO, xmicron_TELO, ymicron_TELO, zmicron_TELO, 
                        sxmicron_TELO, symicron_TELO, szmicron_TELO, 
                        xmicron_TELO_end, ymicron_TELO_end, zmicron_TELO_end,
                        rellen_TELO, maxint_TELO, avint_TELO))

    #converts all pixels into microns, point comparisons in microns, point vs DAPI comparisons in vectors
    num = 0
    for DAPI_vectors in start_end_vectors_DAPI_merged:
        for point in all_H2AX:
            if point[0] == all_H2AX[0][0]:
                totaled_H2AX_count.append(len(filt_H2AX))
            elif DAPI_vectors[0] == start_end_vectors_DAPI_merged[0][0]:
                pass #dapi_vector0 < point0 < dapi_vector3
            elif nuclear_filter(point,point[0],point[1],DAPI_vectors[0],DAPI_vectors[3],DAPI_vectors[1],DAPI_vectors[4]) == True:
                if (not H2AX_size_threshold == 0) and (point[9]       #filters by H2AX foci size
                     or point[10] or point[11] >= H2AX_size_threshold):
                    pass
                    filt_H2AX.append(point)
                elif H2AX_size_threshold == 0:
                    filt_H2AX.append(point)
            else:
                pass
    #totaled_H2AX_count.remove(0)
        
    for DAPI_vectors in start_end_vectors_DAPI_merged:
        for point in all_TELO:
            if point[0] == all_TELO[0][0]:
                totaled_TELO_count.append(len(filt_TELO))
            elif DAPI_vectors[0] == start_end_vectors_DAPI_merged[0][0]:
                pass 
            elif nuclear_filter(point,point[0],point[1],DAPI_vectors[0],DAPI_vectors[3],DAPI_vectors[1],DAPI_vectors[4]) == True:
                if (not TELO_size_threshold == 0) and (point[9]       #filters by H2AX foci size
                     or point[10] or point[11] >= TELO_size_threshold):
                    pass
                    filt_H2AX.append(point)
                elif TELO_size_threshold == 0:
                    filt_TELO.append(point)
            else:
                pass
    #totaled_TELO_count.remove(0)

    dict_nuclei_H2AX = {}
    dict_nuclei_TELO = {}
    dict_H2AX_count = {}
    dict_TELO_count = {}
    for i in range(len(start_end_vectors_DAPI_merged)-1):
        dict_nuclei_H2AX["Nucleus no. " + str(i)] = (filt_H2AX[totaled_H2AX_count[max(i-1,0)]:totaled_H2AX_count[i]])
        dict_nuclei_TELO["Nucleus no. " + str(i)] = (filt_TELO[totaled_TELO_count[max(i-1,0)]:totaled_TELO_count[i]])
        dict_H2AX_count["Nucleus no. " + str(i)] = len((filt_H2AX[totaled_H2AX_count[max(i-1,0)]:totaled_H2AX_count[i]]))
        dict_TELO_count["Nucleus no. " + str(i)] = len((filt_TELO[totaled_TELO_count[max(i-1,0)]:totaled_TELO_count[i]]))
    #nuclei 0 and 1 are showing 0 count, but 2 is showing 10
        
    @jit
    def left(x1,x2):
        return max(x1,x2)
    
    @jit
    def right(w1,w2):
        return min(w1,w2)
    
    @jit
    def bottom(y1,y2):
        return max(y1,y2)
    
    @jit
    def top(h1,h2):
        return min(h1,h2)
    
    @jit
    def colocalisation(x1,y1,w1,h1,x2,y2,w2,h2):
        left = max(x1,x2)
        right = min(w1,w2)
        bottom = max(y1,y2)
        top = min(h1,h2)
        Aoverlap = (bottom - top)*(right - left)
        ratio_to_TELO = Aoverlap/((w1-x1)*(h1-y1))
        #ratio_to_H2AX = Aoverlap/((w2-x2)*(h2-y2))
        if ratio_to_TELO > 0 and ratio_to_TELO < 5:
            return ratio_to_TELO
        else:
            return False
    
    TTAF = {}
    TTAF_len = {}
    HTAF = {}
    n_TAF = {}
    TAF_TELO = []
    TAF_H2AX = []
    TELO_length = []
    n_TAF_TELO = []
    TAF_positive_nuclei = []
    num = 0
    #colocalisation of H2AX and TELO
    
    z_stacks_per_TAF = float(TAF_size_threshold) / float(z_size)
    for (Tkey, Tval), (Hkey, Hval) in zip(dict_nuclei_TELO.items(), dict_nuclei_H2AX.items()):
        if Hkey == Tkey:
            for Tval2 in Tval:
                for Hval2 in Hval:
                    coloc = colocalisation(Tval2[6],Tval2[7],Tval2[12],Tval2[13],Hval2[6],Hval2[7],Hval2[12],Hval2[13])
                    if positive(Hval2[2] - Tval2[2]) > float(z_stacks_per_TAF):
                        pass
                    elif coloc > float(bottom_overlap_ratio) and coloc < float(top_overlap_ratio):
                        TAF_TELO.append(Tval2[0:3])
                        TELO_length.append(Tval2[15]) #15 is relative telomere length
                        TAF_H2AX.append(Hval2[0:3])
        n_TAF_TELO.append(len(TAF_TELO))
        TTAF[Tkey] = TAF_TELO[:]
        TTAF_len[Tkey] = TELO_length[:]
        HTAF[Tkey] = TAF_H2AX[:]
        n_TAF[Tkey] = len(TAF_TELO[:])
        TAF_TELO.clear()
        TAF_H2AX.clear()
        TELO_length.clear()
        #print(TAF_TELO)
        TAF_TELO.clear()
        num += 1
        #print("Nucleus " + str(num) + " of " + str(len(start_end_vectors_DAPI_merged)) + " done.")
        if (len(TTAF.get(Tkey)) < int(TAF_positive_threshold)) or (len(TTAF.get(Tkey)) > int(upper_TAF_positive_threshold)):
            pass
        else:
            TAF_positive_nuclei.append("1")
        #print("Runtime = %s" % (end - start))
    
    TAF_percent_positive = []
    TAF_percent_positive.append(len(TAF_positive_nuclei)/len(TTAF)*100)
    TAF_percent_positive.append(len(TAF_positive_nuclei))
    TAF_percent_positive.append(len(TTAF))

    dffoci = pd.DataFrame.from_dict(dict_H2AX_count, orient='index', columns=['H2A.X count'])
    dffoci_2 = pd.DataFrame.from_dict(dict_TELO_count, orient='index', columns=['Telomere count'])
    dfTTAF = pd.DataFrame.from_dict(TTAF, orient='index')
    dfTTAF_len = pd.DataFrame.from_dict(TTAF_len, orient='index')
    dfHTAF = pd.DataFrame.from_dict(HTAF, orient='index')
    dfnTAF = pd.DataFrame.from_dict(n_TAF, orient='index', columns=['TAF count'])
    dfTAFpercent = pd.DataFrame(data = TAF_percent_positive, index=['Percent TAF positive','TAF positive count','Total nuclei count'])
    allTdf = pd.concat([dfTAFpercent,dffoci,dffoci_2,dfnTAF,dfTTAF_len,dfTTAF], axis=1, sort=False)
    allHdf = pd.concat([dfTAFpercent,dffoci,dffoci_2,dfnTAF,dfHTAF], axis=1, sort=False)

    return allTdf, allHdf, dfTTAF_len, dfTAFpercent

all_images = {}

def retrieve_index(df,index_list):
    index_list = []
    num = 1
    for index_2, obj in enumerate(df):
        if index_2 == 0:
            index_list.append(0)
        elif obj == x_H2AX[0]:
            index_list.append(index_2)
            image_num = "image_" + str(num)
            index_1 = index_list[num-1] if num > 1 else 0
            all_images[image_num] = full_analysis(index_1,index_2-1)
            num += 1
        elif index_2 + 1 == len(x_H2AX):
            image_num = "image_" + str(num)
            index_1 = index_list[num-1]
            all_images[image_num] = full_analysis(index_1,index_2-1)
            break
        else:
            pass
#have to set separate indices for DAPI, H2AX and TELO since they all have diff lengths
    
DAPI_index = retrieve_index(x_DAPI)
H2AX_index = retrieve_index(x_H2AX)
TELO_index = retrieve_index(x_TELO)

        
        
# =============================================================================
# allTdf.to_csv(output_file, index=True)
# allHdf.to_csv(output_file_2, index=True)
# dfTTAF_len.to_csv(output_file_3, index=True)
# 
# =============================================================================
end = time.time()
print("Runtime = %s" % (end - start))