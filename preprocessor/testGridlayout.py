from util.umep_suewsss_export_component import write_GridLayout_file, create_GridLayout_dict
import matplotlib.pylab as plt
import numpy as np


#TODO make into new function later



test = 4
outputDir = "c:/Users/xlinfr/Documents/PythonScripts"
_name = 'testar'
ssVect = np.loadtxt(r'c:\temp\TEstSS\ss_IMPGrid_SS_5.txt',skiprows=1) # input from gui
imp_5 = np.loadtxt(r'c:\temp\TEstSS\ss_IMPGrid_isotropic.txt',skiprows=1) # input from gui
heightMethod = 2 # input from gui
vertheights = [0, 15, 40]
nlayer = 3 # input from gui
skew = 2 #TODO input from gui. linear or shewed shewed = 1 or 2

#zHmax = imp_5[5,4]

ssDict = create_GridLayout_dict()

if heightMethod == 1: # static levels (taken from interface). Last value > max height
    ssDict['height'] = vertheights.append(ssVect[:,0].max())
    ssDict['nlayer'] = len(ssDict['height']) - 1
elif heightMethod == 2: # always nlayers layer based on percentiles
    ssDict['nlayer'] = nlayer
elif heightMethod == 3: # vary number of layers based on height variation. Lowest no of nlayers always 3
    nlayer = 3
    if ssVect[:,0].max() > 40: nlayer = 4
    if ssVect[:,0].max() > 60: nlayer = 5
    if ssVect[:,0].max() > 80: nlayer = 6
    if ssVect[:,0].max() > 120: nlayer = 7
    ssDict['nlayer'] = nlayer

intervals = np.ceil(ssVect[:,0].max() / nlayer)
ssDict['height'] = []
ssDict['height'].append(.0)
for i in range(1, nlayer):
    ssDict['height'].append(float(round((intervals * i) / skew)))
ssDict['height'].append(float(ssVect[:,0].max()))

ssDict['building_frac'] = []
ssDict['veg_frac'] = []
ssDict['building_scale'] = []
ssDict['veg_sale'] = []

ssDict['building_frac'].append(ssVect[0,1]) # first is plan area index of buildings 
ssDict['veg_frac'].append(ssVect[0,3]) # first is plan area index of trees 

index = int(0)
for i in range(1,len(ssDict['height'])):
    #np.where(ssVect[:,0] == ssDict['height'][ind]) #Can be used if ssVect is not in meter intervals
    print(i)
    index += 1
    ssDict['building_frac'].append(np.mean(ssVect[int(ssDict['height'][index-1]):int(ssDict['height'][index]), 1])) # intergrated pai_build mean in ith vertical layer
    ssDict['veg_frac'].append(np.mean(ssVect[int(ssDict['height'][index -1]):int(ssDict['height'][index]), 3])) # intergrated pai_veg mean in ith vertical layer
    ssDict['building_scale'].append(np.mean(ssVect[int(ssDict['height'][index -1]):int(ssDict['height'][index]), 2])) # intergrated bscale mean in ith vertical layer
    ssDict['veg_scale'].append(np.mean(ssVect[int(ssDict['height'][index -1]):int(ssDict['height'][index]), 4])) # intergrated vscale mean in ith vertical layer

#TODO here we nned to add other parameters based on typology

write_GridLayout_file(ssDict, outputDir + '/', _name)