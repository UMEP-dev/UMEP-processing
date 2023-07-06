
'''
This module is used to transform UMEP pre-processing data to fit SuewsSS

input: 
ss_object: dictonary with all relevant inputs
refdir = outputfolder
fname = filename       
'''

import os

def create_GridLayout_dict():

    ssDict = {}

    # dim
    ssDict['nlayer'] = 3 # number of vertical layers

    # geom
    ssDict['height'] = [0., 11., 15., 22.] #TODO  height of top of each layer (start with 0 i.e. one more than nlayers)
    ssDict['building_frac'] = [0.43, 0.38, .2] #TODO  fraction of building coverage of each layer; the first one is plan area index of buildings
    ssDict['veg_frac'] = [0.01, 0.02, .01] #TODO  fraction of vegetation coverage of each layer
    ssDict['building_scale'] =[50., 50., 50] #TODO building scale of each layer [m]
    ssDict['veg_scale'] = [10., 10., 10] #TODO  vegetation scale of each layer [m]

    # roof
    ssDict['sfr_roof'] = [.3, .3, .4] #TODO how? fraction of roofs of each layer (sum should be 1)
    ssDict['tin_roof'] = [5, 5, 6] #TODO? how?  initial temperatures of roofs [degC]
    ssDict['alb_roof'] = [.5, .5, .2]  #TODO albedo of roofs
    ssDict['emis_roof'] = [.95, .95, .95] #TODO emissivity of roofs
    ssDict['state_roof'] = [.0, .0, .0]  # initial surface water depth state of roofs [mm]
    ssDict['statelimit_roof'] = [5, 5, 5] # surface water depth state limit of roofs [mm]
    ssDict['wetthresh_roof'] = [5, 5, 5] # surface water depth threshold of roofs (used in latent heat flux calculation) [mm]
    ssDict['soilstore_roof'] = [20, 20, 20] # soil water store of roofs [mm]
    ssDict['soilstorecap_roof'] = [120, 120, 120] # soil water store capacity of roofs [mm]
    ssDict['roof_albedo_dir_mult_fact'] = [1.,1.,1.] # initial surface water depth state of roofs [mm]

    # The following parameters are used to calculate the heat flux from the roof
    # first roof facet
    ssDict['dz_roof1'] = [.2, .1, .1, .01, .01] #TODO thickness of each layer (strictly five lyaers in total) [m]
    ssDict['k_roof1'] = [1.2, 1.2, 1.2, 1.2, 1.2] #TODO thermal conductivity of each layer [W/m/K]
    ssDict['cp_roof1'] = [2e6, 2e6, 2e6, 2e6, 2e6] #TODO specific heat capacity of each layer [J/kg/K]

    ssDict['dz_roof2'] = [.2, .1, .1, .01, .01] #TODO
    ssDict['k_roof2'] = [2.2, 1.2, 1.2, 1.2, 1.2] #TODO
    ssDict['cp_roof2'] = [2e6, 2e6, 2e6, 2e6, 2e6] #TODO

    ssDict['dz_roof3'] = [.2, .1, .1, .01, .01] #TODO
    ssDict['k_roof3'] = [2.2, 1.2, 1.2, 1.2, 1.2] #TODO
    ssDict['cp_roof3'] = [2e6, 2e6, 2e6, 2e6, 2e6] #TODO

    # wall # similarly to roof parameters but for walls
    ssDict['sfr_wall'] = [.3, .3, .4] #TODO # (sum should be 1)
    ssDict['tin_wall'] = [5, 5, 5]
    ssDict['alb_wall'] = [.5, .5, .5]#TODO
    ssDict['emis_wall'] = [.95, .95, .95]#TODO
    ssDict['state_wall'] = [.0, .0, .0]
    ssDict['statelimit_wall'] = [5, 5, 5]
    ssDict['wetthresh_wall'] = [5, 5, 5]
    ssDict['soilstore_wall'] = [20, 20, 20]
    ssDict['soilstorecap_wall'] = [120, 120, 120]
    ssDict['wall_specular_frac'] = [0.,0.,0.]

    ssDict['dz_wall1'] = [.2,  .1,  .1,  .01, .01]#TODO
    ssDict['k_wall1'] = [1.2, 1.2, 1.2, 1.2, 1.2]#TODO
    ssDict['cp_wall1'] = [3e6, 2e6, 2e6, 2e6, 2e6]#TODO

    ssDict['dz_wall2'] = [.2,  .1,  .1,  .01, .01]#TODO
    ssDict['k_wall2'] = [2.2, 1.2, 1.2, 1.2, 1.2]#TODO
    ssDict['cp_wall2'] = [2e6, 3e6, 2e6, 2e6, 2e6]#TODO

    ssDict['dz_wall3'] = [.2,  .1,  .1,  .01, .01]#TODO
    ssDict['k_wall3'] = [2.2, 1.2, 1.2, 1.2, 1.2]#TODO
    ssDict['cp_wall3'] = [2e6, 3e6, 2e6, 2e6, 2e6]#TODO

    # surf # for generic SUEWS surfaces, used in storage heat flux calculations
    ssDict['tin_surf'] = [2, 2, 2, 2, 2, 2, 2] # intitial temperature

    ssDict['dz_surf_paved'] = [.2,    .15,   .01,   .01,   .01]
    ssDict['k_surf_paved'] = [1.1,   1.1,   1.1,   1.1,   1.1]
    ssDict['cp_surf_paved'] = [2.2e6, 2.2e6, 2.2e6, 2.2e6, 2.6e6]

    ssDict['dz_surf_buildings'] = [.2,    .1,    .1,    .5,    1.6]
    ssDict['k_surf_buildings'] = [1.2,   1.1,   1.1,   1.5,   1.6]
    ssDict['cp_surf_buildings'] = [1.2e6, 1.1e6, 1.1e6, 1.5e6, 1.6e6]

    ssDict['dz_surf_evergreen'] = [.2,    .15,   .01,   .01,   .01]
    ssDict['k_surf_evergreen'] = [1.1,   1.1,   1.1,   1.1,   1.1]
    ssDict['cp_surf_evergreen'] = [3.2e6, 1.1e6, 1.1e6, 1.5e6, 1.6e6]

    ssDict['dz_surf_decid'] = [.2,    .1,    .1,    .1,    2.2]
    ssDict['k_surf_decid'] = [1.2,   1.1,   1.1,   1.5,   1.6]
    ssDict['cp_surf_decid'] = [3.2e6, 1.1e6, 1.1e6, 1.5e6, 1.6e6]

    ssDict['dz_surf_grass'] = [.2,    .05,   .1,    .1,    2.2]
    ssDict['k_surf_grass'] = [1.2,   1.1,   1.1,   1.5,   1.6]
    ssDict['cp_surf_grass'] = [1.6e6, 1.1e6, 1.1e6, 1.5e6, 1.6e6]

    ssDict['dz_surf_baresoil'] = [.2,    .05,   .1,    .1,    2.2]
    ssDict['k_surf_baresoil'] = [1.2,   1.1,   1.1,   1.5,   1.6]
    ssDict['cp_surf_baresoil'] = [1.9e6, 1.1e6, 1.1e6, 1.5e6, 1.6e6]

    ssDict['dz_surf_water'] = [.2,    .05,   .1,    .1,    2.2]
    ssDict['k_surf_water'] = [1.2,   1.1,   1.1,   1.5,   1.6]
    ssDict['cp_surf_water'] = [1.9e6, 1.1e6, 1.1e6, 1.5e6, 1.6e6]

    return ssDict


def write_GridLayout_file(ss_object, refdir, fname):
    ss_file_path = os.path.join(refdir,fname+".nml")
    f = open(ss_file_path, "w")

    f.write("&dim\n")
    f.write("nlayer = {}\n".format(ss_object['nlayer']))
    f.write("/\n")
    f.write("\n")

    f.write("&geom\n")
    f.write("height = {}\n".format(str(ss_object['height'])[1:-1]))
    f.write("building_frac = {}\n".format(str(ss_object['building_frac'])[1:-1]))
    f.write("veg_frac = {}\n".format(str(ss_object['veg_frac'])[1:-1]))
    f.write("building_scale = {}\n".format(str(ss_object['building_scale'])[1:-1]))
    f.write("veg_scale = {}\n".format(str(ss_object['veg_scale'])[1:-1]))
    f.write("/\n")
    f.write("\n")

    f.write("&roof\n")
    f.write("sfr_roof = {}\n".format(str(ss_object['sfr_roof'])[1:-1]))
    f.write("tin_roof = {}\n".format(str(ss_object['tin_roof'])[1:-1]))
    f.write("alb_roof = {}\n".format(str(ss_object['alb_roof'])[1:-1]))
    f.write("emis_roof = {}\n".format(str(ss_object['emis_roof'])[1:-1]))
    f.write("state_roof = {}\n".format(str(ss_object['state_roof'])[1:-1]))   
    f.write("statelimit_roof = {}\n".format(str(ss_object['statelimit_roof'])[1:-1]))      
    f.write("wetthresh_roof = {}\n".format(str(ss_object['wetthresh_roof'])[1:-1]))  
    f.write("soilstore_roof = {}\n".format(str(ss_object['soilstore_roof'])[1:-1]))     
    f.write("soilstorecap_roof = {}\n".format(str(ss_object['soilstorecap_roof'])[1:-1])) 
    f.write("roof_albedo_dir_mult_fact(1,:) = {}\n".format(str(ss_object['roof_albedo_dir_mult_fact'])[1:-1]))    
    
    
    # TODO make a nested loop based on nlayers
    for j in range(1,ss_object['nlayer'] + 1):
        f.write("dz_roof(" + str(j) + ",:) = {}\n".format(str(ss_object['dz_roof' + str(j)])[1:-1])) 
        f.write("k_roof(" + str(j) + ",:) = {}\n".format(str(ss_object['k_roof' + str(j)])[1:-1])) 
        f.write("cp_roof(" + str(j) + ",:) = {}\n".format(str(ss_object['cp_roof' + str(j)])[1:-1])) 

    f.write("/\n")
    f.write("\n")

    f.write("&wall\n") 
    f.write("sfr_wall = {}\n".format(str(ss_object['sfr_wall'])[1:-1])) 
    f.write("tin_wall = {}\n".format(str(ss_object['tin_wall'])[1:-1]))   
    f.write("alb_wall = {}\n".format(str(ss_object['alb_wall'])[1:-1])) 
    f.write("emis_wall = {}\n".format(str(ss_object['emis_wall'])[1:-1]))   
    f.write("state_wall = {}\n".format(str(ss_object['state_wall'])[1:-1])) 
    f.write("statelimit_wall = {}\n".format(str(ss_object['statelimit_wall'])[1:-1]))      
    f.write("wetthresh_wall = {}\n".format(str(ss_object['wetthresh_wall'])[1:-1])) 
    f.write("soilstore_wall = {}\n".format(str(ss_object['soilstore_wall'])[1:-1]))    
    f.write("soilstorecap_wall = {}\n".format(str(ss_object['soilstorecap_wall'])[1:-1])) 
    f.write("wall_specular_frac(1,:) = {}\n".format(str(ss_object['wall_specular_frac'])[1:-1]))    
    
    # TODO make a nested loop based on nlayers
    for j in range(1,ss_object['nlayer'] + 1):
        f.write("dz_wall(" + str(j) + ",:) = {}\n".format(str(ss_object['dz_wall' + str(j)])[1:-1])) 
        f.write("k_wall(" + str(j) + ",:) = {}\n".format(str(ss_object['k_wall' + str(j)])[1:-1])) 
        f.write("cp_wall(" + str(j) + ",:) = {}\n".format(str(ss_object['cp_wall' + str(j)])[1:-1])) 

    f.write("/\n")
    f.write("\n")

    f.write("&surf\n") 
    f.write("tin_surf = {}\n".format(str(ss_object['tin_surf'])[1:-1]))

    f.write("dz_surf(1,:) = {}\n".format(str(ss_object['dz_surf_paved'])[1:-1]))
    f.write("k_surf(1,:) = {}\n".format(str(ss_object['k_surf_paved'])[1:-1]))
    f.write("cp_surf((1,:) = {}\n".format(str(ss_object['cp_surf_paved'])[1:-1]))

    f.write("dz_surf(2,:) = {}\n".format(str(ss_object['dz_surf_buildings'])[1:-1]))
    f.write("k_surf(2,:) = {}\n".format(str(ss_object['k_surf_buildings'])[1:-1]))
    f.write("cp_surf((2,:) = {}\n".format(str(ss_object['cp_surf_buildings'])[1:-1]))

    f.write("dz_surf(3,:) = {}\n".format(str(ss_object['dz_surf_evergreen'])[1:-1]))
    f.write("k_surf(3,:) = {}\n".format(str(ss_object['k_surf_evergreen'])[1:-1]))
    f.write("cp_surf((3,:) = {}\n".format(str(ss_object['cp_surf_evergreen'])[1:-1]))

    f.write("dz_surf(4,:) = {}\n".format(str(ss_object['dz_surf_decid'])[1:-1]))
    f.write("k_surf(4,:) = {}\n".format(str(ss_object['k_surf_decid'])[1:-1]))
    f.write("cp_surf((4,:) = {}\n".format(str(ss_object['cp_surf_decid'])[1:-1]))

    f.write("dz_surf(5,:) = {}\n".format(str(ss_object['dz_surf_grass'])[1:-1]))
    f.write("k_surf(5,:) = {}\n".format(str(ss_object['k_surf_grass'])[1:-1]))
    f.write("cp_surf((5,:) = {}\n".format(str(ss_object['cp_surf_grass'])[1:-1]))

    f.write("dz_surf(6,:) = {}\n".format(str(ss_object['dz_surf_baresoil'])[1:-1]))
    f.write("k_surf(6,:) = {}\n".format(str(ss_object['k_surf_baresoil'])[1:-1]))
    f.write("cp_surf((6,:) = {}\n".format(str(ss_object['cp_surf_baresoil'])[1:-1]))

    f.write("dz_surf(7,:) = {}\n".format(str(ss_object['dz_surf_water'])[1:-1]))
    f.write("k_surf(7,:) = {}\n".format(str(ss_object['k_surf_water'])[1:-1]))
    f.write("cp_surf((7,:) = {}\n".format(str(ss_object['cp_surf_water'])[1:-1]))

    f.write("/\n")
    f.write("\n")
    f.close()

    return ss_file_path

# def read_GridLayout_file(refdir, fname):
#     ss_file_path = os.path.join(refdir, fname + ".uwg")
#     # f = open(uwg_file_path, "r")

#     ssDict = {}
#     skiptype = 0
#     skipcount = 0
#     trafficlist = []
#     bldlist = []
#     l1 = []
#     l2 = []
#     l3 = []

#     with open(ss_file_path) as file:
#         # next(file)
#         for line in file:
#             if line[0:7] == 'SchTraf':
#                 skiptype = 1
#             if line[0:4] == 'bld,':
#                 skiptype = 2
#             if skiptype == 0:
#                 if line[0] == '#' or line == '\n': # empty line or comment
#                     test = 4 
#                 else: # regular input
#                     a = line.find(',')
#                     if line[0:a] == 'zone':
#                         ssDict[line[0:a]] = line[a +1: len(line) - 2]
#                     elif line[-3:] == ',,\n':
#                         ssDict[line[0:a]] = None
#                     else:
#                         ssDict[line[0:a]] = float(line[a +1: len(line) - 2])
#             elif skiptype == 1: # Traffic
#                 if skipcount >= 1:
#                     letter_list = line.split(",")
#                     floats_list = []
#                     for item in letter_list:
#                         if item == '\n':
#                             test = 4
#                         else:
#                             floats_list.append(float(item))
#                     trafficlist.append(floats_list)
#                 skipcount += 1
#                 if skipcount == 4:
#                     skipcount = 0
#                     skiptype = 0
#                     ssDict['SchTraffic'] = trafficlist
#             elif skiptype == 2: #Buildings
#                 if skipcount >= 1:
#                     letter_list = line.split(",")
#                     l1.append(letter_list[0])
#                     l2.append(letter_list[1])
#                     l3.append(float(letter_list[2]))
                    
#                     # if skipcount < 3:
#                     #     for item in letter_list:
#                     #         if item == '\n':
#                     #             test = 4
#                     #         else:
#                     #             floats_list.append(item)
#                     #     bldlist.append(floats_list)
#                     # else:
#                     #     for item in letter_list:
#                     #         if item == '\n':
#                     #             test = 4
#                     #         else:
#                     #             floats_list.append(float(item))
#                     #     bldlist.append(floats_list)
#                 skipcount += 1
#                 if skipcount == 17:
#                     skipcount = 0
#                     skiptype = 0
#                     bldlist.append(l1)
#                     bldlist.append(l2)
#                     bldlist.append(l3)
#                     ssDict['bld'] = bldlist
                    
#     return ssDict







