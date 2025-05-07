import numpy as np
# import pandas as pd
from osgeo import gdal, osr
from osgeo.gdalconst import *

def wallOfInterest(woilyr, woi_field, minx, miny, scale, rows, outputDir, dsm_minx, x_size, dsm_miny, y_size):

    # header = 'yyyy id   it imin dectime altitude azimuth    Ta' 

    idx = woilyr.fields().indexFromName(woi_field[0])
    numfeat = woilyr.featureCount()
    poiname = []
    poisxy = np.zeros((numfeat, 3)) - 999
    ind = 0
    for f in woilyr.getFeatures():  # looping through each POI
        y = f.geometry().centroid().asPoint().y()
        x = f.geometry().centroid().asPoint().x()

        poiname.append(f.attributes()[idx])
        poisxy[ind, 0] = ind
        # poisxy[ind, 1] = np.round((x - minx) * scale)

        xcoordinate = np.floor((x - dsm_minx) * scale)
        ycoordinate = np.floor((dsm_miny - y) * scale)

        # if miny >= 0:
        #     poisxy[ind, 2] = np.round((miny + rows * (1. / scale) - y) * scale)
        # else:
        #     poisxy[ind, 2] = np.round((miny + rows * (1. / scale) - y) * scale)
            # print('y = ' + str(np.round((miny + rows * (1. / scale) - y) * scale)))

        poisxy[ind, 1] = xcoordinate
        poisxy[ind, 2] = ycoordinate

        ind += 1

    # for k in range(0, poisxy.shape[0]):
    #     poi_save = []  # np.zeros((1, 33))
    #     data_out = outputDir + '/POI_' + str(poiname[k]) + '.txt'
    #     np.savetxt(data_out, poi_save,  delimiter=' ', header=header, comments='')

    return poisxy, poiname

# def woiHeader():

def fillWallOfInterest(voxelTable, voxelHeight, woisxy, woiname, outputDir, i):

    if not woisxy is None:
        for k in range(0, woisxy.shape[0]):
            tempTable = voxelTable.loc[((voxelTable['ypos'] == woisxy[k, 2]) & (voxelTable['xpos'] == woisxy[k, 1]) & voxelTable['voxelHeight'] == voxelHeight)].copy()
            output_vars = np.array(['wallTemperature', 'K_in', 'L_in', 'wallTemperatureSolweigOld', 
                                    'Lwallsun', 'Lwallsh', 'Lrefl', 'Lveg', 'Lground', 'Lsky',
                                    'esky', 'F_sh', 'wallSun',
                                    'svfbu', 'svfveg', 'svfaveg'])
            
            counter = 0
            for temp_var in output_vars:
                temp_out = tempTable[temp_var].to_numpy()
                if counter == 0:
                    temp_all = temp_out.copy()
                else:
                    temp_all = np.concatenate([temp_all, temp_out])

            # temp_wall = tempTable['wallTemperature'].to_numpy()
            # K_in = tempTable['K_in'].to_numpy()
            # L_in = tempTable['L_in'].to_numpy()
            # Ts_SolweigOld = tempTable['wallTemperatureSolweigOld'].to_numpy()

            # Lwallsun = tempTable['Lwallsun'].to_numpy()
            # Lwallsh = tempTable['Lwallsh'].to_numpy()
            # Lrefl = tempTable['Lrefl'].to_numpy()
            # Lveg = tempTable['Lveg'].to_numpy()
            # Lground = tempTable['Lground'].to_numpy()
            # Lsky = tempTable['Lsky'].to_numpy()
            # esky_wall = tempTable['esky'].to_numpy()
            # F_sh_wall = tempTable['F_sh'].to_numpy()
            # wallSun_wall = tempTable['wallSun'].to_numpy()

            # svfbu_wall = tempTable['svfbu'].to_numpy()
            # svfveg_wall = tempTable['svfveg'].to_numpy()
            # svfaveg_wall = tempTable['svfaveg'].to_numpy()

            # # temp_all = np.concatenate([temp_wall, K_in, L_in, Ts_Erik, Ts_SolweigOld])
            # temp_all = np.concatenate([temp_wall, K_in, L_in, Ts_SolweigOld, Lwallsun, Lwallsh, Lrefl, Lveg, Lground, Lsky, esky_wall, F_sh_wall, wallSun_wall,
            #                             svfbu_wall, svfveg_wall, svfaveg_wall])
                    
            wall_data = np.zeros((1, 7 + temp_all.shape[0]))
            # Part of file name (wallid), i.e. WOI_wallid.txt
            woiwallId = voxelTable.loc[((voxelTable['ypos'] == woisxy[k, 2]) & (voxelTable['xpos'] == woisxy[k, 1])), 'wallId'].to_numpy()[0]                    
            # print('wall id = ' + str(woiwallId) + ' and direction = ' + str(woiname[k]))
            data_out = outputDir + '/WOI_' + str(woiname[k]) + '_height' + str(voxelHeight) + '.txt'                    
            if i == 0:
                print('wall id = ' + str(woiwallId) + ' and direction = ' + str(woiname[k]))
                # Output file header
                header = 'yyyy id   it imin dectime Ta  SVF '
                voxelHeader = ''
                for temp_header in output_vars:
                    header += '    ' + temp_header
                #header = header + voxelHeader

                # Part of file name (wallid), i.e. WOI_wallid.txt
                # woiname = voxelTable.loc[((voxelTable['ypos'] == woisxy[k, 2]) & (voxelTable['xpos'] == woisxy[k, 1])), 'wallId'].to_numpy()[0]
                woi_save = []  # 
                # data_out = outputDir + '/WOI_' + str(woiname) + '.txt'
                np.savetxt(data_out, woi_save,  delimiter=' ', header=header, comments='')                        
            # Fill wall_data with variables
            wall_data[0, 0] = YYYY[0][i] 
            wall_data[0, 1] = jday[0][i]
            wall_data[0, 2] = hours[i]
            wall_data[0, 3] = minu[i]
            wall_data[0, 4] = dectime[i]
            wall_data[0, 5] = Ta[i]
            wall_data[0, 6] = svf[int(woisxy[k, 2]), int(woisxy[k, 1])]
            wall_data[0, 7:] = temp_all
            # wall_data[0, 7:] = temp_wall

            # Num format for output file data
            # woi_numformat = '%d %d %d %d %.5f %.2f %.2f' + ' %.2f' * temp_wall.shape[0]
            woi_numformat = '%d %d %d %d %.5f %.2f %.2f' + ' %.2f' * temp_all.shape[0]
            # Open file, add data, save
            f_handle = open(data_out, 'ab')
            np.savetxt(f_handle, wall_data, fmt=woi_numformat)
            f_handle.close()    