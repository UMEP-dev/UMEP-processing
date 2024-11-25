import os

'''
input: 
TargetDict: dictonary with all relevant inputs
refdir = outputfolder      
'''

def write_config_file(targetDict, refdir):
    config_file_path = os.path.join(refdir, "config.ini")
    f = open(config_file_path, "w")

    f.write("\n")
    f.write("[DEFAULT]\n")
    f.write("\n")
    f.write("# =================================================\n")
    f.write("# CONTROL FILE\n")
    f.write("# =================================================\n")
    f.write("### INPUTS ###\n")
    f.write("# output path\n")
    f.write("work_dir={}\n".format(targetDict['work_dir']))
    f.write("# parameters json file\n")
    f.write("para_json_path={}\n".format(targetDict['para_json_path']))
    f.write("# site name (string)\n")
    f.write("site_name={}\n".format(targetDict['site_name']))
    f.write("# run name (string)\n")
    f.write("run_name={}\n".format(targetDict['run_name']))
    f.write("# input meteorolgical file (i.e. forcing file)\n")
    f.write("inpt_met_file={}\n".format(targetDict['inpt_met_file']))    
    f.write("# input land cover data file\n")
    f.write("inpt_lc_file={}\n".format(targetDict['inpt_lc_file'])) 
    f.write("# format of datetime in input met files\n")
    f.write("date_fmt={}\n".format(targetDict['date_fmt'])) 
    f.write("# time step (minutes)\n")
    f.write("timestep={}\n".format(targetDict['timestep'])) 
    f.write("\n")
    f.write("include roofs={}\n".format(targetDict['include roofs']))
    f.write("mod_ldwn={}\n".format(targetDict['mod_ldwn']))  
    f.write("domainDim={}\n".format(targetDict['domaindim']))     
    f.write("latEdge={}\n".format(targetDict['latedge']))       
    f.write("lonEdge={}\n".format(targetDict['lonedge']))      
    f.write("latResolution={}\n".format(targetDict['latresolution']))   
    f.write("lonResolution={}\n".format(targetDict['lonresolution']))     
    f.write("\n")
    f.write("#----------------------------------\n")
    f.write("# dates\n")
    f.write("#----------------------------------\n")
    f.write("# year,month,day,hour	# start date for simulation (should be a minimum of 24 hours prior to date1)\n")
    f.write("date1a={}\n".format(targetDict['date1a'])) 
    f.write("# year,month,day,hour	# the date/time for period of interest (i.e. before this will not be saved)\n") 
    f.write("date1={}\n".format(targetDict['date1']))
    f.write("# year,month,day,hour	# end date for validation period\n") 
    f.write("date2={}\n".format(targetDict['date2']))       
    f.write("\n")

    f.close()
