#!/usr/bin/env python

import os
import h5py
import subprocess

new_profile_filename_base = 'input.profiles.h5'
old_profile_filename_base = 'input.profiles.h5_old'

new_grid_filename_base = 'psiAHat.h5'
old_grid_filename_base = 'psiAHat.h5_old'




abs_precision = "0.000001"
#abs_precision = str(abs_precision)

cmd = "h5diff"
difftype = "-d"

seperator = "==============="

differences = 0
dirs_with_differences = []

dir_of_this_file = os.path.dirname(__file__)
if dir_of_this_file =='':
    dir_of_this_file = '.'
subdirs = next(os.walk(dir_of_this_file))[1]

for dir in subdirs:
    dir = os.path.join(dir_of_this_file,dir)
    new_profile_filename = os.path.join(dir,new_profile_filename_base)
    old_profile_filename = os.path.join(dir,old_profile_filename_base)

    new_grid_filename = os.path.join(dir,new_grid_filename_base)
    old_grid_filename = os.path.join(dir,old_grid_filename_base)
    
    print seperator + " dir: " + dir + " " + seperator

    #check input.profiles.h5
    print ">removing " + new_profile_filename
    subprocess.call(["rm",new_profile_filename])
    print ">removing " + new_grid_filename
    subprocess.call(["rm",new_grid_filename])

    print ">generating " + new_profile_filename
    subprocess.call(["python","generate_profiles.py"],cwd=dir)
    
    print cmd +" " + difftype + " " + abs_precision + " " + new_profile_filename_base + " " + old_profile_filename_base +":"
    out = subprocess.call([cmd, difftype,abs_precision,new_profile_filename,old_profile_filename])
    if out != 0:
        print "*DIFFERENCES FOUND!*"
        differences += 1
        dirs_with_differences.append(dir)

    #check psiAHat.h5
    if os.path.isfile(old_grid_filename):
        print cmd + " " + difftype + " " + abs_precision + " " + new_grid_filename_base + " " + old_grid_filename_base +":"
        out = subprocess.call([cmd, difftype,abs_precision,new_grid_filename,old_grid_filename])
        if out != 0:
            print "*DIFFERENCES FOUND!*"
            differences += 1
            dirs_with_differences.append(dir)
        

print seperator

if differences > 0:
    print str(differences)  + " differences found, using " + difftype + " " + abs_precision + " ,in directories " + repr(dirs_with_differences)
else:
    print "No differences found!"
