#!/usr/local/bin/env python3
import os
import shutil

#creating directories
cwd = os.getcwd()
try:
    os.mkdir(os.path.join(cwd, '01_Submitted'))
except:
    print('Directory already exists')
try:
    os.mkdir(os.path.join(cwd, '02_Analyzed'))
except:
    print('Directory already exists')
try:
    os.mkdir(os.path.join(cwd, '03_CalculatedData'))
except:
    print('Directory already exists')
try:
    os.mkdir(os.path.join(cwd, 'Threshold2'))
except:
    print('Directory already exists')

#copying code
dir = '/Users/rowelab/Desktop/kevin/matlab/src/'
for filename in os.listdir(dir): 
    try:
        shutil.copyfile(dir+filename, os.path.join(cwd, filename))
    except:
        print('File already exists')