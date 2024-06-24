#!/usr/local/bin/python3
import os
import shutil

#creating directories
cwd = os.getcwd()

shutil.copytree('/Users/rowelab/Desktop/kevin/KOMP_matlab/src', os.path.join(cwd, 'CCC_EX_hF'), copy_function= shutil.copy, dirs_exist_ok = True)

try:
    os.mkdir(os.path.join(cwd, 'CCC_EX_hF/01_Submitted'))
except:
    print(cwd)
    print('Directory already exists')
try:
    os.mkdir(os.path.join(cwd, 'CCC_EX_hF/02_Analyzed'))
except:
    print('Directory already exists')
try:
    os.mkdir(os.path.join(cwd, 'CCC_EX_hF/03_CalculatedData'))
except:
    print('Directory already exists')
try:
    os.mkdir(os.path.join(cwd, 'CCC_EX_hF/ManualAnalysis'))
except:
    print('Directory already exists')
try:
    os.mkdir(os.path.join(cwd, 'CCC_EX_hF/ManualAnalysis/Threshold2'))
except:
    print('Directory already exists')
try:
    os.mkdir(os.path.join(cwd, 'CCC_EX_hF/ManualAnalysis/Registration'))
except:
    print('Directory already exists')

try:
    os.mkdir(os.path.join(cwd, 'CCC_EX_hF/ManualAnalysis/ROI'))
except:
    print('Directory already exists')

#copying code
# dir = '/Users/rowelab/Desktop/kevin/KOMP_matlab/src/'
# for filename in os.listdir(dir): 
#     try:
#         shutil.copyfile(dir+filename, os.path.join(cwd, filename))
#     except:
#         print('File already exists')

