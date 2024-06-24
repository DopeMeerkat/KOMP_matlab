import os
import shutil

cwd = os.getcwd()
dir = os.path.join(cwd, '02_Analyzed')


for filename in os.listdir(dir):
    f = os.path.join(dir, filename)
    if filename.find('_ro') != -1 : #filter by
        #copy file to other
        shutil.copyfile(f, os.path.join(os.path.join(cwd, 'ManualAnalysis/ROI'), filename))

