import os
import shutil

cwd = os.getcwd()
dir = os.path.join(cwd, '01_Submitted/Layers/Images')


for filename in os.listdir(dir):
    f = os.path.join(dir, filename)
    if f.find('2_NoDAPI') != -1 : #filter by
        #copy file to other
        shutil.copyfile(f, os.path.join(os.path.join(cwd, 'ManualAnalysis/Registration'), filename))
