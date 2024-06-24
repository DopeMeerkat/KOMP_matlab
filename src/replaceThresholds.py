import os
import shutil

cwd = os.getcwd()
dir = os.path.join(cwd, 'ManualAnalysis/Threshold2')


for filename in os.listdir(dir):
    f = os.path.join(dir, filename)
    # f = fixed
    # r = redo
    # u = unusable
    if filename.find('shift3f') != -1 : #filter by
        #copy file to other
        # remove f
        newname = filename[:-5] + filename[-4:]
        # print(newname)


        shutil.move(f, os.path.join(os.path.join(cwd, '01_Submitted/Layers/Threshold'), newname))

        # needs testing
        # destinationPath = os.path.join(os.path.join(cwd, '01_Submitted/Layers/Threshold'), newname)
        # if os.path.exists(destinationPath):
        #     os.remove(destinationPath)
        # shutil.copyfile(f, destinationPath)

