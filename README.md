# Dependencies
- python3 version > 3.11
- matlab
- matlabengine `$ python -m pip install matlabengine==24.1.2`

# Instructions for Mac

## Step 1:
- run `setup.command`, this will create a directory in the current folder named `CCC_EX_hF`
- rename the folder so that it matches the images. i.e. images are named `CCC_E11_hF_ML3_s2_2_shift3.jpg` so folder should be named `CCC_E11_hF`
- open the `00_Scripts` folder found in the directory just created
- run `01_registration.command`

## Step 2:
- after registration is finished running, open the folder `ManualAnalysis/Registration` and verify that the samples are aligned correctly
- run `02_thresholding.command`

## Step 3:
- after thresholding is finished, open `ManualAnalysis/Threshold2` folder and look through the images, make changes to any that needed touch ups and add a f to the filename, mark any unusable ones with u, and mark ones that need to be redone with r
i.e. `CCC_E11_hF_ML3_s2_2_shift3.jpg` -> `CCC_E11_hF_ML3_s2_2_shift3f.jpg`
- run `03_rotate.command`

## Step 4:



