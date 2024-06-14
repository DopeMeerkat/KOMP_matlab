# Dependencies
- python3 version > 3.11
- matlab
- matlabengine `$ python -m pip install matlabengine==24.1.2`

# Instructions for Mac
- change the path for where the source is located in `setup.py`
- run `setup.command`, this will create a directory in the current folder named `CCC_EX_hF`
- rename the folder so that it matches the images. i.e. images are named `CCC_E11_hF_ML3_s2_2_shift3.jpg` so folder should be named `CCC_E11_hF`
- open the `00_Scripts` folder found in the directory just created
- run `01_registration.command`
- after it is done running, run `02_filterRegistration.command`
- open the folder `ManualAnalysis/Registration` and verify that the samples are aligned correctly
- run `03_thresholding.command`
- after it is done running, run `04_filterThresholds.command`
- open `ManualAnalysis/Threshold2` folder and look through the images, make changes to any that needed touch ups and add a f to the filename, mark any unusable ones with u, and mark ones that need to be redone with r
i.e. `CCC_E11_hF_ML3_s2_2_shift3.jpg` -> `CCC_E11_hF_ML3_s2_2_shift3r.jpg`


