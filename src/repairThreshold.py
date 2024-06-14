import matlab.engine

eng = matlab.engine.start_matlab()

eng.repairThreshold(nargout=0)


# a = eng.imread('/Users/rowelab/Desktop/kevin/KOMP_matlab/thresholdTest/CCC_E11_hF_ML3_s2_2_shift3_original.jpg')
# #  eng.imshow(a)

# # S = eng.substruct('()',{':',':',1})
# # results in a python dict without duplicates {':',1} ????????? someone fix this pls
# b = eng.subsref(a,eng.substruct('()',{':',':',1}))

# # th = eng.graythresh(b)*255
# th = eng.graythresh(a)*255
# c = eng.gt(b,th)
# eng.imwrite(c,'/Users/rowelab/Desktop/kevin/KOMP_matlab/thresholdTest/CCC_E11_hF_ML3_s2_2_shift3_new.jpg', nargout=0)