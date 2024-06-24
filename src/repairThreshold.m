function repairThreshold()
offset = 0;
    a = imread('/Users/rowelab/Desktop/kevin/KOMP_matlab/thresholdTest/CCC_E11_hF_ML3_s2_2_shift3_original.jpg');
b = a(:,:,1);
th = graythresh(b)*255;
th = th + offset;
c = b>th;
imwrite(c,'/Users/rowelab/Desktop/kevin/KOMP_matlab/thresholdTest/CCC_E11_hF_ML3_s2_2_shift3_new.jpg');
end