function control_threshold_batch(bt, bn, section)

a=textread('dir_info.txt','%s', 2);
home_dir=a{1,:};
phr_temp=a{2,:};

delimeter='\';
exp_type='F';

direct=[home_dir,'01_Submitted',delimeter,'Layers',delimeter];        % segment tetraycline with 2*std_dev 
exp_name={[phr_temp,'_h',exp_type,'_F'];[phr_temp,'_h',exp_type,'_M']};        % gene name, exp_name, (female, male)
labels={'green','green';'red','red'};
GFP=0;

addpath(home_dir);
%test(dire)
[bt, bn, section]

find_thresholding_values_KOMP_test_for_BioHead_1_batch(direct, exp_name, exp_type, labels, GFP, bt, bn, section, delimeter)

% exit
