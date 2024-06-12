function control_rotate


a=textread('dir_info.txt','%s', 2);
home_dir=a{1,:};
phr_temp=a{2,:};
if ~isempty(strfind(home_dir,'\'))
    delimeter='\';
elseif ~isempty(strfind(home_dir,'/'))
    delimeter='/';
end

exp_type='F';
bone_type=home_dir(end-1);
% delimeter='\';

direct=[home_dir,'01_Submitted',delimeter,'Layers',delimeter,'Images',delimeter];
gene_exp=[phr_temp,'_h',exp_type];




addpath(home_dir);


rotate_main_axis_1(direct,gene_exp,exp_type)

% exit