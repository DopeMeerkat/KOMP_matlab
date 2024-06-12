function control_analysis(bt, bn, section)

a=textread('dir_info.txt','%s', 2);
home_dir=a{1,:};
phr_temp=a{2,:};
if ~isempty(strfind(home_dir,'\'))
    delimeter='\';
elseif ~isempty(strfind(home_dir,'/'))
    delimeter='/';
end

% bone_type={'F'};
bone_type={home_dir(end-1)};
% delimeter='\';

root_di=[home_dir,'01_Submitted',delimeter,'Layers',delimeter];
direct=[home_dir,'01_Submitted',delimeter,'Layers',delimeter,'Images',delimeter];        % segment tetraycline with 2*std_dev 
out_dir=[home_dir,'02_Analyzed',delimeter];       % segment tetraycline with 2*std_dev 

phr={[phr_temp,'_h',cell2mat(bone_type),'_F'];[phr_temp,'_h',cell2mat(bone_type),'_M']};

%%%%%%%%%%%%%%%%%%%%%
%
%	June, 7th 2018
%
% labels={'green','green';'red','red'};

RG_label=true; % 1st label : green,  2nd label : red --> Noemal labels
% RG_label=false; % 1st label : red,  2nd label : green --> Reversed labels

%
%	June, 7th 2018
%
%%%%%%%%%%%%%%%%%%%%%

GFP=0;
auto_exp=0;
mouse = 1; % --> mouse(default)
%       0; % --> rat 

addpath(home_dir);
%test(dire)
[bt, bn, section]

analize_KOMP_for_BioHead_6(root_di, direct, out_dir, phr, delimeter, bt, bn, section, GFP, RG_label, mouse)

% exit