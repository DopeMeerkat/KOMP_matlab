function find_growth_plate_2%(direct, exp_name, exp_type)

% direct='E:\seh00004\Het_KOMP\CCC_E11_hF\01_Submitted\Layers\';
% exp_name={'CCC_E11_hF_F';'CCC_E11_hF_M'};
% exp_type='F';

a=textread('dir_info.txt','%s', 2);
home_dir=a{1,:};
phr_temp=a{2,:};
if ~isempty(strfind(home_dir,'\'))
    delimeter='\';
elseif ~isempty(strfind(home_dir,'/'))
    delimeter='/';
end
exp_type=home_dir(end-1);
direct=[home_dir,'01_Submitted',delimeter,'Layers',delimeter];
exp_name={[phr_temp,'_h',exp_type,'_F'];[phr_temp,'_h',exp_type,'_M']};        % gene name, exp_name, (female, male)

%  add_ind=[                  % add [bt, bn, section] to the existing growth_plate
%      1 7 3
%     2 2 3;
%     2 6 2
%     2 7 1;
%     2 7 2;
%      ];

% labels={'green','green';'red','red'};      % first label, second label
% GFP=0;            % no GFP
%     1;            % there is GFP
% thresh_growth=0;  % find threshold + find growth_plate (default)
%               1;  % find threshold only
%               2;  % find growth_plate only



if strcmp(direct(end),delimeter)~=1
    direct=[direct,delimeter];
end
direc=[direct,'Images',delimeter];

eval(['load ''',direct, 'info'''])
sample_per_bone=0;
if exist('add_ind')
    bone_type=unique(add_ind(:,1));
    no_bone_type=length(bone_type);
    for i=1:no_bone_type
        bt=bone_type(i);
        eval(['bone_number',num2str(bt),'=add_ind(add_ind(:,1)==',num2str(bt),',2);'])
        eval(['bone_number',num2str(bt),'=num2str(bone_number',num2str(bt),');'])
        eval(['bone_number',num2str(bt),'=mat2cell(bone_number',num2str(bt),',[ones(size(bone_number',num2str(bt),',1),1)], [1])'])
    end
else
    no_bone_type=length(exp_name);
    bone_type=[1:no_bone_type];
    for i=1:no_bone_type
        bt=bone_type(i);
        eval(['test=isempty(info',num2str(bt),');'])
        if test==1
            continue
        end
        eval(['temp=info',num2str(bt),'.im;'])
    %     temp1=sum(cell2mat(temp(:,2:end)),2);
        bone_no=temp(:,1);
        sample_per_bone(bt)=max(max(sample_per_bone,size(temp,2)-1));
        eval(['bone_number',num2str(bt),'=bone_no;'])
    end
end


for i=1:no_bone_type
    bt=bone_type(i);
    eval(['temp=exp_name{',num2str(bt),'};'])
    if strcmp(temp(end),'_')~=1
        temp=[temp,'_'];
    end
    eval(['p',num2str(bt),'=temp;'])
end

%%%%
%
%   July 30, 2013
%
%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          growth_plate_1           --------        growth_plate_n
%          |--------------| section 1~3        
%        |--------------| |                        
%      |--------------| | |
%      |b    g     b  | | |                       
%      |o    r     o  | | |                       
%      |n    o     t  | | |                       
%      |e    w     t  | | |                       
%      |#    t     o  | | |                        
%      |     h     m  | | |                        
%      |              | | |                        
%      |     p        | | |                        
%      |     l        | | |                       
%      |     a        | | |                      
%      |     t        | |-|
%      |     e        |-|
%      |--------------| 
%

for i=1:no_bone_type                % bone type (AAA_E01_F, AAA_E01_M)
    bt=bone_type(i);
    eval(['test=exist(''bone_number',num2str(bt),''');'])
    if test==0
        continue
    end
    eval(['bone_number=bone_number',num2str(bt),';']);
    phr1=exp_name{bt,:};             %['AAA_E01_F'];
    
    if ~isempty(dir([direct,'growth_plate',num2str(bt),'.mat']))
        load([direct,'growth_plate',num2str(bt),'.mat'])
    else
        eval(['growth_plate',num2str(bt),'=bone_number;']);   % bottom of growth plate
    end

    for j=1:length(bone_number)     % bone number
        if exist('add_ind')
            temp=add_ind(add_ind(:,1)==bt,:);
            sample_per_bone_add_ind=temp(temp(:,2)==str2num(bone_number{j}),3);
            sample_per_bone(bt)=length(sample_per_bone_add_ind);
        end
        for s=1:sample_per_bone(bt)                   % section number
            if exist('add_ind')
                ss=sample_per_bone_add_ind(s);
            else
                ss=s;
            end
            if strcmp(exp_type,'V')
                eval(['a=dir(''',direc,phr1, bone_number{j},'_h',exp_type,'_s',num2str(ss),'_shift3_NoDAPI.jpg'');']);
            elseif strcmp(exp_type,'F')
                eval(['a=dir(''',direc,phr1, 'L', num2str(ss),'_s',bone_number{j},'_shift3_NoDAPI.jpg'');']);
            end
            if isempty(a)
                continue
            end
            if strcmp(exp_type,'V')
                eval(['a=imread(''',direc,phr1, bone_number{j},'_h',exp_type,'_s',num2str(ss),'_shift3_NoDAPI.jpg'');']);
            elseif strcmp(exp_type,'F')
                eval(['a=imread(''',direc,phr1, 'L', num2str(ss),'_s',bone_number{j},'_shift3_NoDAPI.jpg'');']);
            end
%             imshow(imrotate(a,90-83.5528))
            imshow(a)
            title(['',phr1,bone_number{j},'-section',num2str(ss),''])
            axis on;grid
            pan on
            [x y]=ginput(2);
%             eval(['growth_plate',num2str(i),'{',num2str(j),',1,',num2str(s),'}=''',bone_number{j},''';']);   % bottom of growth plate
%             eval(['growth_plate',num2str(i),'{',num2str(j),',2,',num2str(ss),'}=round(y(1));']);   % bottom of growth plate
%             eval(['growth_plate',num2str(i),'{',num2str(j),',3,',num2str(ss),'}=round(y(2));']);   % bottom of growth plate
            eval(['growth_plate',num2str(bt),'{',bone_number{j},',2,',num2str(ss),'}=round(y(1));']);   % bottom of growth plate
            eval(['growth_plate',num2str(bt),'{',bone_number{j},',3,',num2str(ss),'}=round(y(2));']);   % bottom of growth plate
        end
    end
    eval(['save ''',direct, 'growth_plate',num2str(bt),''' growth_plate',num2str(bt)'])
end
close all