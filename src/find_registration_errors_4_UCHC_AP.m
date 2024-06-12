function find_registration_errors_4_UCHC_AP(varargin)

%
% input : [bt, bn, section]
%          :    :     :
%

[home_dir, Exp_name]=produce_directory_information_file;

cd (home_dir)
tissue_type=home_dir(end-1);
delimeter='\';

if nargin==0    % first Round
    er=c_registration;
else
    er=c_registration(varargin{1});
end

direct=[home_dir,'01_Submitted',delimeter,'Layers',delimeter,'Images',delimeter];
reg_info_files=[direct,'*_reg.mat'];
aa=dir(reg_info_files);
% result=cell(length(aa),2);
for i=1:999
    if i<10
        c=['00',num2str(i)];
    elseif i>=10 & i<100
        c=['0',num2str(i)];
    elseif i>=100 & i<1000
        c=num2str(i);
    end
    
    ResultFileName=['AP_RegistrationResult',c,'.txt'];
    aa=dir(ResultFileName);
    if ~isempty(aa)
        continue
    else
        break
    end
end
fileID = fopen(ResultFileName,'w');
fprintf(fileID,'%s\n\n',[[Exp_name,home_dir(end-3:end-1)], ' AP Registration Results']);
aa=dir([home_dir,'01_Submitted\Layers\Images\*_reg.mat']);
for i=1:length(aa)
    n=aa(i).name;
    ind=strfind(n,'_reg.mat');
    temp=MakeRegistrationReport([home_dir,'01_Submitted\Layers\Images\',n]);
%     result(i,1)=mat2cell(n(1:ind-1),1,length(n(1:ind-1)));
%     result(i,2)=mat2cell(temp,1,length(temp));

    fprintf(fileID,'%s\t',n(ind-6:ind-1));
    fprintf(fileID,'%s\n',temp);
end

fclose all;

function [direct, Exp_name]=produce_directory_information_file(varargin)
%
% input : overwrite 
%         1 --> overwrite (defalut)
%         0 --> do not overwrite
%
if nargin==0
    overwrite=1;
else
    overwrite=varargin{1};
end

aa=dir('dir_info.txt');
if ~isempty(aa) & overwrite==0
    disp('''dir_info.txt'' exists and do not overwrite')
    return
elseif ~isempty(aa) & overwrite==1
    disp('overwrite ''dir_info.txt''')
elseif isempty(aa) & overwrite==0
    disp('need to produce ''dir_info.txt''.')
    disp('Program will automatically produce ''dir_info.txt''')
else
    disp('produce ''dir_info.txt''')
end

aa=dir;
n=aa(1).folder;
del=strfind(aa(1).folder,'\');
hV=strfind(aa(1).folder,'_hV');
hF=strfind(aa(1).folder,'_hF');
if isempty(hV) & isempty(hF)
    disp('File name should have _hV(Vertebra) or _hF(Femur)')
    return
elseif isempty(hV) & ~isempty(hF)
    Exp_name=n(del(end)+1:hF-1);
elseif ~isempty(hV) & isempty(hF)
    Exp_name=n(del(end)+1:hV-1);
else
    disp('File name should have one of _hV or _hF')
    return
end
% outdata=cell(2,1);
% outdata(1,1)=mat2cell([n,'\'],ones(1,1),ones(1,1));
% outdata(2,1)=mat2cell(Exp_name,ones(1,1),ones(1,1));
% outdata
fileID = fopen('dir_info.txt','w');
direct=[n,'\'];
fprintf(fileID,'%s\n',direct);
fprintf(fileID,'%s',Exp_name);
fclose all;

function [er]=c_registration(varargin)

if nargin==0    % first Round
    n=0;
    for i=1:2
        for j=1:8
            for k=1:3
                n=n+1;
                bt(n)=i;
                bn(n)=j;
                section(n)=k;
            end
        end
    end
else
    redo=varargin{1};
    n=size(redo,1);
    for i=1:n
        bt(i)=redo(i,1);
        bn(i)=redo(i,2);
        section(i)=redo(i,3);
    end
end


% par
% er=zeros(8,3,2);
for i=1:n
    er=control_registration_2_v2channel(bt(i), bn(i), section(i));
    ['Gender=',num2str(bt(i)),', Bone section=', num2str(bn(i)), ', bone number=',num2str(section(i))]
end

function [er]=control_registration_2_v2channel(bt, bn, section)

a=textread('dir_info.txt','%s', 2);
home_dir=a{1,:};
phr_temp=a{2,:};

bone_type=home_dir(end-1);
delimeter='\';

direct=[home_dir,'01_Submitted',delimeter,'Layers',delimeter];
out_direct1=[home_dir,'01_Submitted',delimeter,'Layers',delimeter,'Images',delimeter];
% exp_type={'Female';'Male'};
exp_type={'Female';'Male'};
image_type='jpg';
labels={'green','green';'red','red'};   % first label, second label
name={[phr_temp,'_h',bone_type,'_F'];[phr_temp,'_h',bone_type,'_M']};         % gene name, exp_name, (female, male)

GFP=0;              % if there is GFP, then 1. else 0
mouse=1;    % mouse = 1 (default)                
gutta=0;                                % no pre-registration (gutta is not used) 

addpath(home_dir);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remove excessive empty space
% 2020 June 30

% registration_KOMP_6_v2channel(direct, out_direct1, exp_type, image_type, labels, name, bone_type, delimeter, GFP, bt, bn, section, gutta, mouse)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add Cy5 SafO, Brightfield SafO
% 2020 July 4
% registration_KOMP_7_v2channel(direct, out_direct1, exp_type, image_type, labels, name, bone_type, delimeter, GFP, bt, bn, section, gutta, mouse)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% exit

reg_error_code=registration_KOMP_8_v2channel(direct, out_direct1, exp_type, image_type, labels, name, bone_type, delimeter, GFP, bt, bn, section, gutta, mouse);
if ~isempty(reg_error_code)
    if bt==1
        G=['Female, sample ',num2str(bn),', section ',num2str(section)];
    elseif bt==2
        G=['Male, sample ',num2str(bn),', section ',num2str(section)];
    end
end
if reg_error_code>=50
    S='SO';
    reg_error_code=reg_error_code-50;
elseif reg_error_code>=40 & reg_error_code<50
    S='Cy5';
    reg_error_code=reg_error_code-40;
elseif reg_error_code>=30 & reg_error_code<40
    S='TB';
    reg_error_code=reg_error_code-30;
elseif reg_error_code>=20 & reg_error_code<30
    S='AP';
    reg_error_code=reg_error_code-20;
elseif reg_error_code>=10 & reg_error_code<20
    S='TRAP';
    reg_error_code=reg_error_code-10;
end
switch reg_error_code
    case 0
        er=[G, ': ',S, ': Done']
    case 1
        er=[G, ': ',S, ': Second beads group not found']
    case 2.1
        er=[G, ': ',S, ': beads not found']        
    case 2.2
        er=[G, ': ',S, ': too much rotation']
    case 2.2
        er=[G, ': ',S, ': too big or small']
    case 3
        er=[G, ': ',S, ': beads not found'] 
end

function [reg_error_code]=registration_KOMP_8_v2channel(varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% registration_KOMP_8_v2channel
%
% Add Cy5 SafO, and SafO
% 2020 July 4th
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function registration_KOMP(direct, out_direct1, exp_type, image_type, labels, name, bone_type, delimeter, GFP, bt, bn, section, gutta, mouse)

% direct='Z:\KOMP\Cbln3_E01V\01_Submitted\Layers\';
% out_direct1='Z:\KOMP\Cbln3_E01V\01_Submitted\Layers\Images\';
% exp_type={'Female';'Male'};
% image_type='jpg';
% labels={'green','green';'red','red'};   % first label, second label
% name={'AAA_E01_F';'AAA_E01_M'};         % gene name, exp_name, (female, male)
% bone_type={'V'};
% delimeter='\';
% GFP=0;              % if there is GFP, then 1. else 0
% gutta=0;              % no pre-registration (no gutta used)
% mouse=1;    % mouse = 1 (default)                
%             % rat   = 0


if nargin<14
    mouse=1;
else
    mouse=varargin{14};
end
bt=varargin{10};
bn=varargin{11};
section=varargin{12};
gutta=varargin{13};


% bone_type='F'             % Femur
% bone_type='V'             % Verterbae
% bone_type='C'             % Calvarial bone total
% bone_type='C_LR'          % Calvarial bone Left and Right

% name1='L88-1_C_P0_';
% name2='L88-1_C_P1_';
% name3='L88-1_C_P2_';


matlab_ver='7.10.0.499';
warning('off')


bone_type=varargin{7};
image_type=varargin{4};
exp_type=varargin{3};
delimiter=varargin{8};
labels=varargin{5};
GFP_flag=varargin{9};

rot_angle=0;

threshold_green_TRAP=[40, 40];      % threshold for beads of [green, TRAP];
threshold_red_ap=[70, 120];         % threshold for beads of [red, AP];

direct=varargin{1};
if strcmp(direct(end),delimiter)~=1
    direct=[direct,delimiter];
end

med_bubble=0;                                           % 0 : shift to match beads
                                                        % 1 : shift + median filter
out_direct1=varargin{2};
if strcmp(out_direct1(end),delimiter)~=1
    out_direct1=[out_direct1,delimiter];
end
aa=dir(out_direct1);
if isdir(out_direct1)==0
   mkdir(out_direct1)
end
                                                        
                                                        % 2 : shift + median filter + remove bubbles
if med_bubble==0
    mb_com='_shift.jpg';
elseif med_bubble==1
    mb_com='_med.jpg';
elseif med_bubble==2
    mb_com='_med_bubble.jpg';
end
                                                        
no_bone_type=length(varargin{6});
for i=1:no_bone_type
    eval(['temp=varargin{',num2str(6),'}{',num2str(i),'};'])
    eval(['p',num2str(i),'=temp;'])
end


for i=1:length(varargin{5})
    if strcmp(varargin{5}{i},'red')==1
        red_label=1;
    elseif strcmp(varargin{5}{i},'green')==1
        green_label=1;
    end
end        

% eval(['a=dir(''',direct,'info.mat'');'])
% 
% if isempty(a)==1                                        % read image information including im, ap, trap, H&E
    a=dir(direct);
%     if strcmp(bone_type, 'F')
        [info1, info2]=read_file_info_F2(a,varargin{6});
%     else strcmp(bone_type, 'V')
%         [info1, info2]=read_file_info(a,varargin{6});
%     end
    eval(['save ''',direct,'info'' info*'])
% else
%     eval(['load ''',direct,'info'''])
% end

sample_per_bone=0;
for i=bt
    eval(['temp=info',num2str(i),'.im;'])
%     temp1=sum(cell2mat(temp(:,2:end)),2);
    bone_no=temp(:,1);
    sample_per_bone(i)=max(max(sample_per_bone,size(temp,2)-1));
    eval(['bone_number',num2str(i),'=bone_no;'])
end

if mouse         
    bead_size=500;
    bead_margin=460;
else
    bead_size=205;
    bead_margin=195;
end
beads_special=1;            

if strcmp(bone_type, 'C_LR')==1
    LR1=['L', 'R'];
    direct=out_direct1;
else
    LR1=bone_type;
end

no_reg_info=0;

for i=bt%1:no_bone_type
    eval(['information=info',num2str(i),';'])
    eval(['bone_number=bone_number',num2str(i),';']);
    eval(['phr1=p',num2str(i),';'])

    if i==2
        st_bone=1;
    elseif i==1
        st_bone=1;
    end
    for j=bn%1:length(bone_number)
%         if (bt==1 & (bn==2 | bn==4 | bn==5)) | (bt==2 & (bn==2 | bn==3))
%             rot_angle=180;
%         else
            rot_angle=0;
%         end
        
        for s=section%1:sample_per_bone(i)

            for k=1:length(LR1)
                count_trap(j,s,i)=0;
                count_ap(j,s,i)=0;
                close all
                LR=LR1(k);
                [exp_type{i,:},' ', bone_number{j},',  section ',num2str(s)]

                direc=direct;
%                 nn=num2str(bone_number(j));
                nn=bone_number{j};

                c=[phr1,'L'];
    
                if strcmp(LR,'L')==1 | strcmp(LR,'R')==1             % Calvarial defect with Left and Right
                    c1=[c,nn,'_',num2str(s),'_2_',LR];               % DIC
                    c2=[c,nn,'_',num2str(s),'_3_',LR];               % TRAP
                    c3=[c,nn,'_',num2str(s),'_5_',LR];               % AP
                    c4=[c,nn,'_',num2str(s),'_1_',LR];               % Red label + red beads
                    c5=[c,nn,'_',num2str(s),'_0_',LR];               % host osteoblasts + green beads
                    c6=[c,nn,'_',num2str(s),'_4_',LR];               % Donor Cells
                    c7=[c,nn,'_',num2str(s),'_7_',LR];               % Green beads
                    c8=[c,nn,'_',num2str(s),'_6_',LR];               % DAPI
                    co=[c,nn,'_',num2str(s)];
                elseif strcmp(LR,'F')==1 | strcmp(LR,'V')==1 | strcmp(LR,'C')==1   % Femur, Vertebra, Calvarial defect without Left or Right

%                     c1=[c,nn,'_h',cell2mat(bone_type),'_M_s',num2str(s),'c2.jpg'];               % DIC
%                     c2=[c,nn,'_h',cell2mat(bone_type),'_T_s',num2str(s),'c2.jpg'];        % TRAP
%                     c3=[c,nn,'_h',cell2mat(bone_type),'_A_s',num2str(s),'c3.jpg'];              % AP
%                     if strcmp(labels{2,i},'red')                   % 2nd label
%                         c4=[c,nn,'_h',cell2mat(bone_type),'_M_s',num2str(s),'c4.jpg'];        % Red label + red beads
%                     elseif strcmp(labels{2,i},'green')
%                         c4=[c,nn,'_h',cell2mat(bone_type),'_M_s',num2str(s),'c3.jpg'];           % green label + green beads
%                     end
%                     if strcmp(labels{1,i},'red')                   % 1st label
%                         c5=[c,nn,'_h',cell2mat(bone_type),'_M_s',num2str(s),'c4.jpg'];        % Red label + red beads
%                     elseif strcmp(labels{1,i},'green')
%                         c5=[c,nn,'_h',cell2mat(bone_type),'_M_s',num2str(s),'c3.jpg'];           % green label + green beads
%                     end
%                     c6=[c,nn,'_',num2str(s),'_tr_CFP'];                 % Donor Cells
%                     c7=[c,nn,'_h',cell2mat(bone_type),'_T_s',num2str(s),'c1.jpg'];                % Green beads
%                     c8=[c,nn,'_h',cell2mat(bone_type),'_A_s',num2str(s),'c1.jpg'];                % DAPI
%                     co=[c,nn,'_h',cell2mat(bone_type),'_s',num2str(s)];
%                     c9=[c,nn,'_',num2str(s),'_dic_CFP'];                 % GFP
%                     c10=[c,nn,'_h',cell2mat(bone_type),'_T_s',num2str(s),'c3.jpg'];                % tomato from TRAP channel

%                     if strcmp(bone_type, 'V')
%                         nn=bone_number{j};
%                         c=[phr1];
%                         
%                         c1=[c,nn,'_h',cell2mat(bone_type),'_M_s',num2str(s),'c1.jpg'];               % DIC
%                         c2=[c,nn,'_h',cell2mat(bone_type),'_T_s',num2str(s),'c2.jpg'];        % TRAP
%                         c3=[c,nn,'_h',cell2mat(bone_type),'_A_s',num2str(s),'c3.jpg'];              % AP
%                         if strcmp(labels{2,i},'red')                   % 2nd label
%                             c4=[c,nn,'_h',cell2mat(bone_type),'_M_s',num2str(s),'c3.jpg'];        % Red label + red beads
%                         elseif strcmp(labels{2,i},'green')
%                             c4=[c,nn,'_h',cell2mat(bone_type),'_M_s',num2str(s),'c2.jpg'];           % green label + green beads
%                         end
%                         if strcmp(labels{1,i},'red')                   % 1st label
%                             c5=[c,nn,'_h',cell2mat(bone_type),'_M_s',num2str(s),'c3.jpg'];        % Red label + red beads
%                         elseif strcmp(labels{1,i},'green')
%                             c5=[c,nn,'_h',cell2mat(bone_type),'_M_s',num2str(s),'c2.jpg'];           % green label + green beads
%                         end
%                         c6=[c,nn,'_',num2str(s),'_tr_CFP'];                 % Donor Cells
%                         c7=[c,nn,'_h',cell2mat(bone_type),'_T_s',num2str(s),'c1.jpg'];                % Green beads
%                         c8=[c,nn,'_h',cell2mat(bone_type),'_A_s',num2str(s),'c1.jpg'];                % DAPI
%                         co=[c,nn,'_h',cell2mat(bone_type),'_s',num2str(s)];
%                         c9=[c,nn,'_',num2str(s),'_dic_CFP'];                 % GFP
%                         c10=[c,nn,'_h',cell2mat(bone_type),'_T_s',num2str(s),'c3.jpg'];                % tomato from TRAP channel
%                     elseif strcmp(bone_type, 'F')
                        nn=bone_number{j};
                        c=[phr1,'L'];
                        
                        c1=[c,num2str(s),'_M_s',nn,'c1.jpg'];               % DIC
                        c2=[c,num2str(s),'_T_s',nn,'c2.jpg'];        % TRAP 2020/6/26
%                         c2=[c,num2str(s),'_T_s',nn,'c1.jpg'];        % TRAP
%                        c3=[c,num2str(s),'_A_s',nn,'c3.jpg'];              % AP
                        c3=[c,num2str(s),'_A_s',nn,'c2.jpg'];              % AP
%                        if i==1
%                            c3=[c,num2str(s),'_A_s',nn,'c3.jpg'];              % AP
%                        elseif i==2
%                            c3=[c,num2str(s),'_A_s',nn,'c3.jpg'];              % AP
%                        end
%                        if i==1
%                            c3=[c,num2str(s),'_A_s',nn,'c2.jpg'];              % AP
%                        elseif i==2
%                            c3=[c,num2str(s),'_A_s',nn,'c2.jpg'];              % AP
%                        end
                        if strcmp(labels{2,i},'red')                   % 2nd label
                            c4=[c,num2str(s),'_M_s',nn,'c2.jpg'];        % Red label + red beads
                        elseif strcmp(labels{2,i},'green')
%                             c4=[c,num2str(s),'_M_s',nn,'c3.jpg'];           % green label + green beads
                            c4=[c,num2str(s),'_M_s',nn,'c4.jpg'];           % green label + green beads
                        end
                        if strcmp(labels{1,i},'red')                   % 1st label
                            c5=[c,num2str(s),'_M_s',nn,'c2.jpg'];        % Red label + red beads
                        elseif strcmp(labels{1,i},'green')
%                             c5=[c,num2str(s),'_M_s',nn,'c3.jpg'];           % green label + green beads
                            c5=[c,num2str(s),'_M_s',nn,'c4.jpg'];           % green label + green beads
                        end
                        c6=[c,nn,'_',num2str(s),'_tr_CFP'];                 % Donor Cells
%                         c7=[c,num2str(s),'_T_s',nn,'c2.jpg'];            % Green beads 2020/6/26
%                         c7=[c,num2str(s),'_T_s',nn,'c1.jpg'];                % TRAP 2020/6/26
                        c7=[c,num2str(s),'_T_s',nn,'c3.jpg'];                % TRAP beads 2022/7/26
%                        c8=[c,num2str(s),'_A_s',nn,'c1.jpg'];                % DAPI
                        c8=[c,num2str(s),'_A_s',nn,'c1.jpg'];                % DAPI
                        co=[c,num2str(s),'_s',nn];
                        c9=[c,num2str(s),'_',nn,'_dic_CFP'];                 % GFP
                        c10=[c,num2str(s),'_T_s',nn,'c1.jpg'];                % TRAP Cy5
                        c11=[c,num2str(s),'_B_s',nn,'c3.jpg'];                % tomato from TB channel
%                        c12=[c,num2str(s),'_B_s',nn,'c1.jpg'];                % TB
%                        c13=[c,num2str(s),'_B_s',nn,'c2.jpg'];                % green beads/TB
                        c12=[c,num2str(s),'_B_s',nn,'c1.jpg'];                % TB
                        c13=[c,num2str(s),'_B_s',nn,'c2.jpg'];                % green beads/TB
                        c14=[c,num2str(s),'_A_s',nn,'c3.jpg'];                % green beads/AP

                        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% registration_KOMP_8_v2channel
%
% Add Cy5 SafO, and SafO
% 2020 July 4th

                        c15=[c,num2str(s),'_C_s',nn,'c2.jpg'];                % green beads/Cy5 SafO
                        c16=[c,num2str(s),'_C_s',nn,'c1.jpg'];                % Cy5 SafO
                        c17=[c,num2str(s),'_S_s',nn,'c1.jpg'];                % SafO
                        c18=[c,num2str(s),'_S_s',nn,'c2.jpg'];                % green beads/SafO
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        c19=[c,num2str(s),'_M_s',nn,'c3.jpg'];                  % Mineral Cy5 (pink) : Mineral + background


%                     end
                end

                dis_v=10000; dis_h=10000; ratio_v=2; ratio_h=2; angle_between_green_trap=90; angle_between_green_trap1=90;

                if isfield(information,'im')==1
                    test1=cell2mat(information.im(find(~cellfun('isempty',(strfind(information.im(j:end,1),bone_number{j}))))+j-1,s+1));
                    if isempty(test1)
                        test1=0;
                    end
                else
                    test1=0;
                end
                if isfield(information,'tr')==1
                    test2=cell2mat(information.tr(find(~cellfun('isempty',(strfind(information.tr(j:end,1),bone_number{j}))))+j-1,s+1));
                    if isempty(test2)
                        test2=0;
                    end
                else
                    test2=0;
                end
                if isfield(information,'ap')==1
                    test3=cell2mat(information.ap(find(~cellfun('isempty',(strfind(information.ap(j:end,1),bone_number{j}))))+j-1,s+1));
                    if isempty(test3)
                        test3=0;
                    end
                else
                    test3=0;
                end
                
                if test1==0
                    reg_info(1,3)=10000;
                    no_reg_info=1;
                    continue
                elseif test1==1
%                     eval(['DIC1=imrotate(imread(''',direc,c1,''',''',image_type,'''),',num2str(rot_angle),');']);
%                     if length(size(DIC1))==3
%                         DIC1=DIC1(:,:,3);
%                     end
%                     eval(['red1=imrotate(imread(''',direc,c4,''',''',image_type,'''),',num2str(rot_angle),');']);
%                     red1=red1(:,:,1);
%                     eval(['Host1=imrotate(imread(''',direc,c5,''',''',image_type,'''),',num2str(rot_angle),');']);
%                     if length(size(Host1))==3
%                         Host1=Host1(:,:,2);
%                     end
% %                     eval(['M_Cy5=imrotate(imread(''',direc,c19,''',''',image_type,'''),',num2str(rot_angle),');']);
% %                     if length(size(M_Cy5))==3
% %                         M_Cy5=M_Cy5(:,:,1);
% %                     end
%                     
%                     common_background=uint8((Host1==128)&(red1==128)&(DIC1==128))*128;
%                     DIC1=imsubtract(DIC1,common_background);
%                     red1=imsubtract(red1,common_background);
%                     Host1=imsubtract(Host1,common_background);
% 
%                     if GFP_flag==1
%                         eval(['GFP1=imrotate(imread(''',direc,c9,''',''',image_type,'''),',num2str(rot_angle),');']);
%                         if length(size(GFP1))==3
%                             GFP1=GFP1(:,:,1);
%                         end
%                         GFP1=imsubtract(GFP1,common_background);
%                     else
%                         GFP1=uint8(false(size(DIC1)));
%                     end
                    
                    if test2==0  && test3==0                                                          % labels --> shift2
%                         if strcmp(LR,'C')==1 | strcmp(LR,'F')==1 | strcmp(LR,'V')==1
%                             eval(['imwrite(Host1,''',out_direct1,co,'_0_shift2.jpg'',''jpg'',''quality'',100);']);
%                             eval(['imwrite(red1,''',out_direct1,co,'_1_shift2.jpg'',''jpg'',''quality'',100);']);
%                             eval(['imwrite(DIC1,''',out_direct1,co,'_2_shift2.jpg'',''jpg'',''quality'',100);']);
%                             if GFP_flag==1
%                                 eval(['imwrite(GFP1,''',out_direct1,co,'_8_shift2.jpg'',''jpg'',''quality'',100);']);
%                             end
%                         elseif strcmp(LR,'L')==1 | strcmp(LR,'R')==1
%                             eval(['imwrite(Host1,''',out_direct1,co,'_0_',LR,'_shift2.jpg'',''jpg'',''quality'',100);']);
%                             eval(['imwrite(red1,''',out_direct1,co,'_1_',LR,'_shift2.jpg'',''jpg'',''quality'',100);']);
%                             eval(['imwrite(DIC1,''',out_direct1,co,'_2_',LR,'_shift2.jpg'',''jpg'',''quality'',100);']);
%                             if GFP_flag==1
%                                 eval(['imwrite(GFP1,''',out_direct1,co,'_8_',LR,'_shift2.jpg'',''jpg'',''quality'',100);']);
%                             end
%                         end
%                         aa=imadd(red1,immultiply(DIC1,0.2));
%                         aa(:,:,2)=imadd(Host1,immultiply(DIC1,0.2));
%                         if GFP_flag==1
%                             aa(:,:,3)=imadd(GFP1,immultiply(DIC1,0.2));
%                         else
%                             aa(:,:,3)=immultiply(DIC1,0.2);
%                         end
%                         clear Host1 red1 DIC1
%                         if strcmp(LR,'C')==1 | strcmp(LR,'F')==1 | strcmp(LR,'V')==1
%                             eval(['imwrite(aa,''',out_direct1,co,'_shift2.jpg'',''jpg'',''quality'',100);']);
%                             eval(['imwrite(aa,''',out_direct1,co,'_shift2_NoDAPI.jpg'',''jpg'',''quality'',100);']);
%                         elseif strcmp(LR,'L')==1 | strcmp(LR,'R')==1
%                             eval(['imwrite(aa,''',out_direct1,co,'_',LR,'_shift2.jpg'',''jpg'',''quality'',100);']);
%                             eval(['imwrite(aa,''',out_direct1,co,'_',LR,'_shift2_NoDAPI.jpg'',''jpg'',''quality'',100);']);
%                         end                            
%                         clear aa
%                         reg_info(1,3)=1000;
                    elseif test2==1  && test3==0        % TRAP 
% %                         ['start of TRAP at ', datestr(now)]
% %                         ref=[direc,c4];
% %                         im=[direc,c10];
%                         ref=[direc,c5];
%                         im=[direc,c7];
%                         [reg_info, no_bead, reg_error_code, image_after_rough_reg]=find_main_shift_6(ref, im, gutta, rot_angle, 0);
%                         reg_error_code=reg_error_code+10;
% %                       error_code  12: No beads or too much rotated or too big/small
% %                                       12.1 No beads
% %                                       12.2 too much rotated
% %                                       12.3 too big/small
% %                                   13: No beads
% %                                   11: No second Box after rough registration
%                         if no_bead==1
%                             reg_info(1,3)=101;
%                             eval(['imwrite(image_after_rough_reg,''',out_direct1,co,'_TRAP_error.jpg'',''jpg'',''quality'',100);']);
%                             break
%                         end
%                         reg_info(1,3)=1;                % TRAP only registration routine
%                         eval(['TRAP=imrotate(imread(''',direc,c2,''',''',image_type,'''),',num2str(rot_angle),');']);
%                         TRAP=TRAP(:,:,1);
% %                         eval(['T_Cy5=imrotate(imread(''',direc,c10,''',''',image_type,'''),',num2str(rot_angle),');']);
% %                         T_Cy5=T_Cy5(:,:,1);
%                         dis_v=reg_info(1,1);
%                         dis_h=reg_info(1,2);
%                         max_x=max(size(DIC1,2),size(TRAP,2));
%                         max_y=max(size(DIC1,1),size(TRAP,1));
%                         DIC1(max_y,max_x)=0;
%                         red1(max_y,max_x)=0;
%                         Host1(max_y,max_x)=0;
% %                         M_Cy5(max_y,max_x)=0;
%                         TRAP(max_y,max_x)=0;
% %                         T_Cy5(max_y,max_x)=0;
% %                         [DIC1, red1, Host1, M_Cy5, TRAP, T_Cy5]=shift_TRAP_AP_to_DIC(dis_v, dis_h, DIC1, red1, Host1, M_Cy5, TRAP, T_Cy5);
%                         [DIC1, red1, Host1, TRAP]=shift_TRAP_AP_to_DIC(dis_v, dis_h, DIC1, red1, Host1, TRAP);
%                         
%                         for reg=1:size(reg_info,1)-1
%                             TRAP=imrotate(TRAP,-reg_info(reg+1,1));
%                             max_x=max([size(DIC1,2),size(TRAP,2)]);
%                             max_y=max([size(DIC1,1),size(TRAP,1)]);
%                             DIC1(max_y,max_x)=0;
%                             red1(max_y,max_x)=0;
%                             Host1(max_y,max_x)=0;
% %                             M_Cy5(max_y,max_x)=0;
%                             TRAP(max_y,max_x)=0;
% %                             T_Cy5(max_y,max_x)=0;
%                             
%                             ratio_v=reg_info(reg+1,2);
%                             ratio_h=reg_info(reg+1,3);
%                             dis_v=reg_info(reg+1,4);
%                             dis_h=reg_info(reg+1,5);
%                             
%                             TRAP=imresize(TRAP,[round(ratio_v*size(TRAP,1)) round(ratio_h*size(TRAP,2))]);
%                             max_x=max(size(DIC1,2),size(TRAP,2));
%                             max_y=max(size(DIC1,1),size(TRAP,1));
%                             DIC1(max_y,max_x)=0;
%                             red1(max_y,max_x)=0;
%                             Host1(max_y,max_x)=0;
% %                             M_Cy5(max_y,max_x)=0;
%                             TRAP(max_y,max_x)=0;
% %                             T_Cy5(max_y,max_x)=0;
% %                             [DIC1, red1, Host1, M_Cy5, TRAP, T_Cy5]=shift_TRAP_AP_to_DIC(dis_v, dis_h, DIC1, red1, Host1, M_Cy5, TRAP, T_Cy5);
%                             [DIC1, red1, Host1, TRAP]=shift_TRAP_AP_to_DIC(dis_v, dis_h, DIC1, red1, Host1, TRAP);
%                         end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Remove Excessive outer empty spaces after TB is registered
% % 2020 Jun 30
% %
% %                         [DIC1, red1, Host1, M_Cy5, TRAP, T_Cy5]=remove_excessive_space(DIC1, red1, Host1, M_Cy5, TRAP, T_Cy5); 
%                         [DIC1, red1, Host1, TRAP]=remove_excessive_space(DIC1, red1, Host1, TRAP); 
% % 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                         
%                         eval(['imwrite(Host1,''',out_direct1,co,'_0_shift2.jpg'',''jpg'',''quality'',100);']);
%                         eval(['imwrite(red1,''',out_direct1,co,'_1_shift2.jpg'',''jpg'',''quality'',100);']);
%                         eval(['imwrite(DIC1,''',out_direct1,co,'_2_shift2.jpg'',''jpg'',''quality'',100);']);
%                         eval(['imwrite(TRAP,''',out_direct1,co,'_3_shift2.jpg'',''jpg'',''quality'',100);']);
% %                         eval(['imwrite(M_Cy5,''',out_direct1,co,'_11_shift2.jpg'',''jpg'',''quality'',100);']);
% %                         eval(['imwrite(T_Cy5,''',out_direct1,co,'_12_shift2.jpg'',''jpg'',''quality'',100);']);
% 
% 
%                         aa=imadd(TRAP, imadd(red1,immultiply(DIC1,0.3)));
%                         aa(:,:,2)=imadd(TRAP, imadd(Host1,immultiply(DIC1,0.3)));
%                         aa(:,:,3)=immultiply(DIC1,0.3);
%                         eval(['imwrite(aa,''',out_direct1,co,'_shift2.jpg'',''jpg'',''quality'',100);']);
%                         eval(['imwrite(aa,''',out_direct1,co,'_shift2_NoDAPI.jpg'',''jpg'',''quality'',100);']);
%                         ['end of TRAP at ', datestr(now)]

%                    elseif test2==0  && test3==1        % AP
% %                         ['start of AP at ', datestr(now)]
%                         ref=[direc,c4];
%                         im=[direc,c3];
%                         [reg_info, no_bead, reg_error_code, image_after_rough_reg]=find_main_shift_6(ref, im, gutta, rot_angle, 0);
%                         reg_error_code=reg_error_code+20;
% %                       error_code  22: No beads or too much rotated or too big/small
% %                                       22.1 No beads
% %                                       22.2 too much rotated
% %                                       22.3 too big/small
% %                                   23: No beads
% %                                   21: No second Box after rough registration
%                         if no_bead==1
%                             reg_info(1,3)=102;
%                             eval(['imwrite(image_after_rough_reg,''',out_direct1,co,'_AP_error.jpg'',''jpg'',''quality'',100);']);break
%                             break
%                         end
%                         reg_info(1,3)=2;                % AP only registration routine
%                         eval(['AP=imrotate(imread(''',direc,c3,''',''',image_type,'''),',num2str(rot_angle),');']);
%                         AP=AP(:,:,1);
%                         eval(['DAPI=imrotate(imread(''',direc,c8,''',''',image_type,'''),',num2str(rot_angle),');']);
%                         if length(size(DAPI))==3
%                             DAPI=DAPI(:,:,3);
%                         end
%                         dis_v=reg_info(1,1);
%                         dis_h=reg_info(1,2);
%                         max_x=max(size(DIC1,2),size(AP,2));
%                         max_y=max(size(DIC1,1),size(AP,1));
%                         DIC1(max_y,max_x)=0;
%                         red1(max_y,max_x)=0;
%                         Host1(max_y,max_x)=0;
%                         AP(max_y,max_x)=0;
%                         DAPI(max_y,max_x)=0;
%                         [DIC1, red1, Host1, AP, DAPI]=shift_TRAP_AP_to_DIC(dis_v, dis_h, DIC1, red1, Host1, AP, DAPI);
%                         
%                         for reg=1:size(reg_info,1)-1
%                             AP=imrotate(AP,-reg_info(reg+1,1));
%                             DAPI=imrotate(DAPI,-reg_info(reg+1,1));
%                             max_x=max([size(DIC1,2),size(AP,2)]);
%                             max_y=max([size(DIC1,1),size(AP,1)]);
%                             DIC1(max_y,max_x)=0;
%                             red1(max_y,max_x)=0;
%                             Host1(max_y,max_x)=0;
%                             AP(max_y,max_x)=0;
%                             DAPI(max_y,max_x)=0;
% 
%                             ratio_v=reg_info(reg+1,2);
%                             ratio_h=reg_info(reg+1,3);
%                             dis_v=reg_info(reg+1,4);
%                             dis_h=reg_info(reg+1,5);
%                             
%                             AP=imresize(AP,[round(ratio_v*size(AP,1)) round(ratio_h*size(AP,2))]);
%                             DAPI=imresize(DAPI,[round(ratio_v*size(DAPI,1)) round(ratio_h*size(DAPI,2))]);
%                             max_x=max(size(DIC1,2),size(AP,2));
%                             max_y=max(size(DIC1,1),size(AP,1));
%                             DIC1(max_y,max_x)=0;
%                             red1(max_y,max_x)=0;
%                             Host1(max_y,max_x)=0;
%                             AP(max_y,max_x)=0;
%                             DAPI(max_y,max_x)=0;
%                             [DIC1, red1, Host1, AP, DAPI]=shift_TRAP_AP_to_DIC(dis_v, dis_h, DIC1, red1, Host1, AP, DAPI);
%                         end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Remove Excessive outer empty spaces after TB is registered
% % 2020 Jun 30
% %
%                         [DIC1, red1, Host1, TRAP, AP, DAPI]=remove_excessive_space(DIC1, red1, Host1, TRAP, AP, DAPI); 
% % 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%                         eval(['imwrite(Host1,''',out_direct1,co,'_0_shift2.jpg'',''jpg'',''quality'',100);']);
%                         eval(['imwrite(red1,''',out_direct1,co,'_1_shift2.jpg'',''jpg'',''quality'',100);']);
%                         eval(['imwrite(DIC1,''',out_direct1,co,'_2_shift2.jpg'',''jpg'',''quality'',100);']);
%                         eval(['imwrite(AP,''',out_direct1,co,'_5_shift2.jpg'',''jpg'',''quality'',100);']);
%                         eval(['imwrite(DAPI,''',out_direct1,co,'_6_shift2.jpg'',''jpg'',''quality'',100);']);
%                         TRAP=uint8(false(size(DIC1))); 
%                         eval(['imwrite(TRAP,''',out_direct1,co,'_3_shift2.jpg'',''jpg'',''quality'',100);']);
%                         
%                         aa=imadd(AP, imadd(red1,immultiply(DIC1,0.3)));
%                         aa(:,:,2)=imadd(immultiply(AP,0.5),imadd(Host1,immultiply(DIC1,0.3)));
%                         aa(:,:,3)=imadd(DAPI, immultiply(DIC1,0.3));
%                         eval(['imwrite(aa,''',out_direct1,co,'_shift2.jpg'',''jpg'',''quality'',100);']);
%                         aa(:,:,3)=immultiply(DIC1,0.3);
%                         eval(['imwrite(aa,''',out_direct1,co,'_shift2_NoDAPI.jpg'',''jpg'',''quality'',100);']);
%                         ['end of AP at ', datestr(now)]

                    elseif  test2==1  && test3==1        % TRAP & AP
% %                         ['start of TRAP at ', datestr(now)]
% 
% %                         ref=[direc,c4];
% %                         im=[direc,c10];
%                         ref=[direc,c5];
%                         im=[direc,c7];
%                         [reg_info, no_bead, reg_error_code, image_after_rough_reg]=find_main_shift_6(ref, im, gutta, rot_angle, 0);
%                         reg_error_code=reg_error_code+10;
% %                       error_code  12: No beads or too much rotated or too big/small
% %                                       12.1 No beads
% %                                       12.2 too much rotated
% %                                       12.3 too big/small
% %                                   13: No beads
% %                                   11: No second Box after rough registration
%                         if no_bead==1
%                             reg_info(1,3)=103;
%                             eval(['imwrite(image_after_rough_reg,''',out_direct1,co,'_TRAP_error.jpg'',''jpg'',''quality'',100);']);
%                             break
%                         end
%                         reg_info(1,3)=3;                % TRAP and AP registration routine
%                         eval(['TRAP=imrotate(imread(''',direc,c2,''',''',image_type,'''),',num2str(rot_angle),');']);
%                         TRAP=TRAP(:,:,1);
%                         dis_v=reg_info(1,1);
%                         dis_h=reg_info(1,2);
%                         max_x=max(size(DIC1,2),size(TRAP,2));
%                         max_y=max(size(DIC1,1),size(TRAP,1));
%                         DIC1(max_y,max_x)=0;
%                         red1(max_y,max_x)=0;
%                         Host1(max_y,max_x)=0;
%                         TRAP(max_y,max_x)=0;
%                         [DIC1, red1, Host1, TRAP]=shift_TRAP_AP_to_DIC(dis_v, dis_h, DIC1, red1, Host1, TRAP);
%                         
%                         for reg=1:size(reg_info,1)-1
%                             TRAP=imrotate(TRAP,-reg_info(reg+1,1));
%                             max_x=max([size(DIC1,2),size(TRAP,2)]);
%                             max_y=max([size(DIC1,1),size(TRAP,1)]);
%                             DIC1(max_y,max_x)=0;
%                             red1(max_y,max_x)=0;
%                             Host1(max_y,max_x)=0;
%                             TRAP(max_y,max_x)=0;
% 
%                             ratio_v=reg_info(reg+1,2);
%                             ratio_h=reg_info(reg+1,3);
%                             dis_v=reg_info(reg+1,4);
%                             dis_h=reg_info(reg+1,5);
%                             
%                             TRAP=imresize(TRAP,[round(ratio_v*size(TRAP,1)) round(ratio_h*size(TRAP,2))]);
%                             max_x=max(size(DIC1,2),size(TRAP,2));
%                             max_y=max(size(DIC1,1),size(TRAP,1));
%                             DIC1(max_y,max_x)=0;
%                             red1(max_y,max_x)=0;
%                             Host1(max_y,max_x)=0;
%                             TRAP(max_y,max_x)=0;
%                             [DIC1, red1, Host1, TRAP]=shift_TRAP_AP_to_DIC(dis_v, dis_h, DIC1, red1, Host1, TRAP);
%                         end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Remove Excessive outer empty spaces after TB is registered
% % 2020 Jun 30
% %
%                         [DIC1, red1, Host1, TRAP]=remove_excessive_space(DIC1, red1, Host1, TRAP); 
% % 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%                         eval(['imwrite(Host1,''',out_direct1,co,'_0_shift1.jpg'',''jpg'',''quality'',100);']);
%                         eval(['imwrite(red1,''',out_direct1,co,'_1_shift1.jpg'',''jpg'',''quality'',100);']);
%                         eval(['imwrite(DIC1,''',out_direct1,co,'_2_shift1.jpg'',''jpg'',''quality'',100);']);
%                         eval(['imwrite(TRAP,''',out_direct1,co,'_3_shift1.jpg'',''jpg'',''quality'',100);']);
%                         aa=imadd(TRAP,imadd(red1,immultiply(DIC1,0.2)));
%                         aa(:,:,2)=imadd(TRAP,imadd(Host1,immultiply(DIC1,0.2)));
%                         aa(:,:,3)=immultiply(DIC1,0.2);
%                         eval(['imwrite(aa,''',out_direct1,co,'_shift1.jpg'',''jpg'',''quality'',100);']);
%                         ['end of TRAP at ', datestr(now)]


%   Start of AP

%                         ['start of AP at ', datestr(now)]

                        eval(['Host1=imread(''',out_direct1,co,'_0_shift1.jpg'',''jpg'');']);
                        eval(['red1=imread(''',out_direct1,co,'_1_shift1.jpg'',''jpg'');']);
%                         ref=red1;
                        ref=Host1;
                        eval(['DIC1=imread(''',out_direct1,co,'_2_shift1.jpg'',''jpg'');']);
                        eval(['TRAP=imread(''',out_direct1,co,'_3_shift1.jpg'',''jpg'');']);
                        
%                         im=[direc,c3];
                        im=[direc,c14];
                        [reg_info1, no_bead, reg_error_code, image_after_rough_reg]=find_main_shift_6(ref, im, gutta, rot_angle, 1);
                        reg_error_code=reg_error_code+20;
                        load([out_direct1,co,'_reg.mat'])
%                       error_code  22: No beads or too much rotated or too big/small
%                                       22.1 No beads
%                                       22.2 too much rotated
%                                       22.3 too big/small
%                                   23: No beads
%                                   21: No second Box after rough registration
                        if no_bead==1
                            reg_info1(1,3)=104;
%                             reg_info=[reg_info;reg_info1];
                            reg_info(end+1:end+size(reg_info1,1),1:size(reg_info1,2))=reg_info1;
                            eval(['imwrite(image_after_rough_reg,''',out_direct1,co,'_AP_error.jpg'',''jpg'',''quality'',100);']);
                            break
                        end
                        reg_info1(1,3)=4;
%                         reg_info=[reg_info;reg_info1];
%                         eval(['save ''', out_direct1,co,'_reg.mat'' reg_info;']);
%                         load([out_direct1,co,'_reg.mat'])
                        reg_info(end+1:end+size(reg_info1,1),1:size(reg_info1,2))=reg_info1;
                        eval(['AP=imrotate(imread(''',direc,c3,''',''',image_type,'''),',num2str(rot_angle),');']);
                        AP=AP(:,:,1);
                        eval(['DAPI=imrotate(imread(''',direc,c8,''',''',image_type,'''),',num2str(rot_angle),');']);
                        if length(size(DAPI))==3
                            DAPI=DAPI(:,:,3);
                        end
                        dis_v=reg_info1(1,1);
                        dis_h=reg_info1(1,2);
                        max_x=max(size(DIC1,2),size(AP,2));
                        max_y=max(size(DIC1,1),size(AP,1));
                        DIC1(max_y,max_x)=0;
                        red1(max_y,max_x)=0;
                        Host1(max_y,max_x)=0;
                        TRAP(max_y,max_x)=0;
                        AP(max_y,max_x)=0;
                        DAPI(max_y,max_x)=0;
%                         [DIC1, red1, Host1, TRAP, AP, DAPI]=shift_TRAP_AP_to_DIC(dis_v, dis_h, DIC1, red1, Host1, TRAP, AP, DAPI);
                        [DIC1, red1, Host1, TRAP, AP, DAPI]=shift_TRAP_AP_to_DIC(dis_v, dis_h, DIC1, red1, Host1, TRAP, AP, DAPI);
                        
                        for reg=1:size(reg_info1,1)-1
                            AP=imrotate(AP,-reg_info1(reg+1,1));
                            DAPI=imrotate(DAPI,-reg_info1(reg+1,1));
                            max_x=max([size(DIC1,2),size(AP,2)]);
                            max_y=max([size(DIC1,1),size(AP,1)]);
                            DIC1(max_y,max_x)=0;
                            red1(max_y,max_x)=0;
                            Host1(max_y,max_x)=0;
                            TRAP(max_y,max_x)=0;
                            AP(max_y,max_x)=0;
                            DAPI(max_y,max_x)=0;

                            ratio_v=reg_info1(reg+1,2);
                            ratio_h=reg_info1(reg+1,3);
                            dis_v=reg_info1(reg+1,4);
                            dis_h=reg_info1(reg+1,5);
                            
                            AP=imresize(AP,[round(ratio_v*size(AP,1)) round(ratio_h*size(AP,2))]);
                            DAPI=imresize(DAPI,[round(ratio_v*size(DAPI,1)) round(ratio_h*size(DAPI,2))]);
                            max_x=max(size(DIC1,2),size(AP,2));
                            max_y=max(size(DIC1,1),size(AP,1));
                            DIC1(max_y,max_x)=0;
                            red1(max_y,max_x)=0;
                            Host1(max_y,max_x)=0;
                            TRAP(max_y,max_x)=0;
                            AP(max_y,max_x)=0;
                            DAPI(max_y,max_x)=0;
%                             [DIC1, red1, Host1, TRAP, AP, DAPI]=shift_TRAP_AP_to_DIC(dis_v, dis_h, DIC1, red1, Host1, TRAP, AP, DAPI);
                            [DIC1, red1, Host1, TRAP, AP, DAPI]=shift_TRAP_AP_to_DIC(dis_v, dis_h, DIC1, red1, Host1, TRAP, AP, DAPI);
                        end
                        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remove Excessive outer empty spaces after TB is registered
% 2020 Jun 30
%
%                         [DIC1, red1, Host1, TRAP, AP, DAPI]=remove_excessive_space(DIC1, red1, Host1, TRAP, AP, DAPI); 
                        [DIC1, red1, Host1, TRAP, AP, DAPI]=remove_excessive_space(DIC1, red1, Host1, TRAP, AP, DAPI); 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                        
                        eval(['imwrite(Host1,''',out_direct1,co,'_0_shift2.jpg'',''jpg'',''quality'',100);']);
                        eval(['imwrite(red1,''',out_direct1,co,'_1_shift2.jpg'',''jpg'',''quality'',100);']);
                        eval(['imwrite(DIC1,''',out_direct1,co,'_2_shift2.jpg'',''jpg'',''quality'',100);']);
                        eval(['imwrite(TRAP,''',out_direct1,co,'_3_shift2.jpg'',''jpg'',''quality'',100);']);
                        eval(['imwrite(AP,''',out_direct1,co,'_5_shift2.jpg'',''jpg'',''quality'',100);']);
                        eval(['imwrite(DAPI,''',out_direct1,co,'_6_shift2.jpg'',''jpg'',''quality'',100);']);
                        aa=imadd(AP, imadd(TRAP, imadd(red1,immultiply(DIC1,0.3))));
                        aa(:,:,2)=imadd(TRAP, imadd(immultiply(AP,0.5),imadd(Host1,immultiply(DIC1,0.3))));
                        aa(:,:,3)=imadd(DAPI, immultiply(DIC1,0.3));
                        eval(['imwrite(aa,''',out_direct1,co,'_shift2.jpg'',''jpg'',''quality'',100);']);
                        aa(:,:,3)=immultiply(DIC1,0.3);
                        eval(['imwrite(aa,''',out_direct1,co,'_shift2_NoDAPI.jpg'',''jpg'',''quality'',100);']); clear aa
                        ['end of AP at ', datestr(now)]
                    end
% %   Start of TB
%                     if isfield(information,'tb')==1
%                         test4=cell2mat(information.tb(find(~cellfun('isempty',(strfind(information.tb(j:end,1),bone_number{j}))))+j-1,s+1));
%                         if isempty(test4)
%                             test4=0;
%                         end
%                     else
%                         test4=0;
%                     end
% 
%                     if test4==1
%                         ['start of TB at ', datestr(now)]
% %                             ref=red1;
% %                             im=[direc,c11];
%                         ref=Host1;
%                         im=[direc,c13];
%                         [reg_info, no_bead, reg_error_code, image_after_rough_reg]=find_main_shift_6(ref, im, gutta, rot_angle, 0);
%                         reg_error_code=reg_error_code+30;
% %                       error_code  32: No beads or too much rotated or too big/small
% %                                       32.1 No beads
% %                                       32.2 too much rotated
% %                                       32.3 too big/small
% %                                   33: No beads
% %                                   31: No second Box after rough registration
%                         if no_bead==1
%                             reg_info2(1,3)=105;
% %                                 reg_info=[reg_info;reg_info2];
%                             reg_info(end+1:end+size(reg_info2,1),1:size(reg_info2,2))=reg_info2;
%                             eval(['imwrite(image_after_rough_reg,''',out_direct1,co,'_TB_error.jpg'',''jpg'',''quality'',100);']);
%                             break
%                         end
%                         reg_info2(1,3)=5;
% %                             reg_info=[reg_info;reg_info2];
%                         reg_info(end+1:end+size(reg_info2,1),1:size(reg_info2,2))=reg_info2;
%                         eval(['TB_green=imrotate(imread(''',direc,c13,''',''',image_type,'''),',num2str(rot_angle),');']);
%                         TB_green=TB_green(:,:,2);
%                         eval(['TB=imrotate(imread(''',direc,c12,''',''',image_type,'''),',num2str(rot_angle),');']);
%                         TB1=TB(:,:,1); TB2=TB(:,:,2); TB3=TB(:,:,3);
% 
%                         dis_v=reg_info2(1,1);
%                         dis_h=reg_info2(1,2);
%                         max_x=max(size(DIC1,2),size(TB_green,2));
%                         max_y=max(size(DIC1,1),size(TB_green,1));
%                         DIC1(max_y,max_x)=0;
%                         red1(max_y,max_x)=0;
%                         Host1(max_y,max_x)=0;
%                         TRAP(max_y,max_x)=0;
%                         AP(max_y,max_x)=0;
%                         DAPI(max_y,max_x)=0;
%                         TB_green(max_y,max_x)=0;
%                         TB1(max_y,max_x)=0;
%                         TB2(max_y,max_x)=0;
%                         TB3(max_y,max_x)=0;
% 
% %                             for col=1:3
% %                                 TB1(:,:,col)=imrotate(TB(:,:,col),-reg_info1(reg+1,1));
% %                             end
%                         [DIC1, red1, Host1, TRAP, AP, DAPI, TB_green, TB1, TB2 ,TB3]=shift_TRAP_AP_to_DIC(dis_v, dis_h, DIC1, red1, Host1, TRAP, AP, DAPI, TB_green, TB1, TB2, TB3);
% 
%                         for reg=1:size(reg_info2,1)-1
%                             TB_green=imrotate(TB_green,-reg_info2(reg+1,1));
%                             TB1=imrotate(TB1,-reg_info2(reg+1,1));
%                             TB2=imrotate(TB2,-reg_info2(reg+1,1));
%                             TB3=imrotate(TB3,-reg_info2(reg+1,1));
%                             max_x=max([size(DIC1,2),size(TB_green,2)]);
%                             max_y=max([size(DIC1,1),size(TB_green,1)]);
%                             DIC1(max_y,max_x)=0;
%                             red1(max_y,max_x)=0;
%                             Host1(max_y,max_x)=0;
%                             TRAP(max_y,max_x)=0;
%                             AP(max_y,max_x)=0;
%                             DAPI(max_y,max_x)=0;
%                             TB_green(max_y,max_x)=0;
%                             TB1(max_y,max_x)=0;
%                             TB2(max_y,max_x)=0;
%                             TB3(max_y,max_x)=0;
% 
%                             ratio_v=reg_info2(reg+1,2);
%                             ratio_h=reg_info2(reg+1,3);
%                             dis_v=reg_info2(reg+1,4);
%                             dis_h=reg_info2(reg+1,5);
% 
%                             TB_green=imresize(TB_green,[round(ratio_v*size(TB_green,1)) round(ratio_h*size(TB_green,2))]);
%                             TB1=imresize(TB1,[round(ratio_v*size(TB1,1)) round(ratio_h*size(TB1,2))]);
%                             TB2=imresize(TB2,[round(ratio_v*size(TB2,1)) round(ratio_h*size(TB2,2))]);
%                             TB3=imresize(TB3,[round(ratio_v*size(TB3,1)) round(ratio_h*size(TB3,2))]);
%                             max_x=max(size(DIC1,2),size(TB_green,2));
%                             max_y=max(size(DIC1,1),size(TB_green,1));
%                             DIC1(max_y,max_x)=0;
%                             red1(max_y,max_x)=0;
%                             Host1(max_y,max_x)=0;
%                             TRAP(max_y,max_x)=0;
%                             AP(max_y,max_x)=0;
%                             DAPI(max_y,max_x)=0;
%                             TB_green(max_y,max_x)=0;
%                             TB1(max_y,max_x)=0;
%                             TB2(max_y,max_x)=0;
%                             TB3(max_y,max_x)=0;
% %                                 [DIC1, red1, Host1, TRAP, AP, DAPI]=shift_TRAP_AP_to_DIC(dis_v, dis_h, DIC1, red1, Host1, TRAP, AP, DAPI);
%                             [DIC1, red1, Host1, TRAP, AP, DAPI, TB_green, TB1, TB2 ,TB3]=shift_TRAP_AP_to_DIC(dis_v, dis_h, DIC1, red1, Host1, TRAP, AP, DAPI, TB_green, TB1, TB2, TB3);
%                         end
%                         TB=TB1; TB(:,:,2)=TB2; TB(:,:,3)=TB3;
%                     else
%                         TB1=uint8(false(size(DIC1)));
%                         TB2=TB1; TB3=TB1;
%                         TB=TB1; TB(:,:,2)=TB1; TB(:,:,3)=TB1;
%                         TB_green=TB1;
%                     end
%                     
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Remove Excessive outer empty spaces after TB is registered
% % 2020 Jun 30
% %
%                     [DIC1, red1, Host1, TRAP, AP, DAPI, TB_green, TB]=remove_excessive_space(DIC1, red1, Host1, TRAP, AP, DAPI, TB_green, TB); 
% % 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     eval(['imwrite(Host1,''',out_direct1,co,'_0_shift3.jpg'',''jpg'',''quality'',100);']);
%                     eval(['imwrite(red1,''',out_direct1,co,'_1_shift3.jpg'',''jpg'',''quality'',100);']);
%                     eval(['imwrite(DIC1,''',out_direct1,co,'_2_shift3.jpg'',''jpg'',''quality'',100);']);
%                     eval(['imwrite(TRAP,''',out_direct1,co,'_3_shift3.jpg'',''jpg'',''quality'',100);']);
%                     eval(['imwrite(AP,''',out_direct1,co,'_5_shift3.jpg'',''jpg'',''quality'',100);']);
%                     eval(['imwrite(DAPI,''',out_direct1,co,'_6_shift3.jpg'',''jpg'',''quality'',100);']);
%                     eval(['imwrite(TB,''',out_direct1,co,'_7_shift3.jpg'',''jpg'',''quality'',100);']);
%      
%                     TB1=TB(:,:,1); TB2=TB(:,:,2); TB3=TB(:,:,3);
%                     
%                     aa=imadd(imadd(AP, imadd(TRAP, imadd(red1,immultiply(DIC1,0.7)))), TB1);
%                     aa(:,:,2)=imadd(imadd(TRAP, imadd(immultiply(AP,0.5),imadd(Host1,immultiply(DIC1,0.7)))), TB2);
%                     aa(:,:,3)=imadd(imadd(DAPI, immultiply(DIC1,0.7)), TB3);
%                     eval(['imwrite(aa,''',out_direct1,co,'_shift3_TB.jpg'',''jpg'',''quality'',100);']);
% 
%                     aa=imadd(AP, imadd(TRAP, imadd(red1,immultiply(DIC1,0.3))));
%                     aa(:,:,2)=imadd(TRAP, imadd(immultiply(AP,0.5),imadd(Host1,immultiply(DIC1,0.3))));
%                     aa(:,:,3)=immultiply(DIC1,0.3);
%                     eval(['imwrite(aa,''',out_direct1,co,'_shift3_NoDAPI.jpg'',''jpg'',''quality'',100);']); clear aa
%                     ['end of TB at ', datestr(now)]
% 
%                 
% %   Start of Cy5 SafO
%                     if isfield(information,'cy')==1
%                         test5=cell2mat(information.cy(find(~cellfun('isempty',(strfind(information.cy(j:end,1),bone_number{j}))))+j-1,s+1));
%                         if isempty(test5)
%                             test5=0;
%                         end
%                     else
%                         test5=0;
%                     end
% 
%                     if test5==1
%                         ['start of Cy5 SafO at ', datestr(now)]
% %                             ref=red1;
% %                             im=[direc,c11];
%                         ref=Host1;
%                         im=[direc,c15];
%                         [reg_info, no_bead, reg_error_code, image_after_rough_reg]=find_main_shift_6(ref, im, gutta, rot_angle, 0);
%                         reg_error_code=reg_error_code+40;
% %                       error_code  42: No beads or too much rotated or too big/small
% %                                       42.1 No beads
% %                                       42.2 too much rotated
% %                                       42.3 too big/small
% %                                   43: No beads
% %                                   41: No second Box after rough registration
%                         if no_bead==1
%                             reg_info3(1,3)=106;
% %                                 reg_info=[reg_info;reg_info2];
%                             reg_info(end+1:end+size(reg_info3,1),1:size(reg_info3,2))=reg_info3;
%                             eval(['imwrite(image_after_rough_reg,''',out_direct1,co,'_Cy5_error.jpg'',''jpg'',''quality'',100);']);
%                             break
%                         end
%                         reg_info3(1,3)=6;
% %                             reg_info=[reg_info;reg_info2];
%                         reg_info(end+1:end+size(reg_info3,1),1:size(reg_info3,2))=reg_info3;
%                         eval(['Cy5_green=imrotate(imread(''',direc,c15,''',''',image_type,'''),',num2str(rot_angle),');']);
%                         Cy5_green=Cy5_green(:,:,2);
%                         eval(['Cy=imrotate(imread(''',direc,c16,''',''',image_type,'''),',num2str(rot_angle),');']);
%                         Cy1=Cy(:,:,1); Cy2=Cy(:,:,2); Cy3=Cy(:,:,3);
% 
%                         dis_v=reg_info3(1,1);
%                         dis_h=reg_info3(1,2);
%                         max_x=max(size(DIC1,2),size(Cy5_green,2));
%                         max_y=max(size(DIC1,1),size(Cy5_green,1));
%                         DIC1(max_y,max_x)=0;
%                         red1(max_y,max_x)=0;
%                         Host1(max_y,max_x)=0;
%                         TRAP(max_y,max_x)=0;
%                         AP(max_y,max_x)=0;
%                         DAPI(max_y,max_x)=0;
%                         TB1(max_y,max_x)=0;
%                         TB2(max_y,max_x)=0;
%                         TB3(max_y,max_x)=0;
%                         Cy5_green(max_y,max_x)=0;
%                         Cy1(max_y,max_x)=0;
%                         Cy2(max_y,max_x)=0;
%                         Cy3(max_y,max_x)=0;
% 
% %                             for col=1:3
% %                                 TB1(:,:,col)=imrotate(TB(:,:,col),-reg_info1(reg+1,1));
% %                             end
%                         [DIC1, red1, Host1, TRAP, AP, DAPI, TB1, TB2, TB3, Cy5_green, Cy1, Cy2 ,Cy3]=shift_TRAP_AP_to_DIC(dis_v, dis_h, DIC1, red1, Host1, TRAP, AP, DAPI, TB1, TB2, TB3, Cy5_green, Cy1, Cy2, Cy3);
% 
%                         for reg=1:size(reg_info3,1)-1
%                             Cy5_green=imrotate(Cy5_green,-reg_info3(reg+1,1));
%                             Cy1=imrotate(Cy1,-reg_info3(reg+1,1));
%                             Cy2=imrotate(Cy2,-reg_info3(reg+1,1));
%                             Cy3=imrotate(Cy3,-reg_info3(reg+1,1));
%                             max_x=max([size(DIC1,2),size(Cy5_green,2)]);
%                             max_y=max([size(DIC1,1),size(Cy5_green,1)]);
%                             DIC1(max_y,max_x)=0;
%                             red1(max_y,max_x)=0;
%                             Host1(max_y,max_x)=0;
%                             TRAP(max_y,max_x)=0;
%                             AP(max_y,max_x)=0;
%                             DAPI(max_y,max_x)=0;
%                             TB1(max_y,max_x)=0;
%                             TB2(max_y,max_x)=0;
%                             TB3(max_y,max_x)=0;
%                             Cy5_green(max_y,max_x)=0;
%                             Cy1(max_y,max_x)=0;
%                             Cy2(max_y,max_x)=0;
%                             Cy3(max_y,max_x)=0;
% 
%                             ratio_v=reg_info3(reg+1,2);
%                             ratio_h=reg_info3(reg+1,3);
%                             dis_v=reg_info3(reg+1,4);
%                             dis_h=reg_info3(reg+1,5);
% 
%                             Cy5_green=imresize(Cy5_green,[round(ratio_v*size(Cy5_green,1)) round(ratio_h*size(Cy5_green,2))]);
%                             Cy1=imresize(Cy1,[round(ratio_v*size(Cy1,1)) round(ratio_h*size(Cy1,2))]);
%                             Cy2=imresize(Cy2,[round(ratio_v*size(Cy2,1)) round(ratio_h*size(Cy2,2))]);
%                             Cy3=imresize(Cy3,[round(ratio_v*size(Cy3,1)) round(ratio_h*size(Cy3,2))]);
%                             max_x=max(size(DIC1,2),size(Cy5_green,2));
%                             max_y=max(size(DIC1,1),size(Cy5_green,1));
%                             DIC1(max_y,max_x)=0;
%                             red1(max_y,max_x)=0;
%                             Host1(max_y,max_x)=0;
%                             TRAP(max_y,max_x)=0;
%                             AP(max_y,max_x)=0;
%                             DAPI(max_y,max_x)=0;
%                             TB1(max_y,max_x)=0;
%                             TB2(max_y,max_x)=0;
%                             TB3(max_y,max_x)=0;
%                             Cy5_green(max_y,max_x)=0;
%                             Cy1(max_y,max_x)=0;
%                             Cy2(max_y,max_x)=0;
%                             Cy3(max_y,max_x)=0;
% 
%                             
% %                                 [DIC1, red1, Host1, TRAP, AP, DAPI]=shift_TRAP_AP_to_DIC(dis_v, dis_h, DIC1, red1, Host1, TRAP, AP, DAPI);
%                             [DIC1, red1, Host1, TRAP, AP, DAPI, TB1, TB2, TB3, Cy5_green, Cy1, Cy2 ,Cy3]=shift_TRAP_AP_to_DIC(dis_v, dis_h, DIC1, red1, Host1, TRAP, AP, DAPI, TB1, TB2, TB3, Cy5_green, Cy1, Cy2, Cy3);
%                         end
%                         Cy=Cy1; Cy(:,:,2)=Cy2; Cy(:,:,3)=Cy3;
%                         TB=TB1; TB(:,:,2)=TB2; TB(:,:,3)=TB3;
%                     else
%                         Cy1=uint8(false(size(DIC1)));
%                         Cy2=Cy1; Cy3=Cy1;
%                         Cy=Cy1; Cy(:,:,2)=Cy1; Cy(:,:,3)=Cy1;
%                     end
%                  
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Remove Excessive outer empty spaces after TB is registered
% % 2020 Jun 30
% %
%                     [DIC1, red1, Host1, TRAP, AP, DAPI, TB, Cy]=remove_excessive_space(DIC1, red1, Host1, TRAP, AP, DAPI, TB, Cy); 
% % 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     eval(['imwrite(Host1,''',out_direct1,co,'_0_shift3.jpg'',''jpg'',''quality'',100);']);
%                     eval(['imwrite(red1,''',out_direct1,co,'_1_shift3.jpg'',''jpg'',''quality'',100);']);
%                     eval(['imwrite(DIC1,''',out_direct1,co,'_2_shift3.jpg'',''jpg'',''quality'',100);']);
%                     eval(['imwrite(TRAP,''',out_direct1,co,'_3_shift3.jpg'',''jpg'',''quality'',100);']);
%                     eval(['imwrite(AP,''',out_direct1,co,'_5_shift3.jpg'',''jpg'',''quality'',100);']);
%                     eval(['imwrite(DAPI,''',out_direct1,co,'_6_shift3.jpg'',''jpg'',''quality'',100);']);
%                     eval(['imwrite(TB,''',out_direct1,co,'_7_shift3.jpg'',''jpg'',''quality'',100);']);
%                     eval(['imwrite(Cy,''',out_direct1,co,'_8_shift3.jpg'',''jpg'',''quality'',100);']);
%                     
%                     Cy1=Cy(:,:,1); Cy2=Cy(:,:,2); Cy3=Cy(:,:,3);
%                     TB1=TB(:,:,1); TB2=TB(:,:,2); TB3=TB(:,:,3);
%                     
%                     aa=imadd(imadd(AP, imadd(TRAP, imadd(red1,immultiply(DIC1,0.7)))), TB1);
%                     aa(:,:,2)=imadd(imadd(TRAP, imadd(immultiply(AP,0.5),imadd(Host1,immultiply(DIC1,0.7)))), TB2);
%                     aa(:,:,3)=imadd(imadd(DAPI, immultiply(DIC1,0.7)), TB3);
%                     eval(['imwrite(aa,''',out_direct1,co,'_shift3_TB.jpg'',''jpg'',''quality'',100);']);
%                     
%                     aa=imadd(imadd(AP, imadd(TRAP, imadd(red1,immultiply(DIC1,0.7)))), Cy1);
%                     aa(:,:,2)=imadd(imadd(TRAP, imadd(immultiply(AP,0.5),imadd(Host1,immultiply(DIC1,0.7)))), Cy2);
%                     aa(:,:,3)=imadd(imadd(DAPI, immultiply(DIC1,0.7)), Cy3);
%                     eval(['imwrite(aa,''',out_direct1,co,'_shift3_Cy.jpg'',''jpg'',''quality'',100);']);
% 
%                     aa=imadd(AP, imadd(TRAP, imadd(red1,immultiply(DIC1,0.3))));
%                     aa(:,:,2)=imadd(TRAP, imadd(immultiply(AP,0.5),imadd(Host1,immultiply(DIC1,0.3))));
%                     aa(:,:,3)=immultiply(DIC1,0.3);
%                     eval(['imwrite(aa,''',out_direct1,co,'_shift3_NoDAPI.jpg'',''jpg'',''quality'',100);']); clear aa
%                     ['end of Cy5 SafO at ', datestr(now)]                
%                     
%                     
%                     
% %   Start of SafO
%                     if isfield(information,'so')==1
%                         test6=cell2mat(information.so(find(~cellfun('isempty',(strfind(information.so(j:end,1),bone_number{j}))))+j-1,s+1));
%                         if isempty(test6)
%                             test6=0;
%                         end
%                     else
%                         test6=0;
%                     end
% 
%                     if test6==1
%                         ['start of SafO at ', datestr(now)]
% %                             ref=red1;
% %                             im=[direc,c11];
%                         ref=Host1;
%                         im=[direc,c18];
%                         [reg_info, no_bead, reg_error_code, image_after_rough_reg]=find_main_shift_6(ref, im, gutta, rot_angle, 0);
%                         reg_error_code=reg_error_code+50;
% %                       error_code  52: No beads or too much rotated or too big/small
% %                                       52.1 No beads
% %                                       52.2 too much rotated
% %                                       52.3 too big/small
% %                                   53: No beads
% %                                   51: No second Box after rough registration
%                         if no_bead==1
%                             reg_info4(1,3)=107;
% %                                 reg_info=[reg_info;reg_info2];
%                             reg_info(end+1:end+size(reg_info4,1),1:size(reg_info4,2))=reg_info4;
%                             eval(['imwrite(image_after_rough_reg,''',out_direct1,co,'_SO_error.jpg'',''jpg'',''quality'',100);']);
%                             break
%                         end
%                         reg_info4(1,3)=7;
% %                             reg_info=[reg_info;reg_info2];
%                         reg_info(end+1:end+size(reg_info4,1),1:size(reg_info4,2))=reg_info4;
%                         eval(['So_green=imrotate(imread(''',direc,c18,''',''',image_type,'''),',num2str(rot_angle),');']);
%                         So_green=So_green(:,:,2);
%                         eval(['So=imrotate(imread(''',direc,c17,''',''',image_type,'''),',num2str(rot_angle),');']);
%                         So1=So(:,:,1); So2=So(:,:,2); So3=So(:,:,3);
% 
%                         dis_v=reg_info4(1,1);
%                         dis_h=reg_info4(1,2);
%                         max_x=max(size(DIC1,2),size(So_green,2));
%                         max_y=max(size(DIC1,1),size(So_green,1));
%                         DIC1(max_y,max_x)=0;
%                         red1(max_y,max_x)=0;
%                         Host1(max_y,max_x)=0;
%                         TRAP(max_y,max_x)=0;
%                         AP(max_y,max_x)=0;
%                         DAPI(max_y,max_x)=0;
%                         TB1(max_y,max_x)=0;
%                         TB2(max_y,max_x)=0;
%                         TB3(max_y,max_x)=0;
%                         Cy1(max_y,max_x)=0;
%                         Cy2(max_y,max_x)=0;
%                         Cy3(max_y,max_x)=0;
%                         So_green(max_y,max_x)=0;
%                         So1(max_y,max_x)=0;
%                         So2(max_y,max_x)=0;
%                         So3(max_y,max_x)=0;
% 
% %                             for col=1:3
% %                                 TB1(:,:,col)=imrotate(TB(:,:,col),-reg_info1(reg+1,1));
% %                             end
%                         [DIC1, red1, Host1, TRAP, AP, DAPI, TB1, TB2, TB3, Cy1, Cy2, Cy3, So_green, So1, So2 ,So3]=shift_TRAP_AP_to_DIC(dis_v, dis_h, DIC1, red1, Host1, TRAP, AP, DAPI, TB1, TB2, TB3, Cy1, Cy2, Cy3, So_green, So1, So2, So3);
% 
%                         for reg=1:size(reg_info4,1)-1
%                             So_green=imrotate(So_green,-reg_info4(reg+1,1));
%                             So1=imrotate(So1,-reg_info4(reg+1,1));
%                             So2=imrotate(So2,-reg_info4(reg+1,1));
%                             So3=imrotate(So3,-reg_info4(reg+1,1));
%                             max_x=max([size(DIC1,2),size(So_green,2)]);
%                             max_y=max([size(DIC1,1),size(So_green,1)]);
%                             DIC1(max_y,max_x)=0;
%                             red1(max_y,max_x)=0;
%                             Host1(max_y,max_x)=0;
%                             TRAP(max_y,max_x)=0;
%                             AP(max_y,max_x)=0;
%                             DAPI(max_y,max_x)=0;
%                             TB1(max_y,max_x)=0;
%                             TB2(max_y,max_x)=0;
%                             TB3(max_y,max_x)=0;
%                             Cy1(max_y,max_x)=0;
%                             Cy2(max_y,max_x)=0;
%                             Cy3(max_y,max_x)=0;
%                             So5_green(max_y,max_x)=0;
%                             So1(max_y,max_x)=0;
%                             So2(max_y,max_x)=0;
%                             So3(max_y,max_x)=0;
% 
%                             ratio_v=reg_info4(reg+1,2);
%                             ratio_h=reg_info4(reg+1,3);
%                             dis_v=reg_info4(reg+1,4);
%                             dis_h=reg_info4(reg+1,5);
% 
%                             So_green=imresize(So_green,[round(ratio_v*size(So_green,1)) round(ratio_h*size(So_green,2))]);
%                             So1=imresize(So1,[round(ratio_v*size(So1,1)) round(ratio_h*size(So1,2))]);
%                             So2=imresize(So2,[round(ratio_v*size(So2,1)) round(ratio_h*size(So2,2))]);
%                             So3=imresize(So3,[round(ratio_v*size(So3,1)) round(ratio_h*size(So3,2))]);
%                             max_x=max(size(DIC1,2),size(So_green,2));
%                             max_y=max(size(DIC1,1),size(So_green,1));
%                             DIC1(max_y,max_x)=0;
%                             red1(max_y,max_x)=0;
%                             Host1(max_y,max_x)=0;
%                             TRAP(max_y,max_x)=0;
%                             AP(max_y,max_x)=0;
%                             DAPI(max_y,max_x)=0;
%                             TB1(max_y,max_x)=0;
%                             TB2(max_y,max_x)=0;
%                             TB3(max_y,max_x)=0;
%                             Cy1(max_y,max_x)=0;
%                             Cy2(max_y,max_x)=0;
%                             Cy3(max_y,max_x)=0;
%                             So_green(max_y,max_x)=0;
%                             So1(max_y,max_x)=0;
%                             So2(max_y,max_x)=0;
%                             So3(max_y,max_x)=0;
% %                                 [DIC1, red1, Host1, TRAP, AP, DAPI]=shift_TRAP_AP_to_DIC(dis_v, dis_h, DIC1, red1, Host1, TRAP, AP, DAPI);
%                             [DIC1, red1, Host1, TRAP, AP, DAPI, TB1, TB2, TB3, Cy1, Cy2, Cy3, So_green, So1, So2 ,So3]=shift_TRAP_AP_to_DIC(dis_v, dis_h, DIC1, red1, Host1, TRAP, AP, DAPI, TB1, TB2, TB3, Cy1, Cy2, Cy3, So_green, So1, So2, So3);
%                         end
%                         Cy=Cy1; Cy(:,:,2)=Cy2; Cy(:,:,3)=Cy3;
%                         TB=TB1; TB(:,:,2)=TB2; TB(:,:,3)=TB3;
%                         So=So1; So(:,:,2)=So2; So(:,:,3)=So3;
%                     else
%                         So1=uint8(false(size(DIC1)));
%                         So2=So1; So3=So1;
%                         So=So1; So(:,:,2)=So1; So(:,:,3)=So1;
%                     end
%                     
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Remove Excessive outer empty spaces after TB is registered
% % 2020 Jun 30
% %
%                     [DIC1, red1, Host1, TRAP, AP, DAPI, TB, Cy, So]=remove_excessive_space(DIC1, red1, Host1, TRAP, AP, DAPI, TB, Cy, So); 
% % 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     eval(['imwrite(Host1,''',out_direct1,co,'_0_shift3.jpg'',''jpg'',''quality'',100);']);
%                     eval(['imwrite(red1,''',out_direct1,co,'_1_shift3.jpg'',''jpg'',''quality'',100);']);
%                     eval(['imwrite(DIC1,''',out_direct1,co,'_2_shift3.jpg'',''jpg'',''quality'',100);']);
%                     eval(['imwrite(TRAP,''',out_direct1,co,'_3_shift3.jpg'',''jpg'',''quality'',100);']);
%                     eval(['imwrite(AP,''',out_direct1,co,'_5_shift3.jpg'',''jpg'',''quality'',100);']);
%                     eval(['imwrite(DAPI,''',out_direct1,co,'_6_shift3.jpg'',''jpg'',''quality'',100);']);
%                     eval(['imwrite(TB,''',out_direct1,co,'_7_shift3.jpg'',''jpg'',''quality'',100);']);
%                     eval(['imwrite(Cy,''',out_direct1,co,'_8_shift3.jpg'',''jpg'',''quality'',100);']);
%                     eval(['imwrite(So,''',out_direct1,co,'_9_shift3.jpg'',''jpg'',''quality'',100);']);
%                     
%                     So1=So(:,:,1); So2=So(:,:,2); So3=So(:,:,3);
%                     Cy1=Cy(:,:,1); Cy2=Cy(:,:,2); Cy3=Cy(:,:,3);
%                     TB1=TB(:,:,1); TB2=TB(:,:,2); TB3=TB(:,:,3);
%                     
%                     aa=imadd(imadd(AP, imadd(TRAP, imadd(red1,immultiply(DIC1,0.7)))), TB1);
%                     aa(:,:,2)=imadd(imadd(TRAP, imadd(immultiply(AP,0.5),imadd(Host1,immultiply(DIC1,0.7)))), TB2);
%                     aa(:,:,3)=imadd(imadd(DAPI, immultiply(DIC1,0.7)), TB3);
%                     eval(['imwrite(aa,''',out_direct1,co,'_shift3_TB.jpg'',''jpg'',''quality'',100);']);
%                     
%                     aa=imadd(imadd(AP, imadd(TRAP, imadd(red1,immultiply(DIC1,0.7)))), Cy1);
%                     aa(:,:,2)=imadd(imadd(TRAP, imadd(immultiply(AP,0.5),imadd(Host1,immultiply(DIC1,0.7)))), Cy2);
%                     aa(:,:,3)=imadd(imadd(DAPI, immultiply(DIC1,0.7)), Cy3);
%                     eval(['imwrite(aa,''',out_direct1,co,'_shift3_Cy.jpg'',''jpg'',''quality'',100);']);
%                     
%                     aa=imadd(imadd(AP, imadd(TRAP, imadd(red1,immultiply(DIC1,0.7)))), So1);
%                     aa(:,:,2)=imadd(imadd(TRAP, imadd(immultiply(AP,0.5),imadd(Host1,immultiply(DIC1,0.7)))), So2);
%                     aa(:,:,3)=imadd(imadd(DAPI, immultiply(DIC1,0.7)), So3);
%                     eval(['imwrite(aa,''',out_direct1,co,'_shift3_So.jpg'',''jpg'',''quality'',100);']);
% 
%                     aa=imadd(AP, imadd(TRAP, imadd(red1,immultiply(DIC1,0.3))));
%                     aa(:,:,2)=imadd(TRAP, imadd(immultiply(AP,0.5),imadd(Host1,immultiply(DIC1,0.3))));
%                     aa(:,:,3)=immultiply(DIC1,0.3);
%                     eval(['imwrite(aa,''',out_direct1,co,'_shift3_NoDAPI.jpg'',''jpg'',''quality'',100);']); clear aa
%                     ['end of SafO at ', datestr(now)]                      
%                     
%                     
                end
% 
                ['end of [',num2str(i),', ',bone_number{j},', ',num2str(s),', ',num2str(k),'] at ', datestr(now)]
%                 clear aa DAPI AP
%                 close all
           end     % for k
        end         % while s
    end             % for j
end                 % for no_bone_type
% if test1==1 & (test2==1 | test3==1)
%     eval(['save ''', out_direct1,co,'_reg.mat'' reg_info;']);
% end
% if ~no_reg_info
    eval(['save ''', out_direct1,co,'_reg.mat'' reg_info;']);
% end
return

function varargout=shift_TRAP_AP_to_DIC(varargin) 
% function shift_both_images(dis_v, dis_h, DIC1,
dis_v=varargin{1};
dis_h=varargin{2};
if nargin==6        %  (DIC, red, Green)-->(TRAP)
    no_ref=3;
    no_test=1;
elseif nargin==8    % (DIC, red, green, TRAP)-->(AP & DAPI)
    no_ref=4;
    no_test=2;

% elseif nargin==10    % (TRAP & T_Cy5) --> (AP & DAPI)
%     no_ref=6;
%     no_test=2;

% elseif nargin==12    % (TRAP & T_Cy5) --> (AP & DAPI), (TB_red, TB)
%     no_ref=7;
%     no_test=4;
% elseif nargin==15    % (TRAP & T_Cy5) --> (AP & DAPI), (TB_red, TB)
%     no_ref=10;
%     no_test=4;
% elseif nargin==18    % (TRAP & T_Cy5) --> (AP & DAPI), (TB_red, TB)
%     no_ref=13;
%     no_test=4;

elseif nargin==12    % (TRAP & T_Cy5), (AP & DAPI) -->  (TB_red, TB)
    no_ref=7;
    no_test=4;
elseif nargin==15    % (TRAP & T_Cy5), (AP & DAPI), (TB_red, TB), --> Cy5
    no_ref=10;
    no_test=4;
elseif nargin==18    % (TRAP & T_Cy5), (AP & DAPI), (TB_red, TB), Cy5 --> SO
    no_ref=13;
    no_test=4;
end
if dis_v>0
    for i=1:no_ref
        temp=varargin{2+i};
        temp(1:end+dis_v,end)=0;
        varargout{i}=temp;
    end
    for i=1:no_test
        temp1=uint8(false(size(temp)));
        temp1(dis_v+1:end,:)=varargin{2+no_ref+i};
        varargout{no_ref+i}=temp1;
    end
else
    for i=1:no_test
        temp=varargin{2+no_ref+i};
        temp(1:end-dis_v,end)=0;
        varargout{no_ref+i}=temp;
    end
    for i=1:no_ref
        temp1=uint8(false(size(temp)));
        temp1(-dis_v+1:end,:)=varargin{2+i};
        varargout{i}=temp1;
    end
end
if dis_h>0
    for i=1:no_ref
        temp=varargout{i};
        temp(end,1:end+dis_h)=0;
        varargout{i}=temp;
    end
    for i=1:no_test
        temp1=uint8(false(size(temp)));
        temp1(:,dis_h+1:end)=varargout{no_ref+i};
        varargout{no_ref+i}=temp1;
    end
else
    for i=1:no_test
        temp=varargout{no_ref+i};
        temp(end,1:end-dis_h)=0;
        varargout{no_ref+i}=temp;
    end
    for i=1:no_ref
        temp1=uint8(false(size(temp)));
        temp1(:,-dis_h+1:end)=varargout{i};
        varargout{i}=temp1;
    end
end

function varargout=shift_TB_to_DIC(varargin) 
% function shift_both_images(dis_v, dis_h, DIC1, TB) 
dis_v=varargin{1};
dis_h=varargin{2};

no_ref=1;
no_test=1;

if dis_v>0
    for i=1:no_ref
        temp=varargin{2+i};
        temp(1:end+dis_v,end)=0;
        varargout{i}=temp;
    end
    for i=1:no_test
        temp1=uint8(false([size(temp),3]));
        temp1(dis_v+1:end,:,:)=varargin{2+no_ref+i};
        varargout{no_ref+i}=temp1;
    end
else
    for i=1:no_test
        temp=varargin{2+no_ref+i};
        temp(1:end-dis_v,end,3)=0;
        varargout{no_ref+i}=temp;
    end
    for i=1:no_ref
        [x y z]=size(temp);
        temp1=uint8(false(x, y));
        temp1(-dis_v+1:end,:)=varargin{2+i};
        varargout{i}=temp1;
    end
end
if dis_h>0
    for i=1:no_ref
        temp=varargout{i};
        temp(end,1:end+dis_h)=0;
        varargout{i}=temp;
    end
    for i=1:no_test
        temp1=uint8(false([size(temp),3]));
        temp1(:,dis_h+1:end,:)=varargout{no_ref+i};
        varargout{no_ref+i}=temp1;
    end
else
    for i=1:no_test
        temp=varargout{no_ref+i};
        temp(end,1:end-dis_h,3)=0;
        varargout{no_ref+i}=temp;
    end
    for i=1:no_ref
        [x y z]=size(temp);
        temp1=uint8(false(x,y));
        temp1(:,-dis_h+1:end)=varargout{i};
        varargout{i}=temp1;
    end
end

function [varargout]=remove_excessive_space(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   remove the excessive empty outer space after TB is registered
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function [DIC1, red1, Host1, TRAP1, AP1, DAPI1, TB_green1, TB1]=remove_excessive_space(DIC, red, Host, TRAP, AP, DAPI, TB_green, TB)

% a=imadd(DIC, imadd(red, imadd(Host, imadd(TRAP, imadd(AP, imadd(DAPI, imadd(TB_green, imadd(TB(:,:,1),imadd(TB(:,:,2),TB(:,:,3))))))))));

n=nargin;
list3D=[];
for i=1:n
    aa=varargin{i};
    if length(size(aa))==2
%         eval(['a',num2str(i),'=aa;'])
        if i==1
            a=aa;
        else
            a=imadd(a,aa);
        end
    elseif length(size(aa))==3
        list3D=[list3D,i];
%         eval(['a',num2str(i),'=aa;'])
        for j=1:3
            a=imadd(a,aa(:,:,j));
        end
    end
end
% a1=sum(a,3);
ah=sum(a,2);
av=sum(a,1);
sx=find(ah,1,'first');
ex=find(ah,1,'last');
sy=find(av,1,'first');
ey=find(av,1,'last');

for i=1:n
    temp=varargin{i};
    if ~isempty(find(list3D==i, 1))
        varargout{i}=temp(sx:ex,sy:ey,:);
    else
        varargout{i}=temp(sx:ex,sy:ey);
    end
end

function [reg_info, no_bead, error_code, aa]=find_main_shift_6(ref, im, gutta, rot_angle, ap)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Find the discrepancies of the ref and im using goota
%
%   ref             : reference image (red label image)
%   im              : test image (tomato channel of AP/TRAP image)
%   gutta           : if gutta is used gutta=1 (pre_regist with gutta), else gutta=0
%
%   reg_info        : registration information
%                   : 1st row --> [dis_v, dis_h]
%                   : 2nd row --> [angle, ratio_v, ratio_h, dis_v, dis_h] of AP/TRAP
%                   : if there are more than 2 rows, 3rd row and later rows are the same as 2nd row
%
%   February 4, 2013
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
error_code=0;
init_thre=100;
init_thre1=50;
% init_thre=150;      % July 15, 2014 for Lisa's scanning

if size(ref,1)==1
% uiopen('Z:\KOMP\AAA_E01F\01_Submitted\Layers\AAA_E01_F03_hF_A_s1c3_ORG.jpg',1)
    a=imread(ref);
    a=imrotate(a, rot_angle);
    if ap==1
        a=a(:,:,1);
    else
        a=a(:,:,2);
    end
else
    a=ref;
end
% a=a>graythresh(a)*max(a(:));
% a=a>200;

% a1=(find_circles(a>init_thre, 100, 500)>0) | (find_circles(a>init_thre1, 100, 1000)>0);		% red beads (in Mineral) for AP
% if ap==0 | ap==2 								% green beads (in Mineral) for TRAP or TB
%     a1=a1|(find_circles(a>30, 50, 500)>0) | (find_circles(a>150, 100, 500)>0);
% end
a1=(find_circles(a>init_thre, 100, 500)>0) | (find_circles(a>init_thre1, 100, 1000)>0) | (find_circles(a>20, 100, 5000)>0) | (find_circles(imclose(a>250,strel('disk',5)), 100, 5000)>0);		% red beads (in Mineral) for AP
if ap==0 %| ap==2 								% green beads (in Mineral) for TRAP or TB
    a1=a1|(find_circles(a>30, 50, 500)>0) | (find_circles(a>150, 100, 500)>0);
elseif ap==2
%     a1=(find_circles(a>init_thre, 100, 500)>0) | (find_circles(a>init_thre1, 100, 500)>0) | (find_circles(a>30, 50, 500)>0) | (find_circles(a>150, 100, 500)>0);
    a2=a-medfilt2(a,[51,51]);
%     a1=(find_circles(a>init_thre, 100, 500)>0) | (find_circles(a>init_thre1, 100, 500)>0) | (find_circles(a>30, 50, 500)>0) | (find_circles(a>150, 100, 500)>0) | (find_circles(a>20, 100, 500)>0);
    a1=(find_circles(a2>5, 100, 500)>0) | (find_circles(a2>10, 100, 500)>0);
end
% a1=(find_circles(a>init_thre1, 100, 1000)>0);

if size(im,1)==1
    d=imread(im);
    if ap==1
        d=d(:,:,1);
    else
        d=d(:,:,2);
    end
    d=imrotate(d, rot_angle);
else
    d=im;
end

% d=d>graythresh(d)*max(d(:));
% d=d>200;

if ap==2
%     d1=find_circles(d>graythresh(d)*max(d(:)), 30, 500)>0;
%    d1=d>graythresh(d)*max(d(:));						% green beads in TB
    d2=d-medfilt2(d,[401,401]);
    d1=(find_circles(d2>50, 50, 300)>0) | (find_circles(d2>40, 50, 300)>0) | (find_circles(d2>30, 50, 300)>0) | (find_circles(d2>20, 50, 300)>0) | (find_circles(d2>15, 50, 300)>0) | (find_circles(d2>10, 50, 300)>0) | (find_circles(d2>7, 50, 300)>0);
else
%     d1=(find_circles(d>init_thre, 100, 500)>0) | (find_circles(d>init_thre1, 100, 1000)>0); 	% red beads in AP
%     if ap==0
%         d1=d1|(find_circles(d>30, 50, 500)>0) | (find_circles(d>150, 100, 500)>0);		% green beads in TRAP
%     end
%     d1=(find_circles(d>init_thre, 100, 500)>0) | (find_circles(d>init_thre1, 100, 1000)>0) | (find_circles(d>20, 100, 5000)>0) | (find_circles(imclose(d>250,strel('disk',5)), 100, 5000)>0); 	% red beads in AP
%     d1=(find_circles(d>init_thre, 100, 500)>0) | (find_circles(d>init_thre1, 100, 500)>0) | (find_circles(d>20, 100, 500)>0) | (find_circles(imclose(d>25,strel('disk',5)), 100, 500)>0); 	% red beads in AP
%     if ap==0
%         d1=d1|(find_circles(d>30, 50, 500)>0) | (find_circles(d>150, 100, 500)>0);		% green beads in TRAP
        d1=(find_circles(d>init_thre, 100, 300)>0) | (find_circles(d>init_thre1, 100, 300)>0)|(find_circles(d>50, 50, 300)>0) | (find_circles(d>40, 50, 300)>0) | (find_circles(d>30, 50, 300)>0) | (find_circles(d>20, 50, 300)>0) | (find_circles(d>15, 50, 300)>0) | (find_circles(d>10, 50, 300)>0) | (find_circles(d>7, 50, 300)>0);
        d1=d1|(find_circles(d>50, 50, 1000)>0) | (find_circles(d>40, 50, 2000)>0) | (find_circles(d>30, 50, 3000)>0) | (find_circles(d>20, 50, 4000)>0) | (find_circles(d>15, 50, 5000)>0);
%     end
%     d1=(find_circles(d>init_thre1, 100, 1000)>0);
end
if length(find(d))<1000
    reg_info(1,1:2)=[100000,100000];
    no_bead=1;
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Remove beads in the middle of Vertebra to avoid 'no second BoundingBox'
%
%   August 13, 2016
%
temp=sum(d1,2);
m=[];
for i=2000:100:min(4000,floor(size(d1,1)/100)*100-1000)
    ind=round(i/100)-19;
    m=[m, sum(temp(i:i+1000))];
end
ind=find(m==min(m));
if ~isempty(ind)
    if isempty(find(m==0))
        if length(ind)~=1
            ind=ind(floor(length(ind)/2));
        end
        d1((ind+19)*100:(ind+19)*100+1000,:)=false(1001,size(d1,2));
    end
end
temp=sum(d1,1);
m=[];
% for i=2000:100:min(4000,floor(size(d1,2)/100)*100-1000)
%     ind=round(i/100)-19;
%     m=[m, sum(temp(i:i+1000))];
% end
% ind=find(m==min(m));
% if ~isempty(ind)
%     if isempty(find(m==0))
%         if length(ind)~=1
%             ind=ind(floor(length(ind)/2));
%         end
%         d1(:,(ind+19)*100:(ind+19)*100+1000)=false(size(d1,1),1001);
%     end
% end
%
%   Remove beads in the middle of Vertebra to avoid 'no second BoundingBox'
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure;imshow(a1);
% figure;imshow(d1);
% close all

max_x=max(size(a1,1),size(d1,1));
max_y=max(size(a1,2),size(d1,2));
a1(max_x,max_y)=0;
d1(max_x,max_y)=0;
a(max_x,max_y)=0;
d(max_x,max_y)=0;

size_factor=5;
aa=imresize(a1,1/size_factor);
dd=imresize(d1,1/size_factor);
max_x=max(size(aa,1),size(dd,1));
max_y=max(size(aa,2),size(dd,2));
aa(max_x,max_y)=0;
dd(max_x,max_y)=0;
size_x=size(aa,1);
size_y=size(aa,2);

ref=zeros(size_x*2-1,size_y*2-1);

test=ref;
half_x=round(size_x/2);
half_y=round(size_y/2);
ref(half_x:half_x+size_x-1,half_y:half_y+size_y-1)=double(aa); clear aa
test(half_x:half_x+size_x-1,half_y:half_y+size_y-1)=double(dd); clear dd
R=fft2(ref); clear ref;
T=fft2(test); clear test
aR=abs(R);
pR=angle(R);
aT=abs(T);
pT=angle(T); clear R T
RT=aR.*aT.*(exp(j*(pT-pR))); clear aR aT pR pT
rt=ifft2(RT); clear RT
srt=fftshift(rt);
[x3 y3]=find(srt==max(srt(:)));
dis_v=(size_x-(x3(round(length(x3)/2))))*size_factor;
dis_h=(size_y-(y3(round(length(y3)/2))))*size_factor;

if dis_v>0
    a(1:end+dis_v,end)=0;
    temp=uint8(false(size(a)));
    temp(dis_v+1:end,:)=imadd(temp(dis_v+1:end,:),d); d=temp;
else
    d(1:end-dis_v,end)=0;
    temp=uint8(false(size(d)));
    temp(-dis_v+1:end,:)=imadd(temp(-dis_v+1:end,:),a); a=temp;
end
if dis_h>0
    a(end,1:end+dis_h)=0;
    temp=uint8(false(size(a)));
    temp(:,dis_h+1:end)=imadd(temp(:,dis_h+1:end),d); d=temp;
else
    d(end,1:end-dis_h)=0;
    temp=uint8(false(size(d)));
    temp(:,-dis_h+1:end)=imadd(temp(:,-dis_h+1:end),a); a=temp;
end

if dis_v>0
    a1(1:end+dis_v,end)=0;
    temp=false(size(a1));
    temp(dis_v+1:end,:)=imadd(temp(dis_v+1:end,:),d1); d1=temp;
else
    d1(1:end-dis_v,end)=0;
    temp=false(size(d1));
    temp(-dis_v+1:end,:)=imadd(temp(-dis_v+1:end,:),a1); a1=temp;
end
if dis_h>0
    a1(end,1:end+dis_h)=0;
    temp=false(size(a1));
    temp(:,dis_h+1:end)=imadd(temp(:,dis_h+1:end),d1); d1=temp;
else
    d1(end,1:end-dis_h)=0;
    temp=false(size(d1));
    temp(:,-dis_h+1:end)=imadd(temp(:,-dis_h+1:end),a1); a1=temp;
end
clear temp
reg_info(1,1:2)=[dis_v, dis_h];

% [L, na]=bwlabel(a1);
% [L, nd]=bwlabel(d1);
% 
% dummy=1;
% count=1;
% recalc_a=0;
% recalc_d=0;
% while dummy
%     if count>15
%         break
%     end
%     thre=init_thre-count*10;
% %     if (na/nd>1.5 & count==1) | (na/nd>1.5 & count>1  & recalc_d==1)
%     if (na/nd>1.5 & count==1 & nd<20) | (na/nd>1.5 & count>1  & nd<20 & recalc_d==1)
%         d1=find_circles(d>thre, 30, 200)>0;
%         [L, nd]=bwlabel(d1);
%         count=count+1;
%         recalc_d=1;
%     elseif (na/nd<0.75 & count==1 & na<20) | (na/nd<0.75 & count>1 & na<20 & recalc_a==1)
%         a1=find_circles(a>thre, 30, 200)>0;
%         [L, na]=bwlabel(a1);
%         count=count+1;
%         recalc_a=1;
%     else
%         dummy=0;
%     end
% end
%     
% if recalc_a==1
%     a=a>thre;
% elseif recalc_a==0
%     a=a>init_thre;
% end
% if recalc_d==1
%     d=d>thre;
% elseif recalc_d==0
%     d=d>init_thre;
% end


% if gutta==1
%     [L n]=bwlabel(a);
%     stat_a=regionprops(L, 'BoundingBox','Image','Area'); clear L
%     ar_a=[stat_a.Area];
%     ref_cand=find(ar_a>100000);
%     if isempty(ref_cand)
%         ref_cand=find(ar_r==max(ar_r));
%     end
% end
% 
% if gutta==1
%     [L n]=bwlabel(d);
%     stat_d=regionprops(L, 'BoundingBox','Image','Area'); clear L
%     ar_d=[stat_d.Area];
%     im_cand=find(ar_d>100000);
%     if isempty(im_cand)
%         im_cand=find(ar_d==max(ar_d));
%     end
% end
% 
% if gutta==1
%     for i=1:length(ref_cand)
%         for jj=1:length(im_cand)
%             temp=stat_a(ref_cand(i)).BoundingBox./stat_d(im_cand(jj)).BoundingBox;
%             ratio(i,jj,1:2)=temp(3:4);
%             temp1=stat_a(ref_cand(i)).Area/stat_d(im_cand(jj)).Area;
%             ratio(i,jj,3)=temp1;
%             ratio1(i,jj)=prod(temp)*temp1;
%         end
%     end
%     % test=abs(ratio1-1);
%     % test=abs(sum(ratio,3)-3);
%     % [x y]=find(test==min(test(:)));
%     tt=ratio;
%     tt(tt>1.15 |tt<0.85)=0;
%     tt1=abs(prod(tt,3)-1);
%     [x y]=find(tt1==min(tt1(:)));
% 
%     dis_v=stat_a(ref_cand(x)).BoundingBox(2)-stat_d(im_cand(y)).BoundingBox(2);
%     dis_h=stat_a(ref_cand(x)).BoundingBox(1)-stat_d(im_cand(y)).BoundingBox(1);
% 
%     max_x=max(size(a,2),size(d,2));
%     max_y=max(size(a,1),size(d,1));
%     a(max_y,max_x)=0;
%     d(max_y,max_x)=0;
% 
%     if dis_v>0
%         a(1:end+dis_v,end)=0;
%         temp=false(size(a));
%         temp(dis_v+1:end,:)=d; d=temp;
%     else
%         d(1:end-dis_v,end)=0;
%         temp=false(size(d));
%         temp(-dis_v+1:end,:)=a; a=temp;
%     end
%     if dis_h>0
%         a(end,1:end+dis_h)=0;
%         temp=false(size(a));
%         temp(:,dis_h+1:end)=d; d=temp;
%     else
%         d(end,1:end-dis_h)=0;
%         temp=false(size(d));
%         temp(:,-dis_h+1:end)=a; a=temp;
%     end
% 
%     reg_info(1,1:2)=[dis_v, dis_h];
% else
% %     max_x=max(size(a,2),size(d,2));
% %     max_y=max(size(a,1),size(d,1));
% %     a(max_y,max_x)=0;
% %     d(max_y,max_x)=0;
% %     reg_info(1,1:2)=[0, 0];
% end

% a1=find_circles(a, 30, 200)>0;
% d1=find_circles(d, 30, 200)>0;

ta1=imdilate(a1,strel('disk',100));
td1=imdilate(d1,strel('disk',100));
%temp=find_circles(ta1&td1, 6000, 20000)>0;
temp=ta1&td1;

% disk_size=[250, 300];
disk_size=[100, 150, 200, 250, 300];
for i=1:length(disk_size)
    g2=imdilate(temp,strel('disk',disk_size(i)));
    [L, n]=bwlabel(g2); clear g2
    stat=regionprops(L, 'BoundingBox','Image','Area'); clear L
    ar=[stat.Area];
    sar=sort(ar);
    if length(sar)>1
        if sar(end-1)/sum(sar)>0.1 
            break
        end
    end
end
% a2=false(size(a1));
ind=find(ar==sar(end));
st_x=round(stat(ind).BoundingBox(2))+round(disk_size(i)/2);
st_y=round(stat(ind).BoundingBox(1))+round(disk_size(i)/2);
ed_x=max(3000,stat(ind).BoundingBox(4))+st_x-1-round(disk_size(i)/2);
ed_y=max(3000,stat(ind).BoundingBox(3))+st_y-1-round(disk_size(i)/2);
% a2(st_x:ed_x,st_y:ed_y)=stat(ind).Image;
first_Bounding_Box=[st_x,ed_x, st_y, ed_y];
for count=1:length(sar)-1
    if length(sar)>1
        ind=find(ar==sar(end-count));
%         st_x=round(stat(ind).BoundingBox(2))+round(disk_size(i)/2);
%         st_y=round(stat(ind).BoundingBox(1))+round(disk_size(i)/2);
%         if abs(st_x-first_Bounding_Box(1))<2000 & abs(st_y-first_Bounding_Box(2))<2000
%             continue
%         end
%         ed_x=min(3000,stat(ind).BoundingBox(4))+st_x-1-round(disk_size(i)/2);
%         ed_y=min(3000,stat(ind).BoundingBox(3))+st_y-1-round(disk_size(i)/2);
%         % a2(st_x:ed_x,st_y:ed_y)=stat(ind).Image;
%         second_Bounding_Box=[st_x,ed_x, st_y, ed_y];
%         break
        st_x=round(stat(ind).BoundingBox(2))+round(disk_size(i)/2);
        st_y=round(stat(ind).BoundingBox(1))+round(disk_size(i)/2);
        ed_x=stat(ind).BoundingBox(4)+st_x-1-round(disk_size(i)/2);
        ed_y=stat(ind).BoundingBox(3)+st_y-1-round(disk_size(i)/2);
        if st_x-first_Bounding_Box(2)>500 | first_Bounding_Box(1)-ed_x>500 |st_y-first_Bounding_Box(4)>500 | first_Bounding_Box(3)-ed_y>500
            second_Bounding_Box=[st_x,ed_x, st_y, ed_y];
            break
        else
            continue
        end
    else
        second_Bounding_Box=[size(d1,1)-100,size(d1,1), size(d1,2)-100,size(d1,2)];
    end
end
if ~exist('second_Bounding_Box', 'var')
    no_bead=1;
    ['No second_Bounding_Box']
    
    red=a1;
    AP=d1;

    temp=false(size(a1));
    temp(first_Bounding_Box(1):first_Bounding_Box(2),first_Bounding_Box(3):first_Bounding_Box(4))...
       =true(first_Bounding_Box(2)-first_Bounding_Box(1)+1,first_Bounding_Box(4)-first_Bounding_Box(3)+1);
    temp=logical(imdilate(temp,strel('disk',20))-temp);

    aa=a1|temp;
    aa(:,:,2)=d1|temp;
    aa(:,:,3)=0; clear temp;
    error_code=1;
    return
end   
if first_Bounding_Box(1)>second_Bounding_Box(1)
    temp=second_Bounding_Box;
    second_Bounding_Box=first_Bounding_Box;
    first_Bounding_Box=temp;
end
cont=1;
count=2;
% while cont
%     if second_Bounding_Box(1)-first_Bounding_Box(1)<2000 & n>2
%         ind=find(ar==sar(end-count));
%         count=count+1;
%         st_x=round(stat(ind).BoundingBox(2))+round(disk_size(i)/2);
%         st_y=round(stat(ind).BoundingBox(1))+round(disk_size(i)/2);
%         ed_x=min(3000,stat(ind).BoundingBox(4))+st_x-1-round(disk_size(i)/2);
%         ed_y=min(3000,stat(ind).BoundingBox(3))+st_y-1-round(disk_size(i)/2);
%         if first_Bounding_Box(1)>5000
%             first_Bounding_Box=[st_x,ed_x, st_y, ed_y];
%         elseif second_Bounding_Box(1)<3000
%             second_Bounding_Box=[st_x,ed_x, st_y, ed_y];
%         end
%     else
%         cont=0;
%     end
%     if count==n
%         cont=0;
%     end
% end

red=a1;
AP=d1;

temp=false(size(a1));
temp(first_Bounding_Box(1):first_Bounding_Box(2),first_Bounding_Box(3):first_Bounding_Box(4))...
   =true(first_Bounding_Box(2)-first_Bounding_Box(1)+1,first_Bounding_Box(4)-first_Bounding_Box(3)+1);
temp(second_Bounding_Box(1):second_Bounding_Box(2),second_Bounding_Box(3):second_Bounding_Box(4))...
   =true(second_Bounding_Box(2)-second_Bounding_Box(1)+1,second_Bounding_Box(4)-second_Bounding_Box(3)+1);
temp=logical(imdilate(temp,strel('disk',20))-temp);

aa=a1|temp;
aa(:,:,2)=d1|temp;
aa(:,:,3)=0; clear temp;
aa=uint8(aa)*255;
% figure;imshow(uint8(aa)*255);axis on;grid
% first_Bounding_Box
% second_Bounding_Box

count=1;
trap_done=0;
break_while=0;
angle_difference=90;
angle_before=0;
dis_v=10000; dis_h=10000; ratio_v=2; ratio_h=2; angle_between_green_trap=90; angle_between_green_trap1=90;

while abs(angle_between_green_trap)>0.01 & abs(angle_between_green_trap1)>0.01
    if break_while==1
        break
    end
    if abs(angle_difference)<0.5
        trap_done=1;
        break
    end
    for dummy=1
        [dis_v, dis_h, ratio_v, ratio_h, angle_between_green_trap, no_bead]=match_points_using_fft(red, AP, first_Bounding_Box, second_Bounding_Box);

        if ~isempty(find(no_bead, 1)) | abs(angle_between_green_trap)>10 | ratio_v>1.3 | ratio_h>1.3
            dis_v=0; dis_h=0; ratio_v=1; ratio_h=1; angle_between_green_trap=90;
            break_while=1;
            [no_bead, angle_between_green_trap, ratio_v, ratio_h]
            ['TRAP/AP Registration Error 1']
%             continue
            no_bead=1;
            if ~isempty(find(no_bead, 1))
                error_code=2.1;
            end
            if abs(angle_between_green_trap)>10
                error_code=2.2;
            end
            if ratio_v>1.3 | ratio_h>1.3
                error_code=2.3;
            end
            break
        end
        AP=imrotate(AP,-angle_between_green_trap);
        max_x=max([size(red,2),size(AP,2)]);
        max_y=max([size(red,1),size(AP,1)]);
        red(max_y,max_x)=0;
        AP(max_y,max_x)=0;

        [dis_v, dis_h, ratio_v, ratio_h, angle_between_green_trap1, no_bead]=match_points_using_fft(red, AP, first_Bounding_Box, second_Bounding_Box);
        if ~isempty(find(no_bead, 1))
            dis_v=0; dis_h=0; ratio_v=1; ratio_h=1; angle_between_green_trap1=90;
            break_while=1;
            ['TRAP/AP Registration Error 2']
%             continue
            no_bead=1;
            error_code=3;
            break
        end
        AP=imresize(AP,[round(ratio_v*size(AP,1)) round(ratio_h*size(AP,2))]);
        max_x=max(size(red,2),size(AP,2));
        max_y=max(size(red,1),size(AP,1));
        red(max_y,max_x)=0;
        AP(max_y,max_x)=0;
        if dis_v>0
            red(1:end+dis_v,end)=0;
            temp=false(size(red));
            temp(dis_v+1:end,:)=imadd(temp(dis_v+1:end,:),AP); AP=temp;
        else
            AP(1:end-dis_v,end)=0;
            temp=false(size(AP));
            temp(-dis_v+1:end,:)=imadd(temp(-dis_v+1:end,:),red); red=temp;
        end
        if dis_h>0
            red(end,1:end+dis_h)=0;
            temp=false(size(red));
            temp(:,dis_h+1:end)=imadd(temp(:,dis_h+1:end),AP); AP=temp;
        else
            AP(end,1:end-dis_h)=0;
            temp=false(size(AP));
            temp(:,-dis_h+1:end)=imadd(temp(:,-dis_h+1:end),red); red=temp;
        end
        trap_done=1;
    end
    if trap_done==1
        reg_info(count+1,1:5)=[angle_between_green_trap, ratio_v, ratio_h, dis_v, dis_h];
        count=count+1;
    end
    angle_difference=angle_between_green_trap-angle_before;
    angle_before=angle_between_green_trap;
end % while
% ['end of TRAP at ', datestr(now)]

return

function [dis_x1, dis_y1, ratio_x, ratio_y, angle_r_minus_a, no_bead]=match_points_using_fft(red, AP, first_Bounding_Box, second_Bounding_Box)

% bead_size_red=125;
% bead_margin_red=75;
% bead_size_AP=125;
% bead_margin_AP=75;
max_distance_between_beeads=500;
for i=1:2
    switch i
        case 1
            BB=first_Bounding_Box;
            region=5;
        case 2
            BB=second_Bounding_Box;
            region=6;
    end
    cor_red=red(BB(1):BB(2),BB(3):BB(4));
%     cor_red=(cor_red>threshold(1));
%     [cor_red]=find_circles(cor_red, bead_size_red-bead_margin_red, bead_size_red+bead_margin_red);

    cor_ap=AP(BB(1):BB(2),BB(3):BB(4));
%     cor_ap=(cor_ap>threshold(2));
    %cor_ap=imerode(imdilate(cor_ap, strel('disk',dilation')),strel('disk',dilation'));
%     [cor_ap]=find_circles(cor_ap, bead_size_AP-bead_margin_AP, bead_size_AP+bead_margin_AP);

%     [dis_x(5), dis_y(5), center_xr(5), center_yr(5), center_xa(5), center_ya(5), no_bead(5), n(5)]=find_center(cor_red, cor_ap, max_distance_between_beeads, offset_x, offset_y, region);
    [dis_x(i), dis_y(i), center_xr(i), center_yr(i), center_xa(i), center_ya(i), no_bead(i), n(i)]=find_center(cor_red, cor_ap, max_distance_between_beeads, BB(1), BB(3), region);
end

ref_region=1;
opp_region=2;
if center_xr(opp_region)>center_xr(ref_region)
    center_x_denom=center_xr(opp_region);
    dis_xx=dis_x(opp_region)-dis_x(ref_region);
    dis_x1=dis_x(ref_region);
else
    center_x_denom=center_xr(ref_region);
    dis_xx=dis_x(ref_region)-dis_x(opp_region);
    dis_x1=dis_x(opp_region);
end
if center_yr(opp_region)>center_yr(ref_region)
    center_y_denom=center_yr(opp_region);
    dis_yy=dis_y(opp_region)-dis_y(ref_region);
    dis_y1=dis_y(ref_region);
else
    center_y_denom=center_yr(ref_region);
    dis_yy=dis_y(ref_region)-dis_y(opp_region);
    dis_y1=dis_y(opp_region);
end
ratio_x=1+dis_xx(1)/center_x_denom;
ratio_y=1+dis_yy(1)/center_y_denom;
angle_a=atan((center_ya(opp_region)-center_ya(ref_region))/(center_xa(opp_region)-center_xa(ref_region)));
angle_r=atan((center_yr(opp_region)-center_yr(ref_region))/(center_xr(opp_region)-center_xr(ref_region)));
angle_r_minus_a=(angle_a-angle_r)/pi*180;

return

function [g1]=find_circles(g, min_size, max_size)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
%
%       input : g (input with beads)
%               min_size( minimum size of beads)
%               max_size( maximum size of beads)
%
%       output : g1 (image with found beads)
%
%       date  : January 6, 2010
%       author : Sean Hong
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55

[L n]=bwlabel(g);
stats=regionprops(L,'BoundingBox', 'Image', 'Area'); clear L
a=[stats.Area];
x=find(a>min_size & a<max_size);
g1=zeros(size(g));
for i=1:length(x)
    sx=round(stats(x(i)).BoundingBox(2));
    sy=round(stats(x(i)).BoundingBox(1));
    ex=sx+stats(x(i)).BoundingBox(4)-1;
    ey=sy+stats(x(i)).BoundingBox(3)-1;
    test=a(x(i))/(((max(stats(x(i)).BoundingBox(4), stats(x(i)).BoundingBox(3)))/2)^2*pi);
    radius_ratio=stats(x(i)).BoundingBox(4)/stats(x(i)).BoundingBox(3);
%    if (test>0.75 & test <1.25) & (radius_ratio>0.75 & radius_ratio <1.25)
    if (test>0.8 & test <1.2) & (radius_ratio>0.85 & radius_ratio <1.15)
        g1(sx:ex,sy:ey)=stats(x(i)).Image;
    end
end
return

function [dis_x2, dis_y2, center_xr2, center_yr2, center_xa2, center_ya2, no_bead, n]=find_center(cor_red, cor_ap, max_distance_between_beeads, offset_x, offset_y, region)

min_beads=2;
no_bead=0;
[xr_temp yr_temp]=find(cor_red);
[L n1]=bwlabel(cor_red);
if n1<min_beads
    no_bead=1;
    dis_x2=0;
    dis_y2=0;
    center_xr2=0;
    center_yr2=0;
    center_xa2=0;
    center_ya2=0;
    n=0;
    return
end
stats=regionprops(L, 'Centroid');
for i=1:n1
    yr(i)=round(stats(i).Centroid(1));      % horizontal
    xr(i)=round(stats(i).Centroid(2));      % vertical
end
s_xr=sort(xr);
Ix=find(abs(diff(s_xr))>max_distance_between_beeads);
s_yr=sort(yr);
Iy=find(abs(diff(s_yr))>max_distance_between_beeads);
clear tt no
if isempty(Ix)==0 & isempty(Iy)==1
    Ixx=[0,Ix];
    for i=2:length(Ixx)
        no(i-1)=Ixx(i)-Ixx(i-1);
    end
    no(i)=length(s_xr)-Ixx(i);
    loc=find(no==max(no));
    s_xr2=s_xr(Ixx(loc)+1:Ixx(loc)+no(loc));
    for i=1:length(s_xr2)
        temp=find(xr==s_xr2(i));
        tt(i)=temp(1);
    end
    center_xr2=sum(xr(tt))/length(tt)+offset_x;
    center_yr2=sum(yr(tt))/length(tt)+offset_y;
%     switch region          % 1: left upper, 2: right upper, 3: bottom left, 4: bottom right 
%         case 1
%             center_xr2=sum(xr(tt))/length(tt);
%             center_yr2=sum(yr(tt))/length(tt);
%         case 2
% %            center_xr2=sum(xr(tt))/length(tt)+offset_x;
%             center_xr2=sum(xr(tt))/length(tt);
%             center_yr2=sum(yr(tt))/length(tt)+offset_y;
%         case 3
%             center_xr2=sum(xr(tt))/length(tt)+offset_x;
%             center_yr2=sum(yr(tt))/length(tt);
%         case 4
%             center_xr2=sum(xr(tt))/length(tt)+offset_x;
%             center_yr2=sum(yr(tt))/length(tt)+offset_y;
%     end
elseif isempty(Ix)==1 & isempty(Iy)==0
    Iyy=[0,Iy];
    for i=2:length(Iyy)
        no(i-1)=Iyy(i)-Iyy(i-1);
    end
    no(i)=length(s_yr)-Iyy(i);
    loc=find(no==max(no));
    s_yr2=s_yr(Iyy(loc)+1:Iyy(loc)+no(loc));
    for i=1:length(s_yr2)
        temp=find(yr==s_yr2(i));
        tt(i)=temp(1);
    end
    center_xr2=sum(xr(tt))/length(tt)+offset_x;
    center_yr2=sum(yr(tt))/length(tt)+offset_y;
%     switch region          % 1: left upper, 2: right upper, 3: bottom left, 4: bottom right 
%         case 1
%             center_xr2=sum(xr(tt))/length(tt);
%             center_yr2=sum(yr(tt))/length(tt);
%         case 2
% %            center_xr2=sum(xr(tt))/length(tt)+offset_x;
%             center_xr2=sum(xr(tt))/length(tt);
%             center_yr2=sum(yr(tt))/length(tt)+offset_y;
%         case 3
%             center_xr2=sum(xr(tt))/length(tt)+offset_x;
%             center_yr2=sum(yr(tt))/length(tt);
%         case 4
%             center_xr2=sum(xr(tt))/length(tt)+offset_x;
%             center_yr2=sum(yr(tt))/length(tt)+offset_y;
%     end
elseif isempty(Ix)==0 && isempty(Iy)==0
    Ixx=[0,Ix];
    for i=2:length(Ixx)
        no(i-1)=Ixx(i)-Ixx(i-1);
    end
    no(i)=length(s_xr)-Ixx(i);
    loc=find(no==max(no));
    s_xr2=s_xr(Ixx(loc)+1:Ixx(loc)+no(loc));
    for i=1:length(s_xr2)
        temp=find(xr==s_xr2(i));
        tt(i)=temp(1);
    end
    clear tt1 no

    Iyy=[0,Iy];
    for i=2:length(Iyy)
        no(i-1)=Iyy(i)-Iyy(i-1);
    end
    no(i)=length(s_yr)-Iyy(i);
    loc=find(no==max(no));
    s_yr2=s_yr(Iyy(loc)+1:Iyy(loc)+no(loc));
    for i=1:length(s_yr2)
        temp=find(yr==s_yr2(i));
        tt1(i)=temp(1);
    end
    tt2=intersect(tt, tt1);
    center_xr2=sum(xr(tt2))/length(tt2)+offset_x;
    center_yr2=sum(yr(tt2))/length(tt2)+offset_y;
%     switch region          % 1: left upper, 2: right upper, 3: bottom left, 4: bottom right 
%         case 1
%             center_xr2=sum(xr(tt2))/length(tt2);
%             center_yr2=sum(yr(tt2))/length(tt2);
%         case 2
% %            center_xr2=sum(xr(tt2))/length(tt2)+offset_x;
%             center_xr2=sum(xr(tt2))/length(tt2);
%             center_yr2=sum(yr(tt2))/length(tt2)+offset_y;
%         case 3
%             center_xr2=sum(xr(tt2))/length(tt2)+offset_x;
%             center_yr2=sum(yr(tt2))/length(tt2);
%         case 4
%             center_xr2=sum(xr(tt2))/length(tt2)+offset_x;
%             center_yr2=sum(yr(tt2))/length(tt2)+offset_y;
%     end
else
    center_xr2=mean(xr)+offset_x;
    center_yr2=mean(yr)+offset_y;
%     switch region          % 1: left upper, 2: right upper, 3: bottom left, 4: bottom right 
%         case 1
%             center_xr2=mean(xr);
%             center_yr2=mean(yr);
%         case 2
% %            center_xr2=mean(xr)+offset_x;
%             center_xr2=mean(xr);
%             center_yr2=mean(yr)+offset_y;
%         case 3
%             center_xr2=mean(xr)+offset_x;
%             center_yr2=mean(yr);
%         case 4
%             center_xr2=mean(xr)+offset_x;
%             center_yr2=mean(yr)+offset_y;
%     end
end

[xa_temp ya_temp]=find(cor_ap);
[L n2]=bwlabel(cor_ap);
if n2<min_beads
    no_bead=1;
    dis_x2=0;
    dis_y2=0;
    center_xr2=0;
    center_yr2=0;
    center_xa2=0;
    center_ya2=0;
    n=0;
    return
end


% cor_red(max([xr_temp;xa_temp]),max([yr_temp;ya_temp]))=0;
% cor_red=cor_red(1:max([xr_temp;xa_temp]),1:max([yr_temp;ya_temp]));
% cor_ap(max([xr_temp;xa_temp]),max([yr_temp;ya_temp]))=0;
% cor_ap=cor_ap(1:max([xr_temp;xa_temp]),1:max([yr_temp;ya_temp]));
% 
% cor_ap(size(cor_red,1)*2-1,size(cor_red,2)*2-1)=0;
% cor_red(size(cor_red,1)*2-1,size(cor_red,2)*2-1)=0;
% 
% A=fft2(cor_ap);
% ma_A=abs(A);
% ph_A=angle(A); clear A
% 
% R=fft2(cor_red);
% ma_R=abs(R);
% ph_R=angle(R); clear R
% AR=((ma_A.*ma_R).^0.3).*exp(j*(ph_R-ph_A));
% clear ma* ph*
% 
% ar=fftshift(abs(ifft2(AR))); clear AR
% [x y]=find(ar==max(ar(:)));
% x=x(round(length(x)/2));
% y=y(round(length(y)/2));
% 
% temp_r=abs(xr-x)+abs(yr-y);
% ind_r=find(temp_r==min(temp_r));

[L n1]=bwlabel(cor_ap);
stats=regionprops(L, 'Centroid');
for i=1:n1
    ya(i)=round(stats(i).Centroid(1));      % horizontal
    xa(i)=round(stats(i).Centroid(2));      % vertical
end
% temp_a=abs(xa-x)+abs(ya-y);
% ind_a=find(temp_a==min(temp_a));
% 
% % dis_x2=x-max([xr_temp;xa_temp]);
% % dis_y2=y-max([yr_temp;ya_temp]);
% dis_x2=xr(ind_r)-xa(ind_a);
% dis_y2=yr(ind_r)-ya(ind_a);

% [dis_x2, dis_y2]=find_dis_wo_fft(xr, xa, yr, ya);

% aa=xcorr2(double(cor_red),double(cor_ap));
% [x y]=find(aa==max(aa(:)));
% dis_x2=(x(round(length(x)/2))-size(cor_red,1));
% dis_y2=(y(round(length(y)/2))-size(cor_red,2));
ref=zeros(size(cor_red,1)*2-1,size(cor_red,2)*2-1);
test=ref;
half_x=round(size(cor_red,1)/2);
half_y=round(size(cor_red,2)/2);

ref(half_x:half_x+size(cor_red,1)-1,half_y:half_y+size(cor_red,2)-1)=double(cor_red);
test(half_x:half_x+size(cor_red,1)-1,half_y:half_y+size(cor_red,2)-1)=double(cor_ap);
R=fft2(ref);
T=fft2(test);
aR=abs(R);
pR=angle(R);
aT=abs(T);
pT=angle(T);
RT=aR.*aT.*(exp(j*(pT-pR)));
rt=ifft2(RT);
srt=fftshift(rt);
[x3 y3]=find(srt==max(srt(:)));
dis_x2=size(cor_red,1)-(x3(round(length(x3)/2)));
dis_y2=size(cor_red,2)-(y3(round(length(y3)/2)));


center_xa2=center_xr2-dis_x2;
center_ya2=center_yr2-dis_y2;

n=min(n1,n2);
return

function [dis_x1, dis_y1, ratio_x, ratio_y, angle_r_minus_a, no_bead1]=matching_points_fft(red, AP, trap_ap, DIC, bead_size_red, bead_margin_red, bead_size_AP, bead_margin_AP, threshold, rot_angle, bone_type)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   FInd the rough adjacent matching points between red and AP
%
%   Input  : red1, AP	(red, AP images)
%            trap_ap    (1 --> trap, 0 --> ap)
%            bead_size, bead_size_margin
%   Output : dis_x, dis_y	(number of pixels of AP image shift; vertical, horizontal)
%
%   Date   : January, 27, 2010
%   Author : Sean Hong
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
%   left upper corner
%
no_bead=0;
dis_x1=0;
dis_y1=0;
ratio_x=0;
ratio_y=0;
angle_r_minus_a=0;

ind_special=0;      % find reference and test region from 4 corners
                    % if 1, reference region is one of the middle (upper or lower) and test is one of the 4 corners


if trap_ap==0       % ap
    dilation=15;
elseif trap_ap==1;  % trap
    dilation=3;
end
% if trap_ap==0
%     threshold=70;
% else
%     threshold=12;
% end

%threshold=40;

max_distance_between_beeads=800;
%AP=immultiply(AP,temp);


%if rot_angle~=0     % upper middle, bottom middle, upper left corner, upper right corner, bottom left corner, bottom right corner
if strcmp(bone_type,'F')==1 | strcmp(bone_type,'V')==1     % upper middle, bottom middle, upper left corner, upper right corner, bottom left corner, bottom right corner

    %
    %   upper middle
    %

    %cor_red=red1(1:round(size(red1,1)/4*3),1:round(size(red1,2)/4));
    offset_x=0;
    offset_y=round(size(red,2)/4);
    cor_red=red(1:round(size(red,1)/3),offset_y+1:round(size(red,2)*3/4));
    cor_red=(cor_red>threshold(1));
    [cor_red]=find_circles(cor_red, bead_size_red-bead_margin_red, bead_size_red+bead_margin_red);

    cor_ap=AP(1:round(size(red,1)/3),offset_y+1:round(size(red,2)*3/4));
    cor_ap=(cor_ap>threshold(2));
    %cor_ap=imerode(imdilate(cor_ap, strel('disk',dilation')),strel('disk',dilation'));
    [cor_ap]=find_circles(cor_ap, bead_size_AP-bead_margin_AP, bead_size_AP+bead_margin_AP);

    region=5;
    [dis_x(5), dis_y(5), center_xr(5), center_yr(5), center_xa(5), center_ya(5), no_bead(5), n(5)]=find_center(cor_red, cor_ap, max_distance_between_beeads, offset_x, offset_y, region);

    if n(5)~=0
    %
    %   bottom middle
    %
        offset_x=round(size(red,1)*3/4);
        offset_y=round(size(red,2)/4);
        cor_red=red(offset_x+1:end,offset_y+1:round(size(red,2)*3/4));
        cor_red=(cor_red>threshold(1));
        [cor_red]=find_circles(cor_red, bead_size_red-bead_margin_red, bead_size_red+bead_margin_red);

        cor_ap=AP(offset_x+1:end,offset_y+1:round(size(red,2)*3/4));
        cor_ap=(cor_ap>threshold(2));
        %cor_ap=imerode(imdilate(cor_ap, strel('disk',dilation')),strel('disk',dilation'));
        [cor_ap]=find_circles(cor_ap, bead_size_AP-bead_margin_AP, bead_size_AP+bead_margin_AP);

        region=6;
        [dis_x(6), dis_y(6), center_xr(6), center_yr(6), center_xa(6), center_ya(6), no_bead(6), n(6)]=find_center(cor_red, cor_ap, max_distance_between_beeads, offset_x, offset_y, region);

    else
        dis_x(6)=0; dis_y(6)=0; center_xr(6)=0; center_yr(6)=0; center_xa(6)=0; center_ya(6)=0; no_bead(6)=1; n(6)=0;
    end

    if n(5)~=0 & n(6)~=0 & abs(dis_x(6)-dis_x(5))<100 & abs(dis_y(6)-dis_y(5))<100
%     if n(5)~=0 & n(6)>25 %& abs(dis_x(6)-dis_x(5))<100 & abs(dis_y(6)-dis_y(5))<100
        ref_region=5;
        opp_region=6;
        no_bead1=0;
    else
        a=(DIC>graythresh(DIC)*255);
        a=imerode(imdilate(a, strel('disk',31')),strel('disk',31'));
        [L n1]=bwlabel(a); clear a
        stats=regionprops(L, 'Area','BoundingBox','FIlledImage'); clear L
        area_dic=[stats.Area];
        [I J]=find(area_dic>200);
        temp=zeros(size(DIC));
        for i=1:size(J,2)
            sx=round(stats(J(i)).BoundingBox(2));
            sy=round(stats(J(i)).BoundingBox(1));
            ex=sx+stats(J(i)).BoundingBox(4)-1;
            ey=sy+stats(J(i)).BoundingBox(3)-1;
            temp(sx:ex,sy:ey)=temp(sx:ex,sy:ey)+stats(J(i)).FilledImage;
        end
        temp=1-temp;
        temp=uint8(temp>0);
        red1=immultiply(red,temp); clear temp



        red1=red;
    %
    %   upper left corner
    %

    %cor_red=red1(1:round(size(red1,1)/4*3),1:round(size(red1,2)/4));
        offset_x=0;
        offset_y=0;

        cor_red=red1(1:round(size(red1,1)/2),1:round(size(red1,2)/4));
        cor_red=(cor_red>threshold(1));
        [cor_red]=find_circles(cor_red, bead_size_red-bead_margin_red, bead_size_red+bead_margin_red);

        cor_ap=AP(1:round(size(red1,1)/2),1:round(size(red1,2)/4));
        cor_ap=(cor_ap>threshold(2));
        cor_ap=imerode(imdilate(cor_ap, strel('disk',dilation')),strel('disk',dilation'));
        [cor_ap]=find_circles(cor_ap, bead_size_AP-bead_margin_AP, bead_size_AP+bead_margin_AP);

        region=1;
        [dis_x(1), dis_y(1), center_xr(1), center_yr(1), center_xa(1), center_ya(1), no_bead(1), n(1)]=find_center(cor_red, cor_ap, max_distance_between_beeads, offset_x, offset_y, region);


    %
    %   upper right corner
    %
        offset_x=0;
        offset_y=round(size(red1,2)/4*3);
        cor_red=red1(1:round(size(red1,1)/2),offset_y+1:end);
        cor_red=(cor_red>threshold(1));
        [cor_red]=find_circles(cor_red, bead_size_red-bead_margin_red, bead_size_red+bead_margin_red);

        cor_ap=AP(1:round(size(red1,1)/2),offset_y+1:end);
        cor_ap=(cor_ap>threshold(2));

        cor_ap=imerode(imdilate(cor_ap, strel('disk',dilation')),strel('disk',dilation'));
        [cor_ap]=find_circles(cor_ap, bead_size_AP-bead_margin_AP, bead_size_AP+bead_margin_AP);

        region=2;
        [dis_x(2), dis_y(2), center_xr(2), center_yr(2), center_xa(2), center_ya(2), no_bead(2), n(2)]=find_center(cor_red, cor_ap, max_distance_between_beeads, offset_x, offset_y, region);

    %
    %   bottom left corner
    %
        offset_x=round(size(red1,1)/2);
        offset_y=0;
        cor_red=red1(offset_x+1:end,1:round(size(red1,2)/4));
        cor_red=(cor_red>threshold(1));
        [cor_red]=find_circles(cor_red, bead_size_red-bead_margin_red, bead_size_red+bead_margin_red);

        cor_ap=AP(offset_x+1:end,1:round(size(red1,2)/4));
        cor_ap=(cor_ap>threshold(2));
        cor_ap=imerode(imdilate(cor_ap, strel('disk',dilation')),strel('disk',dilation'));
        [cor_ap]=find_circles(cor_ap, bead_size_AP-bead_margin_AP, bead_size_AP+bead_margin_AP);

        region=3;
        [dis_x(3), dis_y(3), center_xr(3), center_yr(3), center_xa(3), center_ya(3), no_bead(3), n(3)]=find_center(cor_red, cor_ap, max_distance_between_beeads, offset_x, offset_y, region);

    %
    %   bottom right corner
    %
        offset_x=round(size(red1,1)/2);
        offset_y=round(size(red1,2)/4*3);
        cor_red=red1(offset_x+1:end,offset_y+1:end);
        cor_red=(cor_red>threshold(1));
        [cor_red]=find_circles(cor_red, bead_size_red-bead_margin_red, bead_size_red+bead_margin_red);

        cor_ap=AP(offset_x+1:end,offset_y+1:end);
        cor_ap=(cor_ap>threshold(2));
        cor_ap=imerode(imdilate(cor_ap, strel('disk',dilation')),strel('disk',dilation'));
        [cor_ap]=find_circles(cor_ap, bead_size_AP-bead_margin_AP, bead_size_AP+bead_margin_AP);

        region=4;
        [dis_x(4), dis_y(4), center_xr(4), center_yr(4), center_xa(4), center_ya(4), no_bead(4), n(4)]=find_center(cor_red, cor_ap, max_distance_between_beeads, offset_x, offset_y, region);


        len=length(find(n));
        if len<2
            dis_x1=0; dis_y1=0; ratio_x=1; ratio_y=1; angle_r_minus_a=90; no_bead1=1;
            return
        end
        no_bead1=0;

        combi=zeros(1,6);
        if n(1)~=0 & n(4)~=0 & abs(dis_x(1)-dis_x(4))<100 & abs(dis_y(1)-dis_y(4))<100
            combi(1)=3*(n(1)+n(4));
        elseif n(2)~=0 & n(3)~=0 & abs(dis_x(2)-dis_x(3))<100 & abs(dis_y(2)-dis_y(3))<100
            combi(2)=3*(n(2)+n(3));
        elseif n(1)~=0 & n(3)~=0 & abs(dis_x(1)-dis_x(3))<100 & abs(dis_y(1)-dis_y(3))<200
            combi(3)=2*(n(1)+n(3));
        elseif n(1)~=0 & n(2)~=0 & abs(dis_x(1)-dis_x(2))<100 & abs(dis_y(1)-dis_y(2))<100
            combi(4)=n(1)+n(2);
        elseif n(2)~=0 & n(4)~=0 & abs(dis_x(2)-dis_x(4))<100 & abs(dis_y(2)-dis_y(4))<100
            combi(5)=2*(n(2)+n(4));
        elseif (n(5)~=0 & n(3)~=0 & abs(dis_x(5)-dis_x(3))<100 & abs(dis_y(5)-dis_y(3))<100) | (n(5)~=0 & n(4)~=0 & abs(dis_x(5)-dis_x(4))<100 & abs(dis_y(5)-dis_y(4))<100)
            if n(3)>=n(4)
                combi(7)=2*(n(5)+n(3))
            else
                combi(8)=2*(n(5)+n(4))
            end
        elseif (n(6)~=0 & n(1)~=0 & abs(dis_x(6)-dis_x(1))<100 & abs(dis_y(6)-dis_y(1))<100) | (n(6)~=0 & n(2)~=0 & abs(dis_x(6)-dis_x(2))<100 & abs(dis_y(6)-dis_y(2))<100)
            if n(1)>=n(2)
                combi(9)=2*(n(6)+n(1))
            else
                combi(10)=2*(n(6)+n(2))
            end
        elseif n(3)~=0 & n(4)~=0 & abs(dis_x(3)-dis_x(4))<100 & abs(dis_y(3)-dis_y(4))<100
            combi(6)=n(3)+n(4);
        else
            d_x5=dis_x(1:4)-dis_x(5);
            d_y5=dis_y(1:4)-dis_y(5);
            d_x6=dis_x(1:4)-dis_x(6);
            d_y6=dis_y(1:4)-dis_y(6);
            dd_x5=find(abs(d_x5)<100);
            dd_x6=find(abs(d_x6)<100);
            dd_y5=find(abs(d_y5)<100);
            dd_y6=find(abs(d_y6)<100);

            if ~isempty(dd_x5) & ~isempty(dd_y5)
                count=1;
                for loop=1:length(dd_x5)
                    for loop1=1:length(dd_y5)
                        if dd_x5(loop)==dd_y5(loop1)
                            candidate5(count)=dd_x5(loop);
                            count=count+1;
                        end
                    end
                end
            end

            if ~isempty(dd_x6) & ~isempty(dd_y6)
                count=1;
                for loop=1:length(dd_x6)
                    for loop1=1:length(dd_y6)
                        if dd_x6(loop)==dd_y6(loop1)
                            candidate6(count)=dd_x6(loop);
                            count=count+1;
                        end
                    end
                end
            end

            who5=who('candidate5');
            who6=who('candidate6');

            if isempty(who5)
                candidate5=[];
            end
            if isempty(who6)
                candidate6=[];
            end        

            if isempty(who5) & isempty(who6)
                dis_x1=0; dis_y1=0; ratio_x=1; ratio_y=1; angle_r_minus_a=90; no_bead1=1;
                return
            else
%                 test=[abs(d_x5(candidate5))+abs(d_y5(candidate5)),abs(d_x6(candidate6))+abs(d_y6(candidate6))];
%                 candidate=[candidate5,candidate6];
%                 final_candidate=candidate(find(test==min(test)));
                test=[abs(d_x5(candidate5))+abs(d_y5(candidate5)),abs(d_x6(candidate6))+abs(d_y6(candidate6))];
                test1=test(find(test>0));
                ind=min(find(test1==min(test1)));
                ind1=find(test==test1(ind(1)));
                candidate=[candidate5,candidate6];
                final_candidate=candidate(ind1(1));
            end
            if ~isempty(find(dd_x5==final_candidate))
                ref_region=5;
                opp_region=final_candidate;
            elseif ~isempty(find(dd_x6==final_candidate))
                ref_region=6;
                opp_region=final_candidate;
            end
            ind_special=1;
        end


        if ind_special==0
            ind=find(combi==max(combi));

            switch ind(1)
                case 1
                    ref_region=1;
                    opp_region=4;
                case 2
                    ref_region=2;
                    opp_region=3;
                case 3
                    ref_region=1;
                    opp_region=3;
                case 4
                    ref_region=1;
                    opp_region=2;
                case 5
                    ref_region=2;
                    opp_region=4;
                case 6
                    ref_region=3;
                    opp_region=4;
                case 7
                    ref_region=5;
                    opp_region=3;
                case 8
                    ref_region=5;
                    opp_region=4;
                case 9
                    ref_region=6;
                    opp_region=1;
                case 106
                    ref_region=6;
                    opp_region=2;
            end
        end
    end

    % nn=sort(n);
    % ref_region=find(n==nn(4));
    % %ref_region=max(ref_region);
    % ref_region=ref_region(1);
    % n1=n;
    % n1(ref_region)=0;
    % if n1(5-ref_region)~=0
    %     opp_region=5-ref_region;
    % else
    %     test=find(n1);
    %     test1=abs(test-ref_region);
    %     ind=max(test1);
    %     opp_region=test(ind);
    % end
    % opp_region=opp_region(1);
    % 
    % % n1=n;
    % % 
    % % count=1;
    % % i=1;
    % % while i<length(find(n))
    % %     temp1=find(n1==nn(4-i));
    % %     for j=1:length(temp1)
    % %         temp(count)=temp1(j);
    % %         di(count)=(dis_x(ref_region)-dis_x(temp(count)))+(dis_y(ref_region)-dis_y(temp(count)));
    % %         count=count+1;
    % %     end
    % %     i=i+j;
    % % end
    % % if len==4
    % %     comp=find(temp==5-ref_region);
    % %     if  di(comp)<mean(di)+2*abs(std(di)) | di(comp)>mean(di)-2*abs(std(di)) %abs(di(comp))<abs(std(di))%100%mean(abs(di))
    % %         opp_region=5-ref_region;
    % %     else
    % %         n1(5-ref_region)=0;
    % %         nn=sort(n1);
    % %         len=length(find(nn))+1;
    % %         clear di temp
    % %         count=1;
    % %         i=1;
    % %         while i<length(find(n1))+1
    % %             temp1=find(n1==nn(5-i));
    % %             for j=1:length(temp1)
    % %                 temp(count)=temp1(j);
    % %                 di(count)=(dis_x(ref_region)-dis_x(temp(count)))+(dis_y(ref_region)-dis_y(temp(count)));
    % %                 count=count+1;
    % %             end
    % %             i=i+j;
    % %         end
    % %     end
    % % end
    % % if len==3
    % %     if n1(5-ref_region)~=0;
    % %         opp_region=5-ref_region;
    % %     else
    % %         opp_region=find(n1==nn(4));
    % %     end
    % %     comp=find(temp==opp_region);
    % %     if di(comp)>mean(di)+2*abs(std(di)) | di(comp)<mean(di)-2*abs(std(di))%100%abs((di(3-comp)))
    % %         n1(opp_region)=0;
    % %         nn=sort(n1);
    % %         len=length(find(nn))+1;
    % %     end
    % % end
    % % if len==2
    % %     opp_region=find(n1==nn(4));    
    % % end
    % % opp_region=opp_region(1);

    if center_xr(opp_region)>center_xr(ref_region)
        center_x_denom=center_xr(opp_region);
        dis_xx=dis_x(opp_region)-dis_x(ref_region);
        dis_x1=dis_x(ref_region);
    else
        center_x_denom=center_xr(ref_region);
        dis_xx=dis_x(ref_region)-dis_x(opp_region);
        dis_x1=dis_x(opp_region);
    end
    if center_yr(opp_region)>center_yr(ref_region)
        center_y_denom=center_yr(opp_region);
        dis_yy=dis_y(opp_region)-dis_y(ref_region);
        dis_y1=dis_y(ref_region);
    else
        center_y_denom=center_yr(ref_region);
        dis_yy=dis_y(ref_region)-dis_y(opp_region);
        dis_y1=dis_y(opp_region);
    end
    
    
% elseif rot_angle==0     % left, right, middle left, middel right
elseif strcmp(bone_type,'C')==1     % left, right, middle left, middel right

    %
    %   left
    %

    %cor_red=red1(1:round(size(red1,1)/4*3),1:round(size(red1,2)/4));
    offset_x=0;
    offset_y=0;
    cor_red=red(1:end,offset_y+1:round(size(red,2)/10));
    cor_red=(cor_red>threshold(1));
    [cor_red]=find_circles(cor_red, bead_size_red-bead_margin_red, bead_size_red+bead_margin_red);

    cor_ap=AP(1:end,offset_y+1:round(size(red,2)/10));
    cor_ap=(cor_ap>threshold(2));
    %cor_ap=imerode(imdilate(cor_ap, strel('disk',dilation')),strel('disk',dilation'));
    [cor_ap]=find_circles(cor_ap, bead_size_AP-bead_margin_AP, bead_size_AP+bead_margin_AP);

    region=5;
    [dis_x(5), dis_y(5), center_xr(5), center_yr(5), center_xa(5), center_ya(5), no_bead(5), n(5)]=find_center(cor_red, cor_ap, max_distance_between_beeads, offset_x, offset_y, region);

%    if n(5)~=0
    %
    %   right
    %
        offset_x=0;
        offset_y=round(size(red,2)*9/10);
        cor_red=red(offset_x+1:end,offset_y+1:end);
        cor_red=(cor_red>threshold(1));
        [cor_red]=find_circles(cor_red, bead_size_red-bead_margin_red, bead_size_red+bead_margin_red);

        cor_ap=AP(offset_x+1:end,offset_y+1:end);
        cor_ap=(cor_ap>threshold(2));
        %cor_ap=imerode(imdilate(cor_ap, strel('disk',dilation')),strel('disk',dilation'));
        [cor_ap]=find_circles(cor_ap, bead_size_AP-bead_margin_AP, bead_size_AP+bead_margin_AP);

        region=6;
        [dis_x(6), dis_y(6), center_xr(6), center_yr(6), center_xa(6), center_ya(6), no_bead(6), n(6)]=find_center(cor_red, cor_ap, max_distance_between_beeads, offset_x, offset_y, region);

%    end

    if n(5)~=0 & n(6)~=0 & abs(dis_x(6)-dis_x(5))<100 & abs(dis_y(6)-dis_y(5))<100
        ref_region=5;
        opp_region=6;
        no_bead1=0;
    else
        a=(DIC>graythresh(DIC)*255);
        a=imerode(imdilate(a, strel('disk',31')),strel('disk',31'));
        [L n1]=bwlabel(a); clear a
        stats=regionprops(L, 'Area','BoundingBox','FIlledImage'); clear L
        area_dic=[stats.Area];
        [I J]=find(area_dic>200);
        temp=zeros(size(DIC));
        for i=1:size(J,2)
            sx=round(stats(J(i)).BoundingBox(2));
            sy=round(stats(J(i)).BoundingBox(1));
            ex=sx+stats(J(i)).BoundingBox(4)-1;
            ey=sy+stats(J(i)).BoundingBox(3)-1;
            temp(sx:ex,sy:ey)=temp(sx:ex,sy:ey)+stats(J(i)).FilledImage;
        end
        temp=1-temp;
        temp=uint8(temp>0);
        red1=immultiply(red,temp); clear temp



        red1=red;
    %
    %   middle left
    %

    %cor_red=red1(1:round(size(red1,1)/4*3),1:round(size(red1,2)/4));
        offset_x=round(size(red1,1)/4);
        offset_y=round(size(red1,2)/4);

        cor_red=red1(offset_x+1:round(size(red1,1)*3/4),offset_y+1:round(size(red1,2)/2));
        cor_red=(cor_red>threshold(1));
        [cor_red]=find_circles(cor_red, bead_size_red-bead_margin_red, bead_size_red+bead_margin_red);

        cor_ap=AP(offset_x+1:round(size(red1,1)*3/4),offset_y+1:round(size(red1,2)/2));
        cor_ap=(cor_ap>threshold(2));
        cor_ap=imerode(imdilate(cor_ap, strel('disk',dilation')),strel('disk',dilation'));
        [cor_ap]=find_circles(cor_ap, bead_size_AP-bead_margin_AP, bead_size_AP+bead_margin_AP);

        region=1;
        [dis_x(1), dis_y(1), center_xr(1), center_yr(1), center_xa(1), center_ya(1), no_bead(1), n(1)]=find_center(cor_red, cor_ap, max_distance_between_beeads, offset_x, offset_y, region);


    %
    %   middle right
    %
        offset_x=round(size(red1,1)/4);
        offset_y=round(size(red1,2)/2);
        cor_red=red1(offset_x+1:round(size(red1,1)*3/4),offset_y+1:round(size(red1,2)*3/4));
        cor_red=(cor_red>threshold(1));
        [cor_red]=find_circles(cor_red, bead_size_red-bead_margin_red, bead_size_red+bead_margin_red);

        cor_ap=AP(offset_x+1:round(size(red1,1)*3/4),offset_y+1:round(size(red1,2)*3/4));
        cor_ap=(cor_ap>threshold(2));

        cor_ap=imerode(imdilate(cor_ap, strel('disk',dilation')),strel('disk',dilation'));
        [cor_ap]=find_circles(cor_ap, bead_size_AP-bead_margin_AP, bead_size_AP+bead_margin_AP);

        region=2;
        [dis_x(2), dis_y(2), center_xr(2), center_yr(2), center_xa(2), center_ya(2), no_bead(2), n(2)]=find_center(cor_red, cor_ap, max_distance_between_beeads, offset_x, offset_y, region);

%     %
%     %   bottom left corner
%     %
%         offset_x=round(size(red1,1)/2);
%         offset_y=0;
%         cor_red=red1(offset_x+1:end,1:round(size(red1,2)/4));
%         cor_red=(cor_red>threshold(1));
%         [cor_red]=find_circles(cor_red, bead_size_red-bead_margin_red, bead_size_red+bead_margin_red);
% 
%         cor_ap=AP(offset_x+1:end,1:round(size(red1,2)/4));
%         cor_ap=(cor_ap>threshold(2));
%         cor_ap=imerode(imdilate(cor_ap, strel('disk',dilation')),strel('disk',dilation'));
%         [cor_ap]=find_circles(cor_ap, bead_size_AP-bead_margin_AP, bead_size_AP+bead_margin_AP);
% 
%         region=3;
%         [dis_x(3), dis_y(3), center_xr(3), center_yr(3), center_xa(3), center_ya(3), no_bead(3), n(3)]=find_center(cor_red, cor_ap, max_distance_between_beeads, offset_x, offset_y, region);
% 
%     %
%     %   bottom right corner
%     %
%         offset_x=round(size(red1,1)/2);
%         offset_y=round(size(red1,2)/4*3);
%         cor_red=red1(offset_x+1:end,offset_y+1:end);
%         cor_red=(cor_red>threshold(1));
%         [cor_red]=find_circles(cor_red, bead_size_red-bead_margin_red, bead_size_red+bead_margin_red);
% 
%         cor_ap=AP(offset_x+1:end,offset_y+1:end);
%         cor_ap=(cor_ap>threshold(2));
%         cor_ap=imerode(imdilate(cor_ap, strel('disk',dilation')),strel('disk',dilation'));
%         [cor_ap]=find_circles(cor_ap, bead_size_AP-bead_margin_AP, bead_size_AP+bead_margin_AP);
% 
%         region=4;
%         [dis_x(4), dis_y(4), center_xr(4), center_yr(4), center_xa(4), center_ya(4), no_bead(4), n(4)]=find_center(cor_red, cor_ap, max_distance_between_beeads, offset_x, offset_y, region);


        len=length(find(n));
        if len<2
            dis_x1=0; dis_y1=0; ratio_x=1; ratio_y=1; angle_r_minus_a=90; no_bead1=1;
            return
        end
        no_bead1=0;

        combi=zeros(1,6);
        if n(1)~=0 & n(4)~=0 & abs(dis_x(1)-dis_x(4))<100 & abs(dis_y(1)-dis_y(4))<100
            combi(1)=3*(n(1)+n(4));
        elseif n(2)~=0 & n(3)~=0 & abs(dis_x(2)-dis_x(3))<100 & abs(dis_y(2)-dis_y(3))<100
            combi(2)=3*(n(2)+n(3));
        elseif n(1)~=0 & n(3)~=0 & abs(dis_x(1)-dis_x(3))<100 & abs(dis_y(1)-dis_y(3))<100
            combi(3)=2*(n(1)+n(3));
        elseif n(1)~=0 & n(2)~=0 & abs(dis_x(1)-dis_x(2))<100 & abs(dis_y(1)-dis_y(2))<100
            combi(4)=n(1)+n(2);
        elseif n(2)~=0 & n(4)~=0 & abs(dis_x(2)-dis_x(4))<100 & abs(dis_y(2)-dis_y(4))<100
            combi(5)=2*(n(2)+n(4));
        elseif n(3)~=0 & n(4)~=0 & abs(dis_x(3)-dis_x(4))<100 & abs(dis_y(3)-dis_y(4))<100
            combi(6)=n(3)+n(4);
        else
            d_x5=dis_x(1:4)-dis_x(5);
            d_y5=dis_y(1:4)-dis_y(5);
            d_x6=dis_x(1:4)-dis_x(6);
            d_y6=dis_y(1:4)-dis_y(6);
            dd_x5=find(abs(d_x5)<100);
            dd_x6=find(abs(d_x6)<100);
            dd_y5=find(abs(d_y5)<100);
            dd_y6=find(abs(d_y6)<100);

            if ~isempty(dd_x5) & ~isempty(dd_y5)
                count=1;
                for loop=1:length(dd_x5)
                    for loop1=1:length(dd_y5)
                        if dd_x5(loop)==dd_y5(loop1)
                            candidate5(count)=dd_x5(loop);
                            count=count+1;
                        end
                    end
                end
            end

            if ~isempty(dd_x6) & ~isempty(dd_y6)
                count=1;
                for loop=1:length(dd_x6)
                    for loop1=1:length(dd_y6)
                        if dd_x6(loop)==dd_y6(loop1)
                            candidate6(count)=dd_x6(loop);
                            count=count+1;
                        end
                    end
                end
            end

            who5=who('candidate5');
            who6=who('candidate6');

            if isempty(who5)
                candidate5=[];
            end
            if isempty(who6)
                candidate6=[];
            end        

            if isempty(who5) & isempty(who6)
                dis_x1=0; dis_y1=0; ratio_x=1; ratio_y=1; angle_r_minus_a=90; no_bead1=1;
                return
            else
%                 test=[abs(d_x5(candidate5))+abs(d_y5(candidate5)),abs(d_x6(candidate6))+abs(d_y6(candidate6))];
%                 candidate=[candidate5,candidate6];
%                 final_candidate=candidate(find(test==min(test)));
                test=[abs(d_x5(candidate5))+abs(d_y5(candidate5)),abs(d_x6(candidate6))+abs(d_y6(candidate6))];
                test1=test(find(test>0));
                ind=min(find(test1==min(test1)));
                ind1=find(test==test1(ind));
                candidate=[candidate5,candidate6];
                final_candidate=candidate(ind1);
            end
            if ~isempty(find(dd_x5==final_candidate))
                ref_region=5;
                opp_region=final_candidate;
            elseif ~isempty(find(dd_x6==final_candidate))
                ref_region=6;
                opp_region=final_candidate;
            end
            ind_special=1;
        end


        if ind_special==0
            ind=find(combi==max(combi));

            switch ind(1)
                case 1
                    ref_region=1;
                    opp_region=4;
                case 2
                    ref_region=2;
                    opp_region=3;
                case 3
                    ref_region=1;
                    opp_region=3;
                case 4
                    ref_region=1;
                    opp_region=2;
                case 5
                    ref_region=2;
                    opp_region=4;
                case 6
                    ref_region=3;
                    opp_region=4;
            end
        end
    end

    if center_xr(opp_region)>center_xr(ref_region)
        center_x_denom=center_xr(opp_region);
        dis_xx=dis_x(opp_region)-dis_x(ref_region);
        dis_x1=dis_x(ref_region);
    else
        center_x_denom=center_xr(ref_region);
        dis_xx=dis_x(ref_region)-dis_x(opp_region);
        dis_x1=dis_x(opp_region);
    end
    if center_yr(opp_region)>center_yr(ref_region)
        center_y_denom=center_yr(opp_region);
        dis_yy=dis_y(opp_region)-dis_y(ref_region);
        dis_y1=dis_y(ref_region);
    else
        center_y_denom=center_yr(ref_region);
        dis_yy=dis_y(ref_region)-dis_y(opp_region);
        dis_y1=dis_y(opp_region);
    end    
elseif strcmp(bone_type,'L')==1 | strcmp(bone_type,'R')==1    % left, right, middle left, middel right

    %
    %   left
    %

    %cor_red=red1(1:round(size(red1,1)/4*3),1:round(size(red1,2)/4));
    offset_x=0;
    offset_y=0;
    cor_red=red(1:end,offset_y+1:round(size(red,2)/2));
    cor_red=(cor_red>threshold(1));
    [cor_red]=find_circles(cor_red, bead_size_red-bead_margin_red, bead_size_red+bead_margin_red);

    cor_ap=AP(1:end,offset_y+1:round(size(red,2)/2));
    cor_ap=(cor_ap>threshold(2));
    %cor_ap=imerode(imdilate(cor_ap, strel('disk',dilation')),strel('disk',dilation'));
    [cor_ap]=find_circles(cor_ap, bead_size_AP-bead_margin_AP, bead_size_AP+bead_margin_AP);

    region=5;
    [dis_x(5), dis_y(5), center_xr(5), center_yr(5), center_xa(5), center_ya(5), no_bead(5), n(5)]=find_center(cor_red, cor_ap, max_distance_between_beeads, offset_x, offset_y, region);

%    if n(5)~=0
    %
    %   right
    %
        offset_x=0;
        offset_y=round(size(red,2)/2);
        cor_red=red(offset_x+1:end,offset_y+1:end);
        cor_red=(cor_red>threshold(1));
        [cor_red]=find_circles(cor_red, bead_size_red-bead_margin_red, bead_size_red+bead_margin_red);

        cor_ap=AP(offset_x+1:end,offset_y+1:end);
        cor_ap=(cor_ap>threshold(2));
        %cor_ap=imerode(imdilate(cor_ap, strel('disk',dilation')),strel('disk',dilation'));
        [cor_ap]=find_circles(cor_ap, bead_size_AP-bead_margin_AP, bead_size_AP+bead_margin_AP);

        region=6;
        [dis_x(6), dis_y(6), center_xr(6), center_yr(6), center_xa(6), center_ya(6), no_bead(6), n(6)]=find_center(cor_red, cor_ap, max_distance_between_beeads, offset_x, offset_y, region);

%    end

    if n(5)~=0 & n(6)~=0 & abs(dis_x(6)-dis_x(5))<100 & abs(dis_y(6)-dis_y(5))<100
        ref_region=5;
        opp_region=6;
        no_bead1=0;
    else
        dis_x1=0;
        dis_y1=0;
        ratio_x=0;
        ratio_y=0;
        angle_r_minus_a=0;
        no_bead1=1;
        return
%         a=(DIC>graythresh(DIC)*255);
%         a=imerode(imdilate(a, strel('disk',31')),strel('disk',31'));
%         [L n1]=bwlabel(a); clear a
%         stats=regionprops(L, 'Area','BoundingBox','FIlledImage'); clear L
%         area_dic=[stats.Area];
%         [I J]=find(area_dic>200);
%         temp=zeros(size(DIC));
%         for i=1:size(J,2)
%             sx=round(stats(J(i)).BoundingBox(2));
%             sy=round(stats(J(i)).BoundingBox(1));
%             ex=sx+stats(J(i)).BoundingBox(4)-1;
%             ey=sy+stats(J(i)).BoundingBox(3)-1;
%             temp(sx:ex,sy:ey)=temp(sx:ex,sy:ey)+stats(J(i)).FilledImage;
%         end
%         temp=1-temp;
%         temp=uint8(temp>0);
%         red1=immultiply(red,temp); clear temp



%         red1=red;
%     %
%     %   middle left
%     %
% 
%     %cor_red=red1(1:round(size(red1,1)/4*3),1:round(size(red1,2)/4));
%         offset_x=round(size(red1,1)/4);
%         offset_y=round(size(red1,2)/4);
% 
%         cor_red=red1(offset_x+1:round(size(red1,1)*3/4),offset_y+1:round(size(red1,2)/2));
%         cor_red=(cor_red>threshold(1));
%         [cor_red]=find_circles(cor_red, bead_size_red-bead_margin_red, bead_size_red+bead_margin_red);
% 
%         cor_ap=AP(offset_x+1:round(size(red1,1)*3/4),offset_y+1:round(size(red1,2)/2));
%         cor_ap=(cor_ap>threshold(2));
%         cor_ap=imerode(imdilate(cor_ap, strel('disk',dilation')),strel('disk',dilation'));
%         [cor_ap]=find_circles(cor_ap, bead_size_AP-bead_margin_AP, bead_size_AP+bead_margin_AP);
% 
%         region=1;
%         [dis_x(1), dis_y(1), center_xr(1), center_yr(1), center_xa(1), center_ya(1), no_bead(1), n(1)]=find_center(cor_red, cor_ap, max_distance_between_beeads, offset_x, offset_y, region);
% 
% 
%     %
%     %   middle right
%     %
%         offset_x=round(size(red1,1)/4);
%         offset_y=round(size(red1,2)/2);
%         cor_red=red1(offset_x+1:round(size(red1,1)*3/4),offset_y+1:round(size(red1,2)*3/4));
%         cor_red=(cor_red>threshold(1));
%         [cor_red]=find_circles(cor_red, bead_size_red-bead_margin_red, bead_size_red+bead_margin_red);
% 
%         cor_ap=AP(offset_x+1:round(size(red1,1)*3/4),offset_y+1:round(size(red1,2)*3/4));
%         cor_ap=(cor_ap>threshold(2));
% 
%         cor_ap=imerode(imdilate(cor_ap, strel('disk',dilation')),strel('disk',dilation'));
%         [cor_ap]=find_circles(cor_ap, bead_size_AP-bead_margin_AP, bead_size_AP+bead_margin_AP);
% 
%         region=2;
%         [dis_x(2), dis_y(2), center_xr(2), center_yr(2), center_xa(2), center_ya(2), no_bead(2), n(2)]=find_center(cor_red, cor_ap, max_distance_between_beeads, offset_x, offset_y, region);

%     %
%     %   bottom left corner
%     %
%         offset_x=round(size(red1,1)/2);
%         offset_y=0;
%         cor_red=red1(offset_x+1:end,1:round(size(red1,2)/4));
%         cor_red=(cor_red>threshold(1));
%         [cor_red]=find_circles(cor_red, bead_size_red-bead_margin_red, bead_size_red+bead_margin_red);
% 
%         cor_ap=AP(offset_x+1:end,1:round(size(red1,2)/4));
%         cor_ap=(cor_ap>threshold(2));
%         cor_ap=imerode(imdilate(cor_ap, strel('disk',dilation')),strel('disk',dilation'));
%         [cor_ap]=find_circles(cor_ap, bead_size_AP-bead_margin_AP, bead_size_AP+bead_margin_AP);
% 
%         region=3;
%         [dis_x(3), dis_y(3), center_xr(3), center_yr(3), center_xa(3), center_ya(3), no_bead(3), n(3)]=find_center(cor_red, cor_ap, max_distance_between_beeads, offset_x, offset_y, region);
% 
%     %
%     %   bottom right corner
%     %
%         offset_x=round(size(red1,1)/2);
%         offset_y=round(size(red1,2)/4*3);
%         cor_red=red1(offset_x+1:end,offset_y+1:end);
%         cor_red=(cor_red>threshold(1));
%         [cor_red]=find_circles(cor_red, bead_size_red-bead_margin_red, bead_size_red+bead_margin_red);
% 
%         cor_ap=AP(offset_x+1:end,offset_y+1:end);
%         cor_ap=(cor_ap>threshold(2));
%         cor_ap=imerode(imdilate(cor_ap, strel('disk',dilation')),strel('disk',dilation'));
%         [cor_ap]=find_circles(cor_ap, bead_size_AP-bead_margin_AP, bead_size_AP+bead_margin_AP);
% 
%         region=4;
%         [dis_x(4), dis_y(4), center_xr(4), center_yr(4), center_xa(4), center_ya(4), no_bead(4), n(4)]=find_center(cor_red, cor_ap, max_distance_between_beeads, offset_x, offset_y, region);


        len=length(find(n));
        if len<2
            dis_x1=0; dis_y1=0; ratio_x=1; ratio_y=1; angle_r_minus_a=90; no_bead1=1;
            return
        end
        no_bead1=0;

        combi=zeros(1,6);
        if n(1)~=0 & n(4)~=0 & abs(dis_x(1)-dis_x(4))<100 & abs(dis_y(1)-dis_y(4))<100
            combi(1)=3*(n(1)+n(4));
        elseif n(2)~=0 & n(3)~=0 & abs(dis_x(2)-dis_x(3))<100 & abs(dis_y(2)-dis_y(3))<100
            combi(2)=3*(n(2)+n(3));
        elseif n(1)~=0 & n(3)~=0 & abs(dis_x(1)-dis_x(3))<100 & abs(dis_y(1)-dis_y(3))<100
            combi(3)=2*(n(1)+n(3));
        elseif n(1)~=0 & n(2)~=0 & abs(dis_x(1)-dis_x(2))<100 & abs(dis_y(1)-dis_y(2))<100
            combi(4)=n(1)+n(2);
        elseif n(2)~=0 & n(4)~=0 & abs(dis_x(2)-dis_x(4))<100 & abs(dis_y(2)-dis_y(4))<100
            combi(5)=2*(n(2)+n(4));
        elseif n(3)~=0 & n(4)~=0 & abs(dis_x(3)-dis_x(4))<100 & abs(dis_y(3)-dis_y(4))<100
            combi(6)=n(3)+n(4);
        else
            d_x5=dis_x(1:4)-dis_x(5);
            d_y5=dis_y(1:4)-dis_y(5);
            d_x6=dis_x(1:4)-dis_x(6);
            d_y6=dis_y(1:4)-dis_y(6);
            dd_x5=find(abs(d_x5)<100);
            dd_x6=find(abs(d_x6)<100);
            dd_y5=find(abs(d_y5)<100);
            dd_y6=find(abs(d_y6)<100);

            if ~isempty(dd_x5) & ~isempty(dd_y5)
                count=1;
                for loop=1:length(dd_x5)
                    for loop1=1:length(dd_y5)
                        if dd_x5(loop)==dd_y5(loop1)
                            candidate5(count)=dd_x5(loop);
                            count=count+1;
                        end
                    end
                end
            end

            if ~isempty(dd_x6) & ~isempty(dd_y6)
                count=1;
                for loop=1:length(dd_x6)
                    for loop1=1:length(dd_y6)
                        if dd_x6(loop)==dd_y6(loop1)
                            candidate6(count)=dd_x6(loop);
                            count=count+1;
                        end
                    end
                end
            end

            who5=who('candidate5');
            who6=who('candidate6');

            if isempty(who5)
                candidate5=[];
            end
            if isempty(who6)
                candidate6=[];
            end        

            if isempty(who5) & isempty(who6)
                dis_x1=0; dis_y1=0; ratio_x=1; ratio_y=1; angle_r_minus_a=90; no_bead1=1;
                return
            else
%                 test=[abs(d_x5(candidate5))+abs(d_y5(candidate5)),abs(d_x6(candidate6))+abs(d_y6(candidate6))];
%                 candidate=[candidate5,candidate6];
%                 final_candidate=candidate(find(test==min(test)));
                test=[abs(d_x5(candidate5))+abs(d_y5(candidate5)),abs(d_x6(candidate6))+abs(d_y6(candidate6))];
                test1=test(find(test>0));
                ind=min(find(test1==min(test1)));
                ind1=find(test==test1(ind));
                candidate=[candidate5,candidate6];
                final_candidate=candidate(ind1);
            end
            if ~isempty(find(dd_x5==final_candidate))
                ref_region=5;
                opp_region=final_candidate;
            elseif ~isempty(find(dd_x6==final_candidate))
                ref_region=6;
                opp_region=final_candidate;
            end
            ind_special=1;
        end


        if ind_special==0
            ind=find(combi==max(combi));

            switch ind(1)
                case 1
                    ref_region=1;
                    opp_region=4;
                case 2
                    ref_region=2;
                    opp_region=3;
                case 3
                    ref_region=1;
                    opp_region=3;
                case 4
                    ref_region=1;
                    opp_region=2;
                case 5
                    ref_region=2;
                    opp_region=4;
                case 6
                    ref_region=3;
                    opp_region=4;
            end
        end
    end

    if center_xr(opp_region)>center_xr(ref_region)
        center_x_denom=center_xr(opp_region);
        dis_xx=dis_x(opp_region)-dis_x(ref_region);
        dis_x1=dis_x(ref_region);
    else
        center_x_denom=center_xr(ref_region);
        dis_xx=dis_x(ref_region)-dis_x(opp_region);
        dis_x1=dis_x(opp_region);
    end
    if center_yr(opp_region)>center_yr(ref_region)
        center_y_denom=center_yr(opp_region);
        dis_yy=dis_y(opp_region)-dis_y(ref_region);
        dis_y1=dis_y(ref_region);
    else
        center_y_denom=center_yr(ref_region);
        dis_yy=dis_y(ref_region)-dis_y(opp_region);
        dis_y1=dis_y(opp_region);
    end    
end
ratio_x=1+dis_xx(1)/center_x_denom;
ratio_y=1+dis_yy(1)/center_y_denom;

% angle_a=atan((center_xa(opp_region)-center_xa(ref_region))/(center_ya(opp_region)-center_ya(ref_region)));
% angle_r=atan((center_xr(opp_region)-center_xr(ref_region))/(center_yr(opp_region)-center_yr(ref_region)));
% angle_r_minus_a=(angle_r-angle_a)/pi*180;
angle_a=atan((center_ya(opp_region)-center_ya(ref_region))/(center_xa(opp_region)-center_xa(ref_region)));
angle_r=atan((center_yr(opp_region)-center_yr(ref_region))/(center_xr(opp_region)-center_xr(ref_region)));
angle_r_minus_a=(angle_a-angle_r)/pi*180;

%dis_x1=round((dis_x(ref_region)+dis_x(opp_region))/2);
%dis_y1=round((dis_y(ref_region)+dis_y(opp_region))/2);
return

function [result]=MakeRegistrationReport(r_info)
load(r_info);
aa=reg_info(:,3);
result=['Registration is done']; 

if isempty(aa==10000)
    result=['No DIC Image'];
elseif isempty(aa==1000)
    result=['No AP Image'];
end
messages={'TRAP Registration Error : Need to check beads in DIC and TRAP images';
    'AP Registration Error : Need to check beads in DIC and AP images';
    'TB Registration Error : Need to check beads in TB images';
    'Cy5 Registration Error : Need to check beads in Cy5 images';
    'SO Registration Error : Need to check beads in SO images'};
for i=3:7
    testOK=find(aa==i);
    if isempty(testOK)
        testError=find(aa==i+100);
        if ~isempty(testError)
            result=messages{i-2,:};
            break
        end
    end
end

function [info1, info2]=read_file_info_F2(a,name)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% registration_KOMP_8_v2channel
%
% Add Cy5 SafO, and SafO
% 2020 July 4th
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

old_bone_no=0;
info1=init_info;
info2=info1;% info1={};
% info2={};

first=1;
first_filter_m=1;
first_filter_t=1;
row=0;
first_F=1;
first_M=1;
flag=1;
bn=[];
for i=3:length(a)
    if a(i).bytes<500000 | a(i).isdir==1
        continue
    end
    if i==26
        i
    end
    aa=a(i).name;
    gender=aa(regexp(aa,name{1}(1:end-1))+length(name{1})-1);
%     bone_number=aa(regexp(aa,name{1}(1:end-1))+length(name{1}):regexp(aa,name{1}(1:end-1))+length(name{1})+1);
    bone_number=str2num(aa(regexp(aa,'_[A,M,T,B,C,S]_s[1-9]c')+4));

%     section_number=aa(regexp(aa,'_[A,M,T]_s[1-9]c')+4);
    section_number=str2num(aa(regexp(aa,name{1}(1:end-1))+length(name{1})+1));
    filter_type=aa(regexp(aa,'_[A,M,T,B,C,S]_s[1-9]c')+1);
    channel=aa(regexp(aa,'_[A,M,T,B,C,S]_s[1-9]c')+6);
    colors=[filter_type,channel];

%     if ~strcmp(old_bone_no,bone_number)
%         row=row+1;
%         old_bone_no=bone_number;
%     end

%     switch gender
%         case 'F'
%             info=info1;
%         case 'M'
%             info=info2;
%     end
    
    switch colors
       case 'A1'
            switch gender
                case 'F'
                    info1.ap{bone_number,section_number+1}=1;
                case 'M'
                    info2.ap{bone_number,section_number+1}=1;
            end                        
%        case 'A3'
%            switch gender
%                case 'M'
%                    info2.ap{bone_number,section_number+1}=1;
%            end
%        case 'A2'
%            switch gender
%                case 'F'
%                    info1.ap{bone_number,section_number+1}=1;
%            end


        case 'M1'
            switch gender
                case 'F'
                    info1.im{bone_number,section_number+1}=1;
                case 'M'
                    info2.im{bone_number,section_number+1}=1;
            end 
        case 'T2'
            switch gender
                case 'F'
                    info1.tr{bone_number,section_number+1}=1;
                case 'M'
                    info2.tr{bone_number,section_number+1}=1;
            end
        case 'B2'
            switch gender
                case 'F'
                    info1.tb{bone_number,section_number+1}=1;
                case 'M'
                    info2.tb{bone_number,section_number+1}=1;
            end
        case 'C2'
            switch gender
                case 'F'
                    info1.cy{bone_number,section_number+1}=1;
                case 'M'
                    info2.cy{bone_number,section_number+1}=1;
            end
        case 'S1'
            switch gender
                case 'F'
                    info1.so{bone_number,section_number+1}=1;
                case 'M'
                    info2.so{bone_number,section_number+1}=1;
            end
    end

end

function [info]=init_info
for i=1:9
    info1{i,1}=num2str(i);
end
info1{8,4}=[];
info.im=info1;
info.ap=info1;
info.tr=info1;
info.tb=info1;
info.cy=info1;
info.so=info1;



         
