function write_data6_cell(varargin)

% function write_data6_cell(direct, phrase1, phrase2, bone_type, delimeter)

% direct='Z:\DO\DO-AC01_G09\01_Submitted\Layers\Images\data\';
% phrase1={'DO-AC01_G09_hF_F';'DO-AC01_G09_hF_M'};
% phrase2={'Female';                              % Exp name for output excel file
%       'Male'};
% bone_type='F';                                % Femur, Vertebra
% delimeter='\';
%

% bone_number1=[792,793,794,858,859,860,861,862];                          % Control_Female
% bone_number2=[47,48,49,53,54,57,58,59];                          % Control_Male

direct=varargin{1};
phrase1=varargin{2};
phrase2=varargin{3};
bone_type=varargin{4};
delimeter=varargin{5};
n=strfind(direct,'Layers');
root_di=[direct(1:n+5),delimeter];
eval(['load ''',root_di,'info''']);
for i=1:length(phrase1)
    eval(['info=info',num2str(i),';'])
    eval(['bone_number',num2str(i),'=info',num2str(i),'.im;'])
end
% for i=5:nargin
%     eval(['bone_number',num2str(i-4),'=varargin{i};'])
% end
not_include=[];
%[9, 1;
%             22, 2;
%             23, 1];

%clear         
fclose all
section_effect=1;           % section_effect 1 : use mean value of sections of 1 mouce as 1 data
                            % section_effect 0 : use all of the values of sections of all mice
                            
dash=strfind(direct,'\');
temp_dir=direct(1:dash(end-2));
load ([temp_dir,'info.mat'])

include_number_TRAP=0;
% col_name1={'AP/BS';
%     'AP\_R/BS';
%     'AP\_only/BS';
%     'GFP/BS';
%     'GFP\_R/BS';
%     'GFP\_only/BS';
%     'TRAP/BS';
%     'TRAP\_R/BS';
%     'TRAP\_only/BS';
%     'TRAP\_on/TRAP';
%     'TRAP/TV'};
% col_name2={'AP_BS';
%     'AP_R_BS';
%     'AP_only_BS';
%     'GFP_BS';
%     'GFP_R_BS';
%     'GFP_only_BS';
%     'TRAP_BS';
%     'TRAP_R_BS';
%     'TRAP_only_BS';
%     'TRAP_on_TRAP';
%     'TRAP_TV'};
col_name1={'AP\_BS';         % 21
    'AP\_R\_BS';             % 22
    'AP\_G\_BS';             % 23
    'AP\_RG\_BS';            % 24
    'AP\_R\_RG\_BS';         % 25
    'AP\_L\_BS';             % 26
    'AP\_NL\_BS';            % 27
    'AP\_L\_AP';             % 28
    'AP\_NL\_AP';            % 29
    'AP\_R\_R';              % 30
    'AP\_G\_G';              % 31
    'AP\_RG\_RG';            % 32
    'AP\_R_RG\_R';           % 33
    'AP\_TRAP\_BS';          % 34
    'TRAP\_BS';              % 35
    'TRAP\_L\_BS';           % 36
    'TRAP\_NL\_BS';          % 37
    'AP\_TRAP\_R\_RG\_BS';   % 38
    'TRAP\_on\_TRAP';        % 39
    'TRAP\_L\_TRAP\_on';     % 40
    'TRAP\_NL\_TRAP\_on';    % 41
    'GFP\_BS';               % 42
    'GFP\_R\_BS';            % 43
    'GFPonly\_BS';           % 44
    'Height\_Width';         % 45
    'Osteocytes_Density';    % 46
    'Cortex_Width';          % 47
    'AC_Intensity';          % 48
    'RedBeads_Intensity';    % 49
    'AC_RedBeads';           % 50
    'AC_BR_Intensity';       % 51
    'AC_BR_RedBeads';        % 52
    'Calcein_Intensity';     % 53
    'GreenBeads_Intensity';  % 54
    'Calcein_GreenBeads';    % 55
    };
col_name2={'AP_BS';         % 21
    'AP_R_BS';              % 22
    'AP_G_BS';              % 23
    'AP_RG_BS';             % 24
    'AP_R_RG_BS';           % 25
    'AP_L_BS';              % 26
    'AP_NL_BS';             % 27
    'AP_L_AP';              % 28
    'AP_NL_AP';             % 29
    'AP_R_R';               % 30
    'AP_G_G';               % 31
    'AP_RG_RG';             % 32
    'AP_R_RG_R';            % 33
    'AP_TRAP_BS';           % 34
    'TRAP_BS';              % 35
    'TRAP_L_BS';            % 36
    'TRAP_NL_BS';           % 37
    'AP_TRAP_R_RG_BS';      % 38
    'TRAP_on_TRAP';         % 39
    'TRAP_L_TRAP_on';       % 40
    'TRAP_NL_TRAP_on';      % 41
    'GFP_BS';               % 42
    'GFP_R_BS';             % 43
    'GFPonly_BS';           % 44
    'Height_Width';         % 45
    'Osteocytes_Density';   % 46
    'Cortex_Width';         % 47
    'AC_Intensity';         % 48
    'RedBeads_Intensity';   % 49
    'AC_RedBeads';          % 50
    'AC_BR_Intensity';      % 51
    'AC_BR_RedBeads';       % 52
    'Calcein_Intensity';    % 53
    'GreenBeads_Intensity'; % 54
    'Calcein_GreenBeads';   % 55
    };
% if bone_type=='V'
%     ana_n=0;
% else
    ana_n=1;               % ana_n=0 : first analysis, ana_n=1 : second analysis
% end

no_column=length(col_name1);

no_roi=1;

tic

if strcmp(direct(end),delimeter)~=1
    direct=[direct,delimeter];
end

fclose('all')

aaa=phrase1{1,:};

eval(['out_file1=''',direct,aaa(1:end-1),'analysis1_cell_tr1.txt'';'])
fid_g1 = fopen(out_file1, 'w');
if no_roi==2
    eval(['out_file2=''',direct,aaa(1:end-1),'analysis1_cell_tr2.txt'';'])
    fid_g2 = fopen(out_file2, 'w');
end
if include_number_TRAP==0
    fprintf(fid_g1,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n', ...
        'Bone type','Bone #','AP_BS', 'AP_R_BS', 'AP_G_BS', 'AP_RG_BS', 'AP_R_RG_BS', 'AP_L_BS', 'AP_NL_BS', 'AP_L_AP', 'AP_NL_AP', 'AP_R_R', 'AP_G_G', 'AP_RG_RG', 'AP_R_RG_R', 'AP_TRAP_BS', 'TRAP_BS', 'TRAP_L_BS', 'TRAP_NL_BS', 'AP_TRAP_R_RG_BS', 'TRAP_on_TRAP', 'TRAP_L_TRAP_on', 'TRAP_NL_TRAP_on', 'GFP_BS', 'GFP_R_BS', 'GFPonly_BS', 'Height_Width', 'Osteocytes_N', 'Cortex_Width', 'AC_Intensity', 'RedBeads_Intensity', 'AC_RedBeads', 'AC_BR_Intensity', 'AC_BR_RedBeads', 'Calcein_Intensity', 'GreenBeads_Intensity', 'Calcein_GreenBeads');
    if no_roi==2
        fprintf(fid_g2,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n', ...
        'Bone type','Bone #','AP_BS', 'AP_R_BS', 'AP_G_BS', 'AP_RG_BS', 'AP_R_RG_BS', 'AP_L_BS', 'AP_NL_BS', 'AP_L_AP', 'AP_NL_AP', 'AP_R_R', 'AP_G_G', 'AP_RG_RG', 'AP_R_RG_R', 'AP_TRAP_BS', 'TRAP_BS', 'TRAP_L_BS', 'TRAP_NL_BS', 'AP_TRAP_R_RG_BS', 'TRAP_on_TRAP', 'TRAP_L_TRAP_on', 'TRAP_NL_TRAP_on', 'GFP_BS', 'GFP_R_BS', 'GFPonly_BS', 'Height_Width', 'Osteocytes_N', 'Cortex_Width', 'AC_Intensity', 'RedBeads_Intensity', 'AC_RedBeads', 'AC_BR_Intensity', 'AC_BR_RedBeads', 'Calcein_Intensity', 'GreenBeads_Intensity', 'Calcein_GreenBeads');
    end
elseif include_number_TRAP==1
    fprintf(fid_g1,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n', ...
        'Bone type','Bone #','AP_BS', 'AP_R_BS', 'AP_G_BS', 'AP_RG_BS', 'AP_R_RG_BS', 'AP_L_BS', 'AP_NL_BS', 'AP_L_AP', 'AP_NL_AP', 'AP_R_R', 'AP_G_G', 'AP_RG_RG', 'AP_R_RG_R', 'AP_TRAP_BS', 'TRAP_BS', 'TRAP_L_BS', 'TRAP_NL_BS', 'AP_TRAP_R_RG_BS', 'TRAP_on_TRAP', 'TRAP_L_TRAP_on', 'TRAP_NL_TRAP_on', 'GFP_BS', 'GFP_R_BS', 'GFPonly_BS', 'Height_Width', 'Osteocytes_N', 'Cortex_Width', 'AC_Intensity', 'RedBeads_Intensity', 'AC_RedBeads', 'AC_BR_Intensity', 'AC_BR_RedBeads', 'Calcein_Intensity', 'GreenBeads_Intensity', 'Calcein_GreenBeads', '# TRAP');
    if no_roi==2
        fprintf(fid_g2,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n', ...
        'Bone type','Bone #','AP_BS', 'AP_R_BS', 'AP_G_BS', 'AP_RG_BS', 'AP_R_RG_BS', 'AP_L_BS', 'AP_NL_BS', 'AP_L_AP', 'AP_NL_AP', 'AP_R_R', 'AP_G_G', 'AP_RG_RG', 'AP_R_RG_R', 'AP_TRAP_BS', 'TRAP_BS', 'TRAP_L_BS', 'TRAP_NL_BS', 'AP_TRAP_R_RG_BS', 'TRAP_on_TRAP', 'TRAP_L_TRAP_on', 'TRAP_NL_TRAP_on', 'GFP_BS', 'GFP_R_BS', 'GFPonly_BS', 'Height_Width', 'Osteocytes_N', 'Cortex_Width', 'AC_Intensity', 'RedBeads_Intensity', 'AC_RedBeads', 'AC_BR_Intensity', 'AC_BR_RedBeads', 'Calcein_Intensity', 'GreenBeads_Intensity', 'Calcein_GreenBeads', '# TRAP');
    end
end

for i=1:size(phrase1,1)
    eval(['number_tr1_',num2str(i),'=0;'])
    %number_tr1_2=0;
    eval(['number_tr2_',num2str(i),'=0;'])
    %number_tr2_2=0;
    eval(['number_tr1_',num2str(i),'_c=0;'])
    %number_tr1_2_c=0;
    eval(['number_tr2_',num2str(i),'_c=0;'])
    %number_tr2_2_c=0;
    eval(['number_tr1_',num2str(i),'_c1(size(bone_number',num2str(i),',1))=0;'])
    %number_tr1_2_c1(length(bone_number2))=0;
    eval(['number_tr2_',num2str(i),'_c1(size(bone_number',num2str(i),',1))=0;'])
    %number_tr2_2_c1(length(bone_number2))=0;
    eval(['no',num2str(i),'=info',num2str(i),';'])
end

NaN_flag=0;             % if Nan_flag==0 --> there is no NaN, if 1 --> there is NaN

new_bone=0;

for loop=1:no_roi
    for b_t=1:size(phrase1,1)       
        eval(['bone_number=bone_number',num2str(b_t),';']);
        for b_n=1:size(bone_number,1)
            if b_n==3
                b_n
            end

            if ana_n==0
                phr1_1=[phrase1{b_t,:},'_',num2str(bone_number(b_n)),'_s1_analysis_tr',num2str(loop),'.mat'];
                phr1_2=[phrase1{b_t,:},'_',num2str(bone_number(b_n)),'_s2_analysis_tr',num2str(loop),'.mat'];
                phr1_3=[phrase1{b_t,:},'_',num2str(bone_number(b_n)),'_s3_analysis_tr',num2str(loop),'.mat'];
            else
%                 if bone_type=='V'
%                     phr1_1=[phrase1{b_t,:},bone_number{b_n},'_h',bone_type,'_s1_analysis_n_tr',num2str(loop),'.mat'];
%                     phr1_2=[phrase1{b_t,:},bone_number{b_n},'_h',bone_type,'_s2_analysis_n_tr',num2str(loop),'.mat'];
%                     phr1_3=[phrase1{b_t,:},bone_number{b_n},'_h',bone_type,'_s3_analysis_n_tr',num2str(loop),'.mat'];
%                 elseif bone_type=='F'
                    phr1_1=[phrase1{b_t,:},'L1', '_s',bone_number{b_n},'_analysis_n_tr',num2str(loop),'.mat'];
                    phr1_2=[phrase1{b_t,:},'L2', '_s',bone_number{b_n},'_analysis_n_tr',num2str(loop),'.mat'];
                    phr1_3=[phrase1{b_t,:},'L3', '_s',bone_number{b_n},'_analysis_n_tr',num2str(loop),'.mat'];
%                 end
%                 phr1_1=[phrase1{b_t,:},'_',num2str(bone_number(b_n)),'_s1_analysis_n_tr',num2str(loop),'.mat'];
%                 phr1_2=[phrase1{b_t,:},'_',num2str(bone_number(b_n)),'_s2_analysis_n_tr',num2str(loop),'.mat'];
%                 phr1_3=[phrase1{b_t,:},'_',num2str(bone_number(b_n)),'_s3_analysis_n_tr',num2str(loop),'.mat'];
            end
            eval(['test1=dir(''',direct,phr1_1,''');']);
            eval(['test2=dir(''',direct,phr1_2,''');']);
            eval(['test3=dir(''',direct,phr1_3,''');']);

            sample=[];
            if isempty(test1)==0
                sample=[1];
            end
            if isempty(test2)==0
                sample=[sample, 2];
            end
            if isempty(test3)==0
                sample=[sample, 3];
            end
            samples_per_bone(b_t, b_n)=length(sample);

            for s=1:3%samples_per_bone(b_t, bone_number(b_n))
                
                nnn=bone_number{b_n};

                if isempty(not_include)==0
                    temp=find(not_include(:,1)==bone_number(b_n));
                    if isempty(temp)==0
                        if isempty(find(not_include(temp,2)==sample(s)))==0
%                             data=[0,0,0,0,0,0,0,0,0];
                            data=zeros(1,size(col_name1,1)); %[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
                            eval(['bone_tr',num2str(loop),'_',num2str(b_t),'(number_tr',num2str(loop),'_',num2str(b_t),'+1,:)=data;'])
                            eval(['bone_tr',num2str(loop),'_',num2str(b_t),'_c(number_tr',num2str(loop),'_',num2str(b_t),'_c+1,:)=data;'])
                            eval(['number_tr',num2str(loop),'_',num2str(b_t),'=number_tr',num2str(loop),'_',num2str(b_t),'+1;'])
                            eval(['number_tr',num2str(loop),'_',num2str(b_t),'_c=number_tr',num2str(loop),'_',num2str(b_t),'_c+1;'])
                            eval(['number_tr',num2str(loop),'_',num2str(b_t),'_c1(',num2str(b_n),')=number_tr',num2str(loop),'_',num2str(b_t),'_c1(',num2str(b_n),')+1;'])
                            eval(['no',num2str(b_t),'.ap(',num2str(b_n),',',num2str(s),')=0;'])
                            eval(['no',num2str(b_t),'.tr(',num2str(nnn),',',num2str(s),')=0;'])
                            continue
                        end
                    end
                end
%                 eval(['aa=test',num2str(sample(s)),';'])
               eval(['aa=test',num2str(s),';'])
                if isempty(aa)==1
                    eval(['bone_tr',num2str(loop),'_',num2str(b_t),'(number_tr',num2str(loop),'_',num2str(b_t),'+1,:)=ones(1,no_column)*NaN;'])
                    eval(['bone_tr',num2str(loop),'_',num2str(b_t),'_c(number_tr',num2str(loop),'_',num2str(b_t),'_c+1,:)=ones(1,no_column)*NaN;'])
                    eval(['number_tr',num2str(loop),'_',num2str(b_t),'=number_tr',num2str(loop),'_',num2str(b_t),'+1;'])
                    eval(['number_tr',num2str(loop),'_',num2str(b_t),'_c=number_tr',num2str(loop),'_',num2str(b_t),'_c+1;'])
                    eval(['number_tr',num2str(loop),'_',num2str(b_t),'_c1(',num2str(b_n),')=number_tr',num2str(loop),'_',num2str(b_t),'_c1(',num2str(b_n),')+1;'])
                else
%                     eval(['phr=phr_',num2str(sample(s)),';']);
%                     eval(['phr1=phr1_',num2str(sample(s)),';']);
%                     eval(['phr=phr_',num2str(s),';']);
                    eval(['phr1=phr1_',num2str(s),';']);
                    eval(['load(''',direct,phr1,''');']);
                    
                    if ana_n==0
                        a=analysis;
                    else
                        a=analysis_n;
                    end
                    LA_BS=a.surface.rA_BS+a.surface.gA_BS+a.surface.rgA_BS;
                    LT_BS=a.surface.rT_BS+a.surface.gT_BS+a.surface.rgT_BS;
                    
                    if include_number_TRAP==0
                        data=[a.surface.A_BS, a.surface.rA_BS, a.surface.gA_BS, a.surface.rgA_BS, a.surface.rA_BS+a.surface.rgA_BS, LA_BS, a.surface.A_only_BS, LA_BS/a.surface.A_BS*100, a.surface.A_only_BS/a.surface.A_BS*100, a.surface.rA_BS/a.surface.r_BS*100, a.surface.gA_BS/a.surface.g_BS*100, a.surface.rgA_BS/a.surface.rg_BS*100, (a.surface.rA_BS+a.surface.rgA_BS)/a.surface.r_BS*100, a.surface.AT_BS, a.surface.T_BS, LT_BS, a.surface.T_only_BS, (a.surface.rAT_BS+a.surface.rgAT_BS),...
                            a.volume.T_on_T, LT_BS/a.surface.T_BS*100, a.surface.T_only_BS/a.surface.T_BS*100, a.surface.GFP_BS, a.surface.rGFP_BS, a.surface.GFP_only_BS, a.roi.ratio, a.cortex.osteocytes.number_per_cortex_area*10^6, a.cortex.thickness, a.intensity.AC, a.intensity.Red_beads, a.intensity.AC_Red_beads, a.intensity.AC_BR, a.intensity.AC_BR_Red_beads, a.intensity.Calcein, a.intensity.Green_beads, a.intensity.Calcein_Green_beads];
                    elseif include_number_TRAP==1
                        data=[a.surface.A_BS, a.surface.rA_BS, a.surface.gA_BS, a.surface.rgA_BS, a.surface.rA_BS+a.surface.rgA_BS, LA_BS, a.surface.A_only_BS, LA_BS/a.surface.A_BS*100, a.surface.A_only_BS/a.surface.A_BS*100, a.surface.rA_BS/a.surface.r_BS*100, a.surface.gA_BS/a.surface.g_BS*100, a.surface.rgA_BS/a.surface.rg_BS*100, (a.surface.rA_BS+a.surface.rgA_BS)/a.surface.r_BS*100, a.surface.AT_BS, a.surface.T_BS, LT_BS, a.surface.T_only_BS, (a.surface.rAT_BS+a.surface.rgAT_BS),...
                            a.volume.T_on_T, LT_BS/a.surface.T_BS*100, a.surface.T_only_BS/a.surface.T_BS*100, a.surface.GFP_BS, a.surface.rGFP_BS, a.surface.GFP_only_BS, a.roi.ratio, a.cortex.osteocytes.number_per_cortex_area*10^6, a.cortex.thickness, a.intensity.AC, a.intensity.Red_beads, a.intensity.AC_Red_beads, a.intensity.AC_BR, a.intensity.AC_BR_Red_beads, a.intensity.Calcein, a.intensity.Green_beads, a.intensity.Calcein_Green_beads, a.volume_T_number];
                    end
%                    data(isnan(data))=0;
                    data(data==0)=NaN;
                    NaN_flag=0;%1-isempty(regexpi(num2str(data),'NaN'))
                    
                    if NaN_flag==1 % if NaN_flag
                        eval(['bone_tr',num2str(loop),'_',num2str(b_t),'(number_tr',num2str(loop),'_',num2str(b_t),'+1,:)=ones(1,length(data))*NaN;'])
                        eval(['number_tr',num2str(loop),'_',num2str(b_t),'=number_tr',num2str(loop),'_',num2str(b_t),'+1;'])
                    else
                        eval(['bone_tr',num2str(loop),'_',num2str(b_t),'(number_tr',num2str(loop),'_',num2str(b_t),'+1,:)=data;'])
                        eval(['bone_tr',num2str(loop),'_',num2str(b_t),'_c(number_tr',num2str(loop),'_',num2str(b_t),'_c+1,:)=data;'])
                        eval(['number_tr',num2str(loop),'_',num2str(b_t),'=number_tr',num2str(loop),'_',num2str(b_t),'+1;'])
                        eval(['number_tr',num2str(loop),'_',num2str(b_t),'_c=number_tr',num2str(loop),'_',num2str(b_t),'_c+1;'])
                        eval(['number_tr',num2str(loop),'_',num2str(b_t),'_c1(',num2str(b_n),')=number_tr',num2str(loop),'_',num2str(b_t),'_c1(',num2str(b_n),')+1;'])
                        if sum(data(1:3))==0
                            eval(['no',num2str(b_t),'.ap(',nnn,',',num2str(s),')=0;'])
                        elseif sum(data(5:9))==0
                            eval(['no',num2str(b_t),'.tr(',nnn,',',num2str(s),')=0;'])
                        end
                            
                    end  % if NaN_flag
                end % if isempty(aa)
            end    % s
        end   % b_n
    end    % b_t
end   % loop


start=1;
for loop=1:no_roi
    total_number_section=zeros(size(phrase1,1),no_column);
    for b_t=1:size(phrase1,1)
%        eval(['data',num2str(loop),'_',num2str(b_t),'=[];'])
        eval(['bone_number=bone_number',num2str(b_t),';']);
        clear number
        switch b_t
            case 1
                no=no1;
            case 2
                no=no2;
            case 3
                no=no3;
            case 4
                no=no4;
        end
        for b_n=1:size(bone_number,1)
            nnn=bone_number{b_n};
            if nnn==3
                nnn;
            end
            for ap_trap=1:3
                if ap_trap==1 % ap
                    end_column=14;
                    
                    for i=1:size(no.ap,1)
                        a=no.ap{i,1};
                        number(i)=str2num(a(double(a)<58 & double(a)>47));
                    end
%                     eval(['ind=no',num2str(b_t),'.ap(',num2str(nnn),',:)'])
                    for i=2:4
                        temp=no.ap{find(number==str2num(nnn)),i};
                        if ~isempty(temp)
                            ind(i-1)=temp;
                        end
                    end
                    eval(['data_for_average=bone_tr',num2str(loop),'_',num2str(b_t),'_c(start:start+number_tr',num2str(loop),'_',num2str(b_t),'_c1(',num2str(b_n),')-1,1:end_column);'])
                    data_for_average=double(data_for_average(find(ind),:));
                    data_for_average1=data_for_average;
                    temp=~isnan(data_for_average);
%                     data_for_average1=reshape(data_for_average(temp),length(find(temp))/3,3);
% %                     eval(['number_section=[sum(no',num2str(b_t),'.ap(',num2str(nnn),',:))*ones(1,3)];'])
%                     number_section=size(data_for_average1,1)*ones(1,3);
%                     phras=['mean_tr',num2str(loop),'_',num2str(b_t),'(',num2str(b_n),',1:3)'];
%                     phras=[phras,'=sum(data_for_average,1)./number_section;'];
%                     eval(phras)
%                     eval(['average=mean_tr',num2str(loop),'_',num2str(b_t),'(',num2str(b_n),',1:3);'])
%                     [x y]=find(data_for_average==0);
%                     for i=1:length(x)
%                         eval(['data_for_average1(x(',num2str(i),'), y(',num2str(i),'))=average(end,y(',num2str(i),'));'])
%                     end
%                     phras=['std_tr',num2str(loop),'_',num2str(b_t),'(',num2str(b_n),',1:3)'];
%                     phras=[phras,'=sqrt(sum((data_for_average1-ones(size(data_for_average1,1),1)*average).^2,1)./(number_section-1));'];
%                     eval(phras)
%                     phras=['RSE_tr',num2str(loop),'_',num2str(b_t),'(',num2str(b_n),',1:3)'];
%                     phras=[phras,'=std_tr',num2str(loop),'_',num2str(b_t),'(',num2str(b_n),',1:3) ./ mean_tr',num2str(loop),'_',num2str(b_t),'(',num2str(b_n),',1:3)./sqrt(number_section)*100;'];
%                     eval(phras)

                    phras=['mean_tr',num2str(loop),'_',num2str(b_t),'(',num2str(b_n),',1:end_column)'];
                    phras=[phras,'=nanmean(data_for_average1,1);'];
                    eval(phras)
%                     phras=['std_tr',num2str(loop),'_',num2str(b_t),'(',num2str(b_n),',1:end_column)'];
%                     phras=[phras,'=nanstd(data_for_average1,0);'];
%                     eval(phras)

                    phras=['std_tr',num2str(loop),'_',num2str(b_t),'(',num2str(b_n),',1:end_column)'];
                    phras=[phras,'=nanstd(data_for_average1,0);'];
                    eval(phras)

                    nan_col=find(sum(temp)==1);
                    
                    phras=['std_tr',num2str(loop),'_',num2str(b_t),'(',num2str(b_n),',nan_col)'];
                    phras=[phras,'=NaN;'];
                    eval(phras)

                

                    phras=['RSE_tr',num2str(loop),'_',num2str(b_t),'(',num2str(b_n),',1:end_column)'];
                    phras=[phras,'=std_tr',num2str(loop),'_',num2str(b_t),'(',num2str(b_n),',1:end_column) ./ mean_tr',num2str(loop),'_',num2str(b_t),'(',num2str(b_n),',1:end_column)./sqrt(sum(temp))*100;'];
                    eval(phras)                    
                elseif ap_trap==2 % trap
                    start_column=15;
                    end_column=21;
                    
                    for i=1:size(no.tr,1)
                        a=no.tr{i,1};
                        number(i)=str2num(a(double(a)<58 & double(a)>47));
                    end
%                     eval(['ind=no',num2str(b_t),'.tr(',num2str(nnn),',:)'])
                    for i=2:4
                        temp=no.tr{find(number==str2num(nnn)),i};
                        if ~isempty(temp)
                            ind(i-1)=temp;
                        end
                    end
                    eval(['data_for_average=bone_tr',num2str(loop),'_',num2str(b_t),'_c(start:start+number_tr',num2str(loop),'_',num2str(b_t),'_c1(',num2str(b_n),')-1,start_column:end_column);'])
                    data_for_average=double(data_for_average(find(ind),:));
                    data_for_average1=data_for_average;
                    temp=~isnan(data_for_average);
%                     data_for_average1=reshape(data_for_average(temp),length(find(temp))/5,5);
% %                     eval(['number_section=[sum(no',num2str(b_t),'.tr(',num2str(nnn),',:))*ones(1,5)];'])
%                     number_section=size(data_for_average1,1)*ones(1,5);
%                     phras=['mean_tr',num2str(loop),'_',num2str(b_t),'(',num2str(b_n),',5:9)'];
%                     phras=[phras,'=sum(data_for_average1,1)./number_section;'];
%                     eval(phras)
%                     eval(['average=mean_tr',num2str(loop),'_',num2str(b_t),'(',num2str(b_n),',5:9);'])
%                     [x y]=find(data_for_average==0);
%                     for i=1:length(x)
%                         eval(['data_for_average1(x(',num2str(i),'), y(',num2str(i),'))=average(end,y(',num2str(i),'));'])
%                     end
%                     phras=['std_tr',num2str(loop),'_',num2str(b_t),'(',num2str(b_n),',5:9)'];
%                     phras=[phras,'=sqrt(sum((data_for_average1-ones(size(data_for_average1,1),1)*average).^2,1)./(number_section-1));'];
%                     eval(phras)
%                     phras=['RSE_tr',num2str(loop),'_',num2str(b_t),'(',num2str(b_n),',5:9)'];
%                     phras=[phras,'=std_tr',num2str(loop),'_',num2str(b_t),'(',num2str(b_n),',5:9) ./ mean_tr',num2str(loop),'_',num2str(b_t),'(',num2str(b_n),',5:9)./sqrt(number_section)*100;'];
%                     eval(phras)


                    phras=['mean_tr',num2str(loop),'_',num2str(b_t),'(',num2str(b_n),',start_column:end_column)'];
                    phras=[phras,'=nanmean(data_for_average1,1);'];
                    eval(phras);
%                     phras=['std_tr',num2str(loop),'_',num2str(b_t),'(',num2str(b_n),',start_column:end_column)'];
%                     phras=[phras,'=nanstd(data_for_average1,0);'];
%                     eval(phras);

                    phras=['std_tr',num2str(loop),'_',num2str(b_t),'(',num2str(b_n),',start_column:end_column)'];
                    phras=[phras,'=nanstd(data_for_average1,0);'];
                    eval(phras);

                    nan_col=find(sum(temp)==1);
                    
                    phras=['std_tr',num2str(loop),'_',num2str(b_t),'(',num2str(b_n),',nan_col+start_column-1)'];
                    phras=[phras,'=NaN;'];
                    eval(phras)
                    
                    phras=['RSE_tr',num2str(loop),'_',num2str(b_t),'(',num2str(b_n),',start_column:end_column)'];
                    phras=[phras,'=std_tr',num2str(loop),'_',num2str(b_t),'(',num2str(b_n),',start_column:end_column) ./ mean_tr',num2str(loop),'_',num2str(b_t),'(',num2str(b_n),',start_column:end_column)./sqrt(sum(temp))*100;'];
                    eval(phras)                    
                elseif ap_trap==3 % GFP
                    start_column=22;
                    end_column=35;
                    for i=1:size(no.tr,1)
                        a=no.tr{i,1};
                        number(i)=str2num(a(double(a)<58 & double(a)>47));
                    end
%                     eval(['ind=no',num2str(b_t),'.tr(',num2str(nnn),',:)'])
                    for i=2:4
                        temp=no.tr{find(number==str2num(nnn)),i};
                        if ~isempty(temp)
                            ind(i-1)=temp;
                        end
                    end
                    eval(['data_for_average=bone_tr',num2str(loop),'_',num2str(b_t),'_c(start:start+number_tr',num2str(loop),'_',num2str(b_t),'_c1(',num2str(b_n),')-1,start_column:end_column);'])
                    data_for_average=double(data_for_average(find(ind),:));
                    data_for_average1=data_for_average;
                    temp=~isnan(data_for_average);
% %                     eval(['number_section=[sum(no',num2str(b_t),'.ap(',num2str(nnn),',:))*ones(1,1)];'])
%                     number_section=size(data_for_average,1);
%                     phras=['mean_tr',num2str(loop),'_',num2str(b_t),'(',num2str(b_n),',4)'];
%                     phras=[phras,'=sum(data_for_average1,1)./number_section;'];
%                     eval(phras)
%                     eval(['average=mean_tr',num2str(loop),'_',num2str(b_t),'(',num2str(b_n),',4);'])
%                     [x y]=find(data_for_average==0);
%                     for i=1:length(x)
%                         eval(['data_for_average1(x(',num2str(i),'), y(',num2str(i),'))=average(end,y(',num2str(i),'));'])
%                     end
%                     phras=['std_tr',num2str(loop),'_',num2str(b_t),'(',num2str(b_n),',4)'];
%                     phras=[phras,'=sqrt(sum((data_for_average1-ones(size(data_for_average1,1),1)*average).^2,1)./(number_section-1));'];
%                     eval(phras)
%                     phras=['RSE_tr',num2str(loop),'_',num2str(b_t),'(',num2str(b_n),',4)'];
%                     phras=[phras,'=std_tr',num2str(loop),'_',num2str(b_t),'(',num2str(b_n),',4) ./ mean_tr',num2str(loop),'_',num2str(b_t),'(',num2str(b_n),',4)./sqrt(number_section)*100;'];
%                     eval(phras)
                    phras=['mean_tr',num2str(loop),'_',num2str(b_t),'(',num2str(b_n),',start_column:end_column)'];
                    phras=[phras,'=nanmean(data_for_average1,1);'];
                    eval(phras);
%                     phras=['std_tr',num2str(loop),'_',num2str(b_t),'(',num2str(b_n),',start_column:end_column)'];
%                     phras=[phras,'=nanstd(data_for_average1,0);'];
%                     eval(phras);

                    phras=['std_tr',num2str(loop),'_',num2str(b_t),'(',num2str(b_n),',start_column:end_column)'];
                    phras=[phras,'=nanstd(data_for_average1,0);'];
                    eval(phras);

                    nan_col=find(sum(temp)==1);
                    
                    phras=['std_tr',num2str(loop),'_',num2str(b_t),'(',num2str(b_n),',nan_col+start_column-1)'];
                    phras=[phras,'=NaN;'];
                    eval(phras)
                    
                    phras=['RSE_tr',num2str(loop),'_',num2str(b_t),'(',num2str(b_n),',start_column:end_column)'];
                    phras=[phras,'=std_tr',num2str(loop),'_',num2str(b_t),'(',num2str(b_n),',start_column:end_column) ./ mean_tr',num2str(loop),'_',num2str(b_t),'(',num2str(b_n),',start_column:end_column)./sqrt(sum(temp))*100;'];
                    eval(phras)  
                end    % if ap_trap
            end     % for ap_trap
            eval(['start=start+number_tr',num2str(loop),'_',num2str(b_t),'_c1(',num2str(b_n),');'])
        end     % for b_n
        start=1;
    end     % for b_t
end     % for loop




for loop=1:no_roi
    for b_t=1:size(phrase1,1)
        b_t
        eval(['bone_number=bone_number',num2str(b_t),';']);
        count=0;
        minus_count(b_t,length(bone_number))=0;
        
        for b_n=1:size(bone_number,1)
            if b_n==3
                b_n;
            end
            nnn=bone_number{b_n};

            if ana_n==0
                phr1_1=[phrase1{b_t,:},num2str(bone_number(b_n)),'_1_analysis_tr',num2str(loop),'.mat'];
                phr1_2=[phrase1{b_t,:},num2str(bone_number(b_n)),'_2_analysis_tr',num2str(loop),'.mat'];
                phr1_3=[phrase1{b_t,:},num2str(bone_number(b_n)),'_3_analysis_tr',num2str(loop),'.mat'];
            else
%                 if bone_type=='V'
%                     phr1_1=[phrase1{b_t,:},bone_number{b_n},'_h',bone_type,'_s1_analysis_n_tr',num2str(loop),'.mat'];
%                     phr1_2=[phrase1{b_t,:},bone_number{b_n},'_h',bone_type,'_s2_analysis_n_tr',num2str(loop),'.mat'];
%                     phr1_3=[phrase1{b_t,:},bone_number{b_n},'_h',bone_type,'_s3_analysis_n_tr',num2str(loop),'.mat'];
%                 elseif bone_type=='F'
                    phr1_1=[phrase1{b_t,:},'L1', '_s',bone_number{b_n},'_analysis_n_tr',num2str(loop),'.mat'];
                    phr1_2=[phrase1{b_t,:},'L2', '_s',bone_number{b_n},'_analysis_n_tr',num2str(loop),'.mat'];
                    phr1_3=[phrase1{b_t,:},'L3', '_s',bone_number{b_n},'_analysis_n_tr',num2str(loop),'.mat'];
%                 end
%                 phr1_1=[phrase1{b_t,:},num2str(bone_number(b_n)),'_1_analysis_n_tr',num2str(loop),'.mat'];
%                 phr1_2=[phrase1{b_t,:},num2str(bone_number(b_n)),'_2_analysis_n_tr',num2str(loop),'.mat'];
%                 phr1_3=[phrase1{b_t,:},num2str(bone_number(b_n)),'_3_analysis_n_tr',num2str(loop),'.mat'];
            end
            eval(['test1=dir(''',direct,phr1_1,''');']);
            eval(['test2=dir(''',direct,phr1_2,''');']);
            eval(['test3=dir(''',direct,phr1_3,''');']);

            sample=[];
            if isempty(test1)==0
                sample=[1];
            end
            if isempty(test2)==0
                sample=[sample, 2];
            end
            if isempty(test3)==0
                sample=[sample, 3];
            end
%             samples_per_bone(b_t, bone_number(b_n))=length(sample);

            for s=1:3%samples_per_bone(b_t, bone_number(b_n))

                if isempty(not_include)==0
                    temp=find(not_include(:,1)==bone_number(b_n));
                    if isempty(temp)==0
                        if isempty(sample)==1
                            minus_count(b_t,b_n)=minus_count(b_t,b_n)+1;
                            continue
                        else
                            if isempty(find(not_include(temp,2)==sample(s)))==0
                                minus_count(b_t,b_n)=minus_count(b_t,b_n)+1;
                                continue
                            end
                        end
                    end
                end
                
%                 if bone_type=='V'
%                     section_name=[phrase1{b_t,:},num2str(nnn),'_h',bone_type,'_s',num2str(s)];
%                 elseif bone_type=='F'
                    section_name=[phrase1{b_t,:},'L',num2str(s),'_s',num2str(nnn)];
%                 end
                
                if loop==1
                    if (s==1 & b_n==1) %| samples_per_bone==1
                        ph=['fprintf(fid_g1,''%s\t%s\'];
                        for col=1:no_column
                            ph=[ph, 't%8.6f\'];
                        end
                        ph=[ph, 'n'', phrase2{b_t,:}, section_name, bone_tr1_',num2str(b_t),'(count+s-sum(minus_count(b_t,:)),:));'];
                        eval(ph);
                    else
                        ph=['fprintf(fid_g1,''\t%s\'];
                        for col=1:no_column
                            ph=[ph, 't%8.6f\'];
                        end
                        ph=[ph, 'n'', section_name, bone_tr1_',num2str(b_t),'(count+s-sum(minus_count(b_t,:)),:));'];
                        eval(ph);
                    end
                elseif loop==2
                    ph=['fprintf(fid_g2,''\t%s\'];
                    for col=1:no_column
                        ph=[ph, 't%8.6f\'];
                    end
                    ph=[ph, 'n'', section_name, bone_tr2_',num2str(b_t),'(count+s-sum(minus_count(b_t,:)),:));'];
                    eval(ph);
                end
            end
            %        end
            
            if isempty(not_include)==0
                temp=find(not_include(:,1)==bone_number(b_n));
                if isempty(temp)==0
                    if isempty(sample)==1
%                         minus_count(b_n)=minus_count(b_n)+1;
                        continue
                    end
                end
            end
            
%             if bone_type=='V'
%                 section_name_mean=[phrase1{b_t,:},nnn,'_h',bone_type,'_mean'];
%                 section_name_std=[phrase1{b_t,:},nnn,'_h',bone_type,'_std'];
%                 section_name_RSE=[phrase1{b_t,:},nnn,'_h',bone_type,'_RSE(%)'];
%             elseif bone_type=='F'
                section_name_mean=[phrase1{b_t,:},'_s',nnn,'_mean'];
                section_name_std=[phrase1{b_t,:},'_s',nnn,'_std'];
                section_name_RSE=[phrase1{b_t,:},'_s',nnn,'_RSE(%)'];
%             end
            
            ph=['fprintf(fid_g',num2str(loop),',''\t%s\'];
            for col=1:no_column
                ph=[ph, 't%8.6f\'];
            end
%                ph=[ph, 'n'', ''mean'', mean_tr1_',num2str(b_t),'(b_n,:));'];
            ph=[ph, 'n'', ''',section_name_mean,''', mean_tr1_',num2str(b_t),'(b_n,:));'];
            eval(ph);
            ph=['fprintf(fid_g',num2str(loop),',''\t%s\'];
            for col=1:no_column
                ph=[ph, 't%8.6f\'];
            end
%                ph=[ph, 'n'', ''std'', std_tr1_',num2str(b_t),'(b_n,:));'];
            ph=[ph, 'n'', ''',section_name_std,''', std_tr1_',num2str(b_t),'(b_n,:));'];
            eval(ph);
            ph=['fprintf(fid_g',num2str(loop),',''\t%s\'];
            for col=1:no_column
                ph=[ph, 't%8.6f\'];
            end
%                ph=[ph, 'n'', ''RSE(%)'', RSE_tr1_',num2str(b_t),'(b_n,:));'];
            ph=[ph, 'n'', ''',section_name_RSE,''', RSE_tr1_',num2str(b_t),'(b_n,:));'];
            eval(ph);

%             count=count+samples_per_bone(b_t, bone_number(b_n))
            count=count+3;%samples_per_bone(b_t, b_n)
        end % b_n
        
%         if bone_type=='V'
%             section_name_mean=[phrase1{b_t,:},bone_type,'_total mean'];
%             section_name_std=[phrase1{b_t,:},bone_type,'_total std'];
%             section_name_RSE=[phrase1{b_t,:},bone_type,'_total RSE(%)'];
%         elseif bone_type=='F'
            section_name_mean=[phrase1{b_t,:},'_total mean'];
            section_name_std=[phrase1{b_t,:},'_total std'];
            section_name_RSE=[phrase1{b_t,:},'_total RSE(%)'];
%         end

        eval(['data_for_std=mean_tr',num2str(loop),'_',num2str(b_t),;])
        for ppp=1:size(data_for_std,2)
            temp=(strfind(num2str((data_for_std(:,ppp))'),'NaN'))
            number_bone_sample(b_t, ppp)=size(data_for_std,1)-length(temp);
        end
        average=nansum(data_for_std,1)./number_bone_sample(b_t,:);
        standard_deviation=sqrt(nansum((data_for_std-ones(size(data_for_std,1),1)*average).^2)./(number_bone_sample(b_t,:)-1));
        RSE=standard_deviation./average/sqrt(number_bone_sample(b_t))*100;

        eval(['fid_g=fid_g',num2str(loop),';'])
        ph=['fprintf(fid_g,''\t%s\'];
        for col=1:no_column
            ph=[ph, 't%8.6f\'];
        end
        ph=[ph, 'n'', ''',section_name_mean,''', average);'];
        eval(ph);
        ph=['fprintf(fid_g,''\t%s\'];
        for col=1:no_column
            ph=[ph, 't%8.6f\'];
        end
        ph=[ph, 'n'', ''',section_name_std,''', standard_deviation);'];
        eval(ph);
        ph=['fprintf(fid_g,''\t%s\'];
        for col=1:no_column
            ph=[ph, 't%8.6f\'];
        end
        ph=[ph, 'n'', ''',section_name_RSE,''', RSE);'];
        eval(ph);
    end % b_t
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       t-test

if size(phrase1,1)>1
    for j=1:size(phrase1,1)-1
        for k=j+1:size(phrase1,1)
            for loop=1:no_roi
                if section_effect==1
                    eval(['b=(mean_tr',num2str(loop),'_',num2str(j),');'])
                    eval(['c=(mean_tr',num2str(loop),'_',num2str(k),');'])
                else
                    eval(['b=(bone_tr',num2str(loop),'_',num2str(j),'_c);'])
                    eval(['c=(bone_tr',num2str(loop),'_',num2str(k),'_c);'])
                end


                for i=1:size(b,2)

                    bb=b(:,i);
                    cc=c(:,i);

                    b_non_zero=bb(find(bb~=-Inf));
                    c_non_zero=cc(find(cc~=-Inf));
                    eval(['[h_b_c_tr',num2str(loop),'(',num2str(i),'), p_b_c_tr',num2str(loop),'(',num2str(i),'), ci_b_c_tr',num2str(loop),'(:,',num2str(i),'), stats_b_c_tr',num2str(loop),'(',num2str(i),')]=ttest2(b_non_zero,c_non_zero,[],[],''unequal'');'])
                end
            end

%             p1=['t-test',num2str(j),' vs ',num2str(k)];
            for loop=1:no_roi
                eval(['a=phrase2{',num2str(j),',:};'])
                eval(['b=phrase2{',num2str(k),',:};'])
                p1=['t-test : ',a,' vs ',b];
                eval(['pp=p_b_c_tr',num2str(loop),';']);
                eval(['hh=h_b_c_tr',num2str(loop),';']);
                eval(['fid_g=fid_g',num2str(loop),';']);
%                 if loop==1
                    ph=['fprintf(fid_g,''%s\t%s\'];
                    for col=1:no_column
                        ph=[ph, 't%4.2e\'];
                    end
 %                   ph=[ph, 'n'', ''t-test'', ''p_value'', p_b_c_tr1);'];
                    ph=[ph, 'n'', ''',p1,''', ''p_value'', pp);'];
                    eval(ph);
                    ph=['fprintf(fid_g,''%s\t\'];
                    for col=1:no_column
                        ph=[ph, 't%d\'];
                    end
                    ph=[ph, 'n'', ''',p1,''', hh);'];
%                     eval(ph);


%                 elseif loop==2
%                     ph=['fprintf(fid_g,''%s\t%s\'];
%                     for col=1:no_column
%                         ph=[ph, 't%4.2e\'];
%                     end
%                     ph=[ph, 'n'', ',p1,', ''p_value'', p_b_c_tr2);'];
%                     eval(ph);
%                     ph=['fprintf(fid_g,''%s\t\'];
%                     for col=1:no_column
%                         ph=[ph, 't%d\'];
%                     end
%                     ph=[ph, 'n'', ',p1,', h_b_c_tr2);'];
%                     eval(ph);
%                 end % if loop
            end % loop
        end % k
    end % j
end % if

fprintf(fid_g1,'\n');
for loop=1:no_roi
    for b_t=1:size(phrase1,1)   
%         if bone_type=='V'
%             section_name_mean=[phrase1{b_t,:},bone_type,'_total mean'];
%         elseif bone_type=='F'
            section_name_mean=[phrase1{b_t,:},'_total mean'];
%         end
        eval(['data_for_std=mean_tr',num2str(loop),'_',num2str(b_t),;])
        for ppp=1:size(data_for_std,2)
            temp=(strfind(num2str((data_for_std(:,ppp))'),'NaN'));
            number_bone_sample(b_t, ppp)=size(data_for_std,1)-length(temp);
        end
        average=nansum(data_for_std,1)./number_bone_sample(b_t,:);

        eval(['fid_g=fid_g',num2str(loop),';'])
        ph=['fprintf(fid_g,''\t%s\'];
        for col=1:no_column
            ph=[ph, 't%8.6f\'];
        end
        ph=[ph, 'n'', ''',section_name_mean,''', average);'];
        eval(ph);
    end
end


fclose all

for i=1:size(phrase1,1)
    eval(['ma(',num2str(i),')=size(RSE_tr1_',num2str(i),',1);'])
end
max_RSE=max(ma);

for i=1:size(phrase1,1)
    eval(['si=size(RSE_tr1_',num2str(i),',1);'])
    if si<max_RSE
        eval(['RSE_tr1_',num2str(i),'(size(RSE_tr1_',num2str(i),',1),end)=0;'])
    end
end


for i=1:size(phrase1,1)
    eval(['a=isnan(RSE_tr1_',num2str(i),');'])
    eval(['RSE_tr1_',num2str(i),'(a)=0;'])
end


color={'red';
	   'blue';
       'black';
       'magenta';
       'yellow';
       'cyan';
       'green'};

ma=0;
for i=1:size(phrase1,1)
    eval(['ma=max(ma,max(RSE_tr1_',num2str(i),'(:)));'])
end
ma=ma*1.2;

for loop=1:size(RSE_tr1_1,2)
    plot(RSE_tr1_1(:,loop),color{1});
    for i=2:size(phrase1,1)
        hold on;
        eval(['plot(RSE_tr1_',num2str(i),'(:,',num2str(loop),'), color{',num2str(i),'});'])
    end
    aaa=[];
    for i=1:size(bone_number1,1)
        aaa=[aaa,' ',bone_number1{i,1}];
    end
    ph=['legend([phrase2{1,:}, '' (',aaa,')'']'];
    for j=2:size(phrase1,1)
        eval(['a=bone_number',num2str(j),';'])
        aaa=[];
        for i=1:size(a,1)
            aaa=[aaa,' ',a{i,1}];
        end
        ph=[ph,',[phrase2{',num2str(j),',:}, '' (',aaa,')'']'];
    end
    ph=[ph,');'];
    eval(ph)
    
    xlabel('Bone number');
    ylabel('RSE (%)');
    title(['RSE of each bone for ',cell2mat(col_name1(loop))])
    axis([1 size(RSE_tr1_1,1) 0 ma])
    output_image_file=[direct,'RSE_cell_',cell2mat(col_name2(loop))];
    eval(['print -djpeg100 ''',output_image_file,''''])
    close all
end

% fclose all
% 
% if size(RSE_tr1_1,1)>size(RSE_tr1_2,1)
%     RSE_tr1_2(size(RSE_tr1_1,1),end)=0;
% elseif size(RSE_tr1_1,1)<size(RSE_tr1_2,1)
%     RSE_tr1_1(size(RSE_tr1_2,1),end)=0;
% end
% a=isnan(RSE_tr1_1);
% RSE_tr1_1(a)=0;
% a=isnan(RSE_tr1_2);
% RSE_tr1_2(a)=0;
% % if phrase2(1,1)=='F'
% %     color1='r';
% % elseif phrase2(1,1)=='M'
% %     color1='b';
% % else
% %     color1='b';
% % end
% % if phrase2(2,1)=='F'
% %     color2='r';
% % elseif phrase2(2,1)=='M'
% %     color2='b';
% % else
% %     color2='r';
% % end
% 
% color1='r';
% color2='b';
% 
% ma=max(max(RSE_tr1_1(:)),max(RSE_tr1_2(:)))*1.1;
% for loop=1:size(RSE_tr1_1,2)
%     plot(RSE_tr1_1(:,loop),color1);
%     hold on; plot(RSE_tr1_2(:,loop), color2);
%     legend([phrase2(1,:),' (',num2str(bone_number1),')'],[phrase2(2,:),' (',num2str(bone_number2),')']);
%     xlabel('Bone number');
%     ylabel('RSE (%)');
%     title(['RSE of each bone for ',cell2mat(col_name1(loop))])
%     axis([1 size(RSE_tr1_1,1) 0 ma])
%     output_image_file=[direct,'RSE_cell_',cell2mat(col_name2(loop))];
%     eval(['print -djpeg100 ''',output_image_file,''''])
%     close all
% end
