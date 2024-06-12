function control_data_write_2    %(o_direct, exp_name, bone_type)


a=textread('dir_info.txt','%s', 2);
home_dir=a{1,:};
o_direct=home_dir;
phr_temp=a{2,:};
bone_type=home_dir(end-1);
exp_name=a{2,:};
if ~isempty(strfind(home_dir,'\'))
    delimeter='\';
elseif ~isempty(strfind(home_dir,'/'))
    delimeter='/';
end
% delimeter='\';
% o_direct='E:\seh00004\Het_KOMP\CCC_E11_hF\';
% exp_name='CCC_E11';
% bone_type='F';                                % Femur, Vertebra
bone_type=home_dir(end-1);
phrase2={'Female';                              % Exp name for output excel file
      'Male'};
temp=strfind(exp_name,'DO-');
if temp
    DO_flag=1;
else
    DO_flag=0;
end
common=[exp_name,'_h',bone_type,'_'];

% exclude_list=[1 3 1 0 1 0;              % [bt, bn, section, All, AP, TRAP]
%     1 3 2 1 0 0;
%     2 4 3 0 0 1];

% phrase1={'DO-LRD01_G09_hF_F';'DO-LRD01_G09_hF_M'};
% phrase1={[common,phrase2{1}(1)];[common,phrase2{2}(1)]};
phrase1='{';
for i=1:length(phrase2)
    n=['''',common,phrase2{i}(1),''';'];
    eval(['phrase1=[phrase1,n];'])
end
phrase1=[phrase1,'}'];
phrase1=eval(phrase1);

direct=[o_direct,'01_Submitted',delimeter,'Layers',delimeter,'Images',delimeter,'data',delimeter'];
bone_n=[8 8];
if exist('exclude_list')
    write_data5(direct, phrase1, phrase2, bone_type, delimeter, bone_n, exclude_list)
    write_data7_cell(direct, phrase1, phrase2, bone_type, delimeter, bone_n, exclude_list)
else
    write_data5(direct, phrase1, phrase2, bone_type, delimeter, bone_n)
    write_data7_cell(direct, phrase1, phrase2, bone_type, delimeter, bone_n)
end

% bone_excel_file='Z:\DO\DO-LRD01_G09\01_Submitted\Layers\Images\data\DO-LRD01_G09_hF_analysis1_tr1.xls';
% cell_excel_file='Z:\DO\DO-LRD01_G09\01_Submitted\Layers\Images\data\DO-LRD01_G09_hF_analysis1_cell_tr1.xls';
% output_excel_file='Z:\DO\DO-LRD01_G09\04_CalculatedData\DO-LRD02_G09.xls';
bone_excel_file=[direct,common,'analysis1_tr1.xlsx'];
cell_excel_file=[direct,common,'analysis1_cell_tr1.xlsx'];
if DO_flag==0
    output_excel_file=[o_direct,'03_CalculatedData',delimeter,exp_name,bone_type,'.xlsx'];
else
    output_excel_file=[o_direct,'03_CalculatedData',delimeter,exp_name,'.xlsx'];
end

template='\\cse-3m58853\seh00004\Excel_template_withGraphs.xlsx';
if isempty(dir([o_direct,'03_CalculatedData',delimeter]))
    mkdir([o_direct,'03_CalculatedData',delimeter])
end
eval(['copyfile(''',[template],''',''',[output_excel_file],''',''f'');'])
read_excel_data_4(bone_excel_file, cell_excel_file, output_excel_file, bone_n)

function write_data5(varargin)

% function write_data5(direct, phrase1, phrase2, bone_type, delimeter)

% direct='Z:\DO\DO-LRD01_G07\01_Submitted\Layers\Images\data\';
% phrase1={'DO-LRD01_G07_hF_F';'DO-LRD01_G07_hF_M'};
% phrase2={'Female';                              % Exp name for output excel file
%       'Male'};
% bone_type='F';                                % Femur, Vertebra
% delimeter='\';

% bone_number1=[792,793,794,858,859,860,861,862];                          % Control_Female
% bone_number2=[47,48,49,53,54,57,58,59];                          % Control_Male

direct=varargin{1};
phrase1=varargin{2};
phrase2=varargin{3};
bone_type=varargin{4};
delimeter=varargin{5};
number_bones=varargin{6};
if nargin<7
    exclude_list=[0 0 0 0 0 0 0];
else
    exclude_list=varargin{7};
end
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
% not_include=[322,3;328,2];

section_effect=1;           % section_effect 1 : use mean value of sections of 1 mouce as 1 data
                            % section_effect 0 : use all of the values of sections of all mice
col_name1={'R_BS';      % 1
    'G_BS';             % 2
    'R_G';              % 3              %
    'Ronly_BS';        % 4
    'Gonly_BS';        % 5
    'sLS_BS';           % 6
    'dLS_BS';           % 7
    'LS_BS';            % 8
    'R_LS';             % 9             %
    'G_LS';             % 10             %
    'MS_BS';            % 11
    'sLS_LS';           % 12
    'dLS_LS';           % 13
    'dLS_sLS';          % 14
    'MAR_um';          % 15
    'BFR';              % 16
    'BV_TV';            % 17
    'Tb.Th';            % 18
    'Tb.N';             % 19
    'Tb.Sp'};           % 20
col_name2={'R_BS';      % 1
    'G_BS';             % 2
    'R_G';              % 3              %
    'R_nly_BS';        % 4
    'Gonly_BS';        % 5
    'sLS_BS';           % 6
    'dLS_BS';           % 7
    'LS_BS';            % 8
    'R_LS';             % 9             %
    'G_LS';             % 10             %
    'MS_BS';            % 11
    'sLS_LS';           % 12
    'dLS_LS';           % 13
    'dLS_sLS';          % 14
    'MAR_um';          % 15
    'BFR';              % 16
    'BV_TV';            % 17
    'Tb.Th';            % 18
    'Tb.N';             % 19
    'Tb.Sp'};           % 20
fclose all

% if bone_type=='V'
%     ana_n=0;
% else
    ana_n=0;               % ana_n=0 : first analysis, ana_n=1 : second analysis
% end


no_roi=1;

tic

if strcmp(direct(end),delimeter)~=1
    direct=[direct,delimeter];
end

fclose('all')

aaa=phrase1{1,:};

eval(['out_file1=''',direct,aaa(1:end-1),'analysis1_tr1.xlsx'';'])
% eval(['fid_g1 = fopen(''',out_file1,''', ''w'');'])
% fid_g1 = fopen(out_file1,'w');
% fprintf(fid_g1,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n', ...
%     'Bone type','Bone #','R_BS', 'G_BS', 'R_G', 'Ronly_BS', 'Gonly_BS', 'sLS_BS', 'dLS_BS', 'LS_BS', 'R_LS', 'G_LS', 'MS_BS', 'sLS_LS', 'dLS_LS', 'dLS_sLS', 'MAR_um', 'BFR', 'BV_TV', 'Tb.Th', 'Tb.N', 'Tb.Sp');    tem='';

tem(1)={'Bone type'};
tem(2)={'Bone #'};

for k=1:length(col_name1)
    tem(k+2)=col_name1(k,:);
end
% tem={tem};

% if no_roi==2
%     eval(['out_file2=''',direct,phrase1{1,:},'analysis1_tr2.txt'';'])
%     fid_g2 = fopen(out_file2, 'w');
%     fprintf(fid_g2,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n', ...
%         'Bone type','Bone #','R_BS', 'G_BS', 'R_G', 'Ronly_BS', 'Gonly_BS', 'sLS_BS', 'dLS_BS', 'LS_BS', 'R_LS', 'G_LS', 'MS_BS', 'sLS_LS', 'dLS_LS', 'dLS_sLS', 'MAR_um', 'BFR', 'BV_TV', 'Tb.Th', 'Tb.N', 'Tb.Sp');
% end


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
end

NaN_flag=0;             % if Nan_flag==0 --> there is no NaN, if 1 --> there is NaN

new_bone=0;


for loop=1:no_roi
    for b_t=1:size(phrase1,1)
        b_t;
        eval(['bone_number=bone_number',num2str(b_t),';']);
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Remove 9th bone, if empty
%
%   Dec. 23 2020
%
        if isempty(bone_number{9,2}) | number_bones(1)==8
            bone_number=bone_number(1:8,:);
        end
%
%   Dec. 23 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
        for b_n=1:size(bone_number,1)

            bn=str2num(bone_number{b_n,1});


            if ana_n==0
%                 phr1_1=[phrase1{b_t,:},'_',num2str(bone_number(b_n)),'_s1_analysis_tr',num2str(loop),'.mat'];
%                 phr1_2=[phrase1{b_t,:},'_',num2str(bone_number(b_n)),'_s2_analysis_tr',num2str(loop),'.mat'];
%                 phr1_3=[phrase1{b_t,:},'_',num2str(bone_number(b_n)),'_s3_analysis_tr',num2str(loop),'.mat'];
                phr1_1=[phrase1{b_t,:},'L1', '_s',bone_number{b_n},'_analysis_tr',num2str(loop),'.mat'];
                phr1_2=[phrase1{b_t,:},'L2', '_s',bone_number{b_n},'_analysis_tr',num2str(loop),'.mat'];
                phr1_3=[phrase1{b_t,:},'L3', '_s',bone_number{b_n},'_analysis_tr',num2str(loop),'.mat'];
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
%                 phr1_1=[phrase1{b_t,:},'0',num2str(bone_number(b_n)),'_1_analysis_n_tr',num2str(loop),'.mat'];
%                 phr1_2=[phrase1{b_t,:},'0',num2str(bone_number(b_n)),'_2_analysis_n_tr',num2str(loop),'.mat'];
%                 phr1_3=[phrase1{b_t,:},'0',num2str(bone_number(b_n)),'_3_analysis_n_tr',num2str(loop),'.mat'];
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
            
            for s=1:3%samples_per_bone(b_t,bone_number(b_n))

                if isempty(not_include)==0
                    temp=find(not_include(:,1)==bn);
                    if isempty(temp)==0
                        if isempty(find(not_include(temp,2)==sample(s)))==0
                            continue
                        end
                    end
                end

%                 eval(['aa=test',num2str(sample(s)),';'])
                eval(['aa=test',num2str(s),';'])
%                 if isempty(aa)==1 | ~isempty(exclude_list(find(exclude_list(find(exclude_list(find(exclude_list(:,1)==b_t),2)==bn),3)==s)))
                temp_list=exclude_list(exclude_list(:,1)==b_t,:);
                temp_list=temp_list(temp_list(:,2)==bn,:);
                temp_list=temp_list(temp_list(:,3)==s,:);

                if isempty(aa)==1 | ~isempty(temp_list(temp_list(:,4)==1,:))
%                     eval(['bone_tr',num2str(loop),'_',num2str(b_t),'(number_tr',num2str(loop),'_',num2str(b_t),'+1,:)=ones(1,17)*NaN;'])
%                     eval(['number_tr',num2str(loop),'_',num2str(b_t),'=number_tr',num2str(loop),'_',num2str(b_t),'+1;'])
                    eval(['bone_tr',num2str(loop),'_',num2str(b_t),'(number_tr',num2str(loop),'_',num2str(b_t),'+1,:)=ones(1,20)*NaN;'])
                    eval(['bone_tr',num2str(loop),'_',num2str(b_t),'_c(number_tr',num2str(loop),'_',num2str(b_t),'_c+1,:)=ones(1,20)*NaN;'])
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
                    a.surface.LS_BS=a.surface.LS/a.surface.BS*100;
                    if a.surface.LS==0
                        a.surface.sLS_LS=NaN;
                        a.surface.rg_LS=NaN;
                    end
                    data=[a.surface.r_BS, a.surface.g_BS, a.surface.r_BS/a.surface.g_BS*100, a.surface.r_only_BS, a.surface.g_only_BS, a.surface.sLS_BS, a.surface.rg_BS, a.surface.LS_BS, a.surface.r_BS/a.surface.LS_BS*100, a.surface.g_BS/a.surface.LS_BS*100,...
                        a.surface.sLS_BS/2+a.surface.rg_BS, a.surface.sLS_LS, a.surface.rg_LS, a.surface.dLS_sLS, a.thickness.Mar, a.thickness.BFR, a.volume.BV_TV, a.thickness.Tb_Th, a.thickness.Tb_N*1000, a.thickness.Tb_Sp];
                    data(data==Inf)=NaN;
                    if data(3)>1000
                        data(3)=NaN;
                    end
%                     NaN_flag=1-isempty(regexpi(num2str(data),'NaN'));
%                     
%                     if NaN_flag==1
%                         eval(['bone_tr',num2str(loop),'_',num2str(b_t),'(number_tr',num2str(loop),'_',num2str(b_t),'+1,:)=ones(1,17)*NaN;'])
%                         eval(['number_tr',num2str(loop),'_',num2str(b_t),'=number_tr',num2str(loop),'_',num2str(b_t),'+1;'])
% 
%                     else % if NaN_flag
                        eval(['bone_tr',num2str(loop),'_',num2str(b_t),'(number_tr',num2str(loop),'_',num2str(b_t),'+1,:)=data;'])
                        eval(['bone_tr',num2str(loop),'_',num2str(b_t),'_c(number_tr',num2str(loop),'_',num2str(b_t),'_c+1,:)=data;'])
                        eval(['number_tr',num2str(loop),'_',num2str(b_t),'=number_tr',num2str(loop),'_',num2str(b_t),'+1;'])
                        eval(['number_tr',num2str(loop),'_',num2str(b_t),'_c=number_tr',num2str(loop),'_',num2str(b_t),'_c+1;'])
                        eval(['number_tr',num2str(loop),'_',num2str(b_t),'_c1(',num2str(b_n),')=number_tr',num2str(loop),'_',num2str(b_t),'_c1(',num2str(b_n),')+1;'])
%                     end  % if NaN_flag
                end % if isempty(aa)
            end    % s
        end   % b_n
    end    % b_t
end   % loop




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       t-test


start=1;
for loop=1:no_roi
    for b_t=1:size(phrase1,1)       
        eval(['bone_number=bone_number',num2str(b_t),';']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Remove 9th bone, if empty
%
%   Dec. 23 2020
%
        if isempty(bone_number{9,2}) | number_bones(1)==8
            bone_number=bone_number(1:8,:);
        end
%
%   Dec. 23 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        for b_n=1:size(bone_number,1)
            eval(['tt=number_tr',num2str(loop),'_',num2str(b_t),'_c1(',num2str(b_n),');'])
            if tt==0
                eval(['mean_tr',num2str(loop),'_',num2str(b_t),'(',num2str(b_n),',:)=ones(1,17)*NaN;'])
                eval(['std_tr',num2str(loop),'_',num2str(b_t),'(',num2str(b_n),',:)=ones(1,17)*NaN;'])
                eval(['RSE_tr',num2str(loop),'_',num2str(b_t),'(',num2str(b_n),',:)=ones(1,17)*NaN;'])
            else
                eval(['temp=bone_tr',num2str(loop),'_',num2str(b_t),'_c(start:start+number_tr',num2str(loop),'_',num2str(b_t),'_c1(',num2str(b_n),')-1,:);'])
                number=0;
                clear data
%                 for i=1:3
                for i=1:tt
                    %if ~isnan(temp(i,1))       % October 16, 2014
                        number=number+1;
                        data(number,:)=temp(i,:);
                    %end                        % October 16, 2014
                end
                phras=['mean_tr',num2str(loop),'_',num2str(b_t),'(',num2str(b_n),',:)'];
                phras=[phras,'=mean(data,1,''omitnan'');'];
                eval(phras)
                
                if size(data,1)>1
                    phras=['std_tr',num2str(loop),'_',num2str(b_t),'(',num2str(b_n),',:)'];
                    phras=[phras,'=std(data,''omitnan'');'];
                    eval(phras)
                else
                    phras=['std_tr',num2str(loop),'_',num2str(b_t),'(',num2str(b_n),',:)'];
                    phras=[phras,'=NaN*[1:',num2str(size(temp,2)),'];'];
                    eval(phras)
                end
                
                phras=['RSE_tr',num2str(loop),'_',num2str(b_t),'(',num2str(b_n),',:)'];
                phras=[phras,'=std_tr',num2str(loop),'_',num2str(b_t),'(',num2str(b_n),',:) ./ mean_tr',num2str(loop),'_',num2str(b_t),'(',num2str(b_n),',:)/sqrt(number)*100;'];
                eval(phras)
                eval(['start=start+number_tr',num2str(loop),'_',num2str(b_t),'_c1(',num2str(b_n),');'])
            end
        end
        start=1;

    end
end



row=2;
for loop=1:no_roi
    for b_t=1:size(phrase1,1)       
        eval(['bone_number=bone_number',num2str(b_t),';']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Remove 9th bone, if empty
%
%   Dec. 23 2020
%
        if isempty(bone_number{9,2}) | number_bones(1)==8
            bone_number=bone_number(1:8,:);
        end
%
%   Dec. 23 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
count=0;
        minus_count(b_t,length(bone_number))=0;
        
        for b_n=1:size(bone_number,1)
            if b_n==3
                b_n;
            end
            
%             nnn=num2str(bone_number(b_n));
            nnn=bone_number{b_n,1};

            if ana_n==0
%                 phr1_1=[phrase1{b_t,:},num2str(bone_number(b_n)),'_1_analysis_tr',num2str(loop),'.mat'];
%                 phr1_2=[phrase1{b_t,:},num2str(bone_number(b_n)),'_2_analysis_tr',num2str(loop),'.mat'];
%                 phr1_3=[phrase1{b_t,:},num2str(bone_number(b_n)),'_3_analysis_tr',num2str(loop),'.mat'];
                phr1_1=[phrase1{b_t,:},'L1', '_s',bone_number{b_n},'_analysis_tr',num2str(loop),'.mat'];
                phr1_2=[phrase1{b_t,:},'L2', '_s',bone_number{b_n},'_analysis_tr',num2str(loop),'.mat'];
                phr1_3=[phrase1{b_t,:},'L3', '_s',bone_number{b_n},'_analysis_tr',num2str(loop),'.mat'];
            else
%                 phr1_1=[phrase1{b_t,:},num2str(bone_number(b_n)),'_1_analysis_n_tr',num2str(loop),'.mat'];
%                 phr1_2=[phrase1{b_t,:},num2str(bone_number(b_n)),'_2_analysis_n_tr',num2str(loop),'.mat'];
%                 phr1_3=[phrase1{b_t,:},num2str(bone_number(b_n)),'_3_analysis_n_tr',num2str(loop),'.mat'];
%                 phr1_1=[phrase1{b_t,:},'0',num2str(bone_number(b_n)),'_1_analysis_n_tr',num2str(loop),'.mat'];
%                 phr1_2=[phrase1{b_t,:},'0',num2str(bone_number(b_n)),'_2_analysis_n_tr',num2str(loop),'.mat'];
%                 phr1_3=[phrase1{b_t,:},'0',num2str(bone_number(b_n)),'_3_analysis_n_tr',num2str(loop),'.mat'];

%                 if bone_type=='V'
%                     phr1_1=[phrase1{b_t,:},bone_number{b_n},'_h',bone_type,'_s1_analysis_n_tr',num2str(loop),'.mat'];
%                     phr1_2=[phrase1{b_t,:},bone_number{b_n},'_h',bone_type,'_s2_analysis_n_tr',num2str(loop),'.mat'];
%                     phr1_3=[phrase1{b_t,:},bone_number{b_n},'_h',bone_type,'_s3_analysis_n_tr',num2str(loop),'.mat'];
%                 elseif bone_type=='F'
                    phr1_1=[phrase1{b_t,:},'L1', '_s',bone_number{b_n},'_analysis_n_tr',num2str(loop),'.mat'];
                    phr1_2=[phrase1{b_t,:},'L2', '_s',bone_number{b_n},'_analysis_n_tr',num2str(loop),'.mat'];
                    phr1_3=[phrase1{b_t,:},'L3', '_s',bone_number{b_n},'_analysis_n_tr',num2str(loop),'.mat'];
%                 end
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
%             samples_per_bone=length(sample);

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
                eval(['d=bone_tr',num2str(loop),'_',num2str(b_t),'(',num2str(count),'+',num2str(s),'-sum(minus_count(',num2str(b_t),',:)),:);'])
%                 if (s==1 & b_n==1) %| samples_per_bone(b_t, b_n)==1
%                     fprintf(fid_g1,'%s\t%s\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\n',...
%                         phrase2{b_t,:}, section_name, d);
%                 else
%                     fprintf(fid_g1,'\t%s\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\n',...
%                         section_name, d);
%                 end
                if (s==1 & b_n==1) %| samples_per_bone(b_t, b_n)==1
                    tem(row,1)={phrase2{b_t,:}};
                end
                tem(row,2)={section_name};
                for var_loop=1:length(d)
                    tem(row,var_loop+2)={d(var_loop)};
                end
                row=row+1;
                    
            end
            
            if isempty(not_include)==0
                temp=find(not_include(:,1)==bone_number(b_n));
                if isempty(temp)==0
                    if isempty(sample)==1
%                         minus_count(b_n)=minus_count(b_n)+1;
                        continue
                    end
                end
            end
            
%             section_name_mean=[phrase1{b_t,:},num2str(nnn),'_mean'];
%             section_name_std=[phrase1{b_t,:},num2str(nnn),'_std'];
%             section_name_RSE=[phrase1{b_t,:},num2str(nnn),'_RSE(%)'];

%             if bone_type=='V'
%                 section_name_mean=[phrase1{b_t,:},nnn,'_h',bone_type,'_mean'];
%                 section_name_std=[phrase1{b_t,:},nnn,'_h',bone_type,'_std'];
%                 section_name_RSE=[phrase1{b_t,:},nnn,'_h',bone_type,'_RSE(%)'];
%             elseif bone_type=='F'
                section_name_mean=[phrase1{b_t,:},'_s',nnn,'_mean'];
                section_name_std=[phrase1{b_t,:},'_s',nnn,'_std'];
                section_name_RSE=[phrase1{b_t,:},'_s',nnn,'_RSE(%)'];
%             end
            eval(['d_mean=mean_tr',num2str(loop),'_',num2str(b_t),'(',num2str(b_n),',:);'])
            eval(['d_std=std_tr',num2str(loop),'_',num2str(b_t),'(',num2str(b_n),',:);'])
            eval(['d_RSE=RSE_tr',num2str(loop),'_',num2str(b_t),'(',num2str(b_n),',:);'])
%             eval(['fid_g=fid_g',num2str(loop),';'])

%             fprintf(fid_g,'\t%s\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\n',...
%                 section_name_mean, d_mean);
%             fprintf(fid_g,'\t%s\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\n',...
%                 section_name_std, d_std);
%             fprintf(fid_g,'\t%s\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\n',...
%                 section_name_RSE, d_RSE);

            tem(row,2)={section_name_mean};
            for var_loop=1:length(d_mean)
                tem(row,var_loop+2)={d_mean(var_loop)};
            end
            row=row+1;
            tem(row,2)={section_name_std};
            for var_loop=1:length(d_std)
                tem(row,var_loop+2)={d_std(var_loop)};
            end
            row=row+1;                
            tem(row,2)={section_name_RSE};
            for var_loop=1:length(d_RSE)
                tem(row,var_loop+2)={d_RSE(var_loop)};
            end
            row=row+1;
            count=count+3;%samples_per_bone(b_t, b_n);
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
%         eval(['fid_g=fid_g',num2str(loop),';'])
        eval(['len=(size(mean_tr',num2str(loop),'_',num2str(b_t),',1)-minus_count(',num2str(b_n),'));'])
%         eval(['me=nansum(mean_tr',num2str(loop),'_',num2str(b_t),',1)/len;'])
%         eval(['stdev=sqrt(nansum((mean_tr',num2str(loop),'_',num2str(b_t),'-ones(size(mean_tr',num2str(loop),'_',num2str(b_t),',1),1)*me).^2)/(len-1));'])
        eval(['me=mean(mean_tr',num2str(loop),'_',num2str(b_t),',1,''omitnan'');'])
        eval(['stdev=std(mean_tr',num2str(loop),'_',num2str(b_t),',0,''omitnan'');'])
%         fprintf(fid_g1,'\t%s\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\n',...
%             section_name_mean, me);
%         fprintf(fid_g1,'\t%s\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\n',...
%             section_name_std, stdev);
%         fprintf(fid_g1,'\t%s\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\n',...
%             section_name_RSE, stdev./me/sqrt(len)*100);

        tem(row,2)={section_name_mean};
        for var_loop=1:length(me)
            tem(row,var_loop+2)={me(var_loop)};
        end
        row=row+1;
        tem(row,2)={section_name_std};
        for var_loop=1:length(stdev)
            tem(row,var_loop+2)={stdev(var_loop)};
        end
        row=row+1;                
        tem(row,2)={section_name_RSE};
        for var_loop=1:length(d_RSE)
            tem(row,var_loop+2)={(stdev(var_loop)./me(var_loop)/sqrt(len)*100)};
        end
        row=row+1;
    end % b_t
end

   
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
                    eval(['[h_b_c_tr',num2str(loop),'(',num2str(i),'), p_b_c_tr',num2str(loop),'(',num2str(i),'), ci_b_c_tr',num2str(loop),'(:,',num2str(i),'), stats_b_c_tr',num2str(loop),'(',num2str(i),')]=ttest2(c(:,',num2str(i),'),b(:,',num2str(i),'),[],[],''unequal'');'])
                end
            end

            eval(['a=phrase2{',num2str(j),',:};'])
            eval(['b=phrase2{',num2str(k),',:};'])
            p1=['t-test : ',a,' vs ',b];
            eval(['pp=p_b_c_tr',num2str(loop),';']);
            eval(['hh=h_b_c_tr',num2str(loop),';']);
%             eval(['fid_g=fid_g',num2str(loop),';']);
            for loop=1:no_roi
%                 fprintf(fid_g,'%s\t%s\t%6.4e\t%6.4e\t%6.4e\t%6.4e\t%6.4e\t%6.4e\t%6.4e\t%6.4e\t%6.4e\t%6.4e\t%6.4e\t%6.4e\t%6.4e\t%6.4e\t%6.4e\t%6.4e\t%6.4e\t%6.4e\t%6.4e\t%6.4e\n',...
%                     p1,'p_value', pp);
                tem(row,1)={p1};
                tem(row,2)={'p_value'};
                for var_loop=1:length(pp)
                    tem(row,var_loop+2)={pp(var_loop)};
                end
                row=row+1;                
            end

        end % k
    end % j
end % if

% fprintf(fid_g1,'\n');
row=row+1;
for loop=1:no_roi
    for b_t=1:size(phrase1,1)
%         if bone_type=='V'
%             section_name_mean=[phrase1{b_t,:},bone_type,'_total mean'];
%         elseif bone_type=='F'
            section_name_mean=[phrase1{b_t,:},'_total mean'];
%         end
%         eval(['fid_g=fid_g',num2str(loop),';'])
        eval(['me=mean(mean_tr',num2str(loop),'_',num2str(b_t),',1,''omitnan'');'])
%         fprintf(fid_g1,'\t%s\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\n',...
%             section_name_mean, me);
        tem(row,2)={section_name_mean};
        for var_loop=1:length(me)
            tem(row,var_loop+2)={me(var_loop)};
        end
        row=row+1; 
    end
end

% eval(['xlswrite(''',out_file1,''',tem)'])
col_name1={'R_BS';      % 1
    'G_BS';             % 2
    'R_G';              % 3              %
    'Ronly_BS';        % 4
    'Gonly_BS';        % 5
    'sLS_BS';           % 6
    'dLS_BS';           % 7
    'LS_BS';            % 8
    'R_LS';             % 9             %
    'G_LS';             % 10             %
    'MS_BS';            % 11
    'sLS_LS';           % 12
    'dLS_LS';           % 13
    'dLS_sLS';          % 14
    'MAR_um';          % 15
    'BFR';              % 16
    'BV_TV';            % 17
    'Tb_Th';            % 18
    'Tb_N';             % 19
    'Tb_Sp'};          

Bone_type=tem(2:end,1);
Bone_no=tem(2:end,2);
for i=3:size(tem,2)
    eval(['name=cell2mat(col_name1(',num2str(i-2),',1));'])
    eval([name,'=cell2mat(tem(2:end-3,',num2str(i),'));'])
    eval([name,'(end+2:end+3)=cell2mat(tem(end-1:end,',num2str(i),'));'])    
end
T=table(Bone_type, Bone_no, R_BS,G_BS,R_G,Ronly_BS,Gonly_BS,sLS_BS,dLS_BS,LS_BS,R_LS,G_LS,MS_BS,sLS_LS,dLS_LS,dLS_sLS,MAR_um,BFR,BV_TV,Tb_Th,Tb_N,Tb_Sp);
eval(['writetable(T,''',out_file1,''')'])

% fclose all

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


% color={'red';
% 	   'blue';
%        'black';
%        'magenta';
%        'yellow';
%        'cyan';
%        'green'};
% 
% ma=0;
% for i=1:size(phrase1,1)
%     eval(['ma=max(ma,max(RSE_tr1_',num2str(i),'(:)));'])
% end
% 
% for loop=1:size(RSE_tr1_1,2)
%     plot(RSE_tr1_1(:,loop),color{1});
%     for i=2:size(phrase1,1)
%         hold on;
%         eval(['plot(RSE_tr1_',num2str(i),'(:,',num2str(loop),'), color{',num2str(i),'});'])
%     end
%     aaa=[];
%     for i=1:size(bone_number1,1)
%         aaa=[aaa,' ',bone_number1{i,1}];
%     end
%     ph=['legend([phrase2{1,:}, '' (',aaa,')'']'];
%     for j=2:size(phrase1,1)
%         eval(['a=bone_number',num2str(j),';'])
%         aaa=[];
%         for i=1:size(a,1)
%             aaa=[aaa,' ',a{i,1}];
%         end
%         ph=[ph,',[phrase2{',num2str(j),',:}, '' (',aaa,')'']'];
%     end
%     ph=[ph,');'];
%     eval(ph)
%     
%     xlabel('Bone number');
%     ylabel('RSE (%)');
%     title(['RSE of each bone for ',cell2mat(col_name1(loop))])
%     axis([1 size(RSE_tr1_1,1) 0 ma])
%     output_image_file=[direct,'RSE_Histo_',cell2mat(col_name2(loop))];
%     eval(['print -djpeg100 ''',output_image_file,''''])
%     close all
% end

function write_data7_cell(varargin)

% function write_data7_cell(direct, phrase1, phrase2, bone_type, delimeter)

% direct='Z:\KOMP\AAA_E16F\01_Submitted\Layers\Images\data\';
% phrase1={'AAA_E16_hF_F';'AAA_E16_hF_M'};
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
number_bones=varargin{6};
if nargin<7
    exclude_list=[0 0 0 0 0 0 0];
else
    exclude_list=varargin{7};
end
not_include=[];

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

%[9, 1;
%             22, 2;
%             23, 1];

%clear         
fclose all
section_effect=1;           % section_effect 1 : use mean value of sections of 1 mouce as 1 data
                            % section_effect 0 : use all of the values of sections of all mice
                            
dash=strfind(direct,delimeter);
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
    'AP_Intensity';          % 56
    'AP_Background_Intensity'; %57
    'TRAP_Intensity';        % 58
    'TRAP_Background_Intensity'; %59
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
    'AP_Intensity';          % 56
    'AP_Background_Intensity'; %57
    'TRAP_Intensity';        % 58
    'TRAP_Background_Intensity'; %59
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

% eval(['out_file1=''',direct,aaa(1:end-1),'analysis1_cell_tr1.txt'';'])
% fid_g1 = fopen(out_file1, 'w');
eval(['out_file1=''',direct,aaa(1:end-1),'analysis1_cell_tr1.xlsx'';'])
% if no_roi==2
%     eval(['out_file2=''',direct,aaa(1:end-1),'analysis1_cell_tr2.txt'';'])
%     fid_g2 = fopen(out_file2, 'w');
% end

tem(1)={'Bone type'};
tem(2)={'Bone #'};

if include_number_TRAP==0
%     fprintf(fid_g1,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n', ...
%         'Bone type','Bone #','AP_BS', 'AP_R_BS', 'AP_G_BS', 'AP_RG_BS', 'AP_R_RG_BS', 'AP_L_BS', 'AP_NL_BS', 'AP_L_AP', 'AP_NL_AP', 'AP_R_R', 'AP_G_G', 'AP_RG_RG', 'AP_R_RG_R', 'AP_TRAP_BS', 'TRAP_BS', 'TRAP_L_BS', 'TRAP_NL_BS', 'AP_TRAP_R_RG_BS', 'TRAP_on_TRAP', 'TRAP_L_TRAP_on', 'TRAP_NL_TRAP_on', 'GFP_BS', 'GFP_R_BS', 'GFPonly_BS', 'Height_Width', 'Osteocytes_N', 'Cortex_Width', 'AC_Intensity', 'RedBeads_Intensity', 'AC_RedBeads', 'AC_BR_Intensity', 'AC_BR_RedBeads', 'Calcein_Intensity', 'GreenBeads_Intensity', 'Calcein_GreenBeads');
%     if no_roi==2
%         fprintf(fid_g2,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n', ...
%         'Bone type','Bone #','AP_BS', 'AP_R_BS', 'AP_G_BS', 'AP_RG_BS', 'AP_R_RG_BS', 'AP_L_BS', 'AP_NL_BS', 'AP_L_AP', 'AP_NL_AP', 'AP_R_R', 'AP_G_G', 'AP_RG_RG', 'AP_R_RG_R', 'AP_TRAP_BS', 'TRAP_BS', 'TRAP_L_BS', 'TRAP_NL_BS', 'AP_TRAP_R_RG_BS', 'TRAP_on_TRAP', 'TRAP_L_TRAP_on', 'TRAP_NL_TRAP_on', 'GFP_BS', 'GFP_R_BS', 'GFPonly_BS', 'Height_Width', 'Osteocytes_N', 'Cortex_Width', 'AC_Intensity', 'RedBeads_Intensity', 'AC_RedBeads', 'AC_BR_Intensity', 'AC_BR_RedBeads', 'Calcein_Intensity', 'GreenBeads_Intensity', 'Calcein_GreenBeads');
%     end
    for k=1:length(col_name2)
        tem(k+2)=col_name2(k,:);
    end
elseif include_number_TRAP==1
%     fprintf(fid_g1,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n', ...
%         'Bone type','Bone #','AP_BS', 'AP_R_BS', 'AP_G_BS', 'AP_RG_BS', 'AP_R_RG_BS', 'AP_L_BS', 'AP_NL_BS', 'AP_L_AP', 'AP_NL_AP', 'AP_R_R', 'AP_G_G', 'AP_RG_RG', 'AP_R_RG_R', 'AP_TRAP_BS', 'TRAP_BS', 'TRAP_L_BS', 'TRAP_NL_BS', 'AP_TRAP_R_RG_BS', 'TRAP_on_TRAP', 'TRAP_L_TRAP_on', 'TRAP_NL_TRAP_on', 'GFP_BS', 'GFP_R_BS', 'GFPonly_BS', 'Height_Width', 'Osteocytes_N', 'Cortex_Width', 'AC_Intensity', 'RedBeads_Intensity', 'AC_RedBeads', 'AC_BR_Intensity', 'AC_BR_RedBeads', 'Calcein_Intensity', 'GreenBeads_Intensity', 'Calcein_GreenBeads', '# TRAP');
%     if no_roi==2
%         fprintf(fid_g2,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n', ...
%         'Bone type','Bone #','AP_BS', 'AP_R_BS', 'AP_G_BS', 'AP_RG_BS', 'AP_R_RG_BS', 'AP_L_BS', 'AP_NL_BS', 'AP_L_AP', 'AP_NL_AP', 'AP_R_R', 'AP_G_G', 'AP_RG_RG', 'AP_R_RG_R', 'AP_TRAP_BS', 'TRAP_BS', 'TRAP_L_BS', 'TRAP_NL_BS', 'AP_TRAP_R_RG_BS', 'TRAP_on_TRAP', 'TRAP_L_TRAP_on', 'TRAP_NL_TRAP_on', 'GFP_BS', 'GFP_R_BS', 'GFPonly_BS', 'Height_Width', 'Osteocytes_N', 'Cortex_Width', 'AC_Intensity', 'RedBeads_Intensity', 'AC_RedBeads', 'AC_BR_Intensity', 'AC_BR_RedBeads', 'Calcein_Intensity', 'GreenBeads_Intensity', 'Calcein_GreenBeads', '# TRAP');
%     end
    for k=1:length(col_name2)
        tem(k+2)=col_name2(k,:);
    end
    tem(k+3)={'# TRAP'};
end
row=2;

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Remove 9th bone, if empty
%
%   Dec. 23 2020
%
        if isempty(bone_number{9,2}) | number_bones(1)==8
            bone_number=bone_number(1:8,:);
        end
%
%   Dec. 23 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        for b_n=1:size(bone_number,1)
            if b_n==7
                b_n;
            end

            if ana_n==0
%                 phr1_1=[phrase1{b_t,:},'_',num2str(bone_number(b_n)),'_s1_analysis_tr',num2str(loop),'.mat'];
%                 phr1_2=[phrase1{b_t,:},'_',num2str(bone_number(b_n)),'_s2_analysis_tr',num2str(loop),'.mat'];
%                 phr1_3=[phrase1{b_t,:},'_',num2str(bone_number(b_n)),'_s3_analysis_tr',num2str(loop),'.mat'];
                phr1_1=[phrase1{b_t,:},'L1', '_s',bone_number{b_n},'_analysis_tr',num2str(loop),'.mat'];
                phr1_2=[phrase1{b_t,:},'L2', '_s',bone_number{b_n},'_analysis_tr',num2str(loop),'.mat'];
                phr1_3=[phrase1{b_t,:},'L3', '_s',bone_number{b_n},'_analysis_tr',num2str(loop),'.mat'];
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

%                 if isempty(not_include)==0
%                     temp=find(not_include(:,1)==bone_number(b_n));
%                     if isempty(temp)==0
%                         if isempty(find(not_include(temp,2)==sample(s)))==0
% %                             data=[0,0,0,0,0,0,0,0,0];
%                             data=zeros(1,size(col_name1,1)); %[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
%                             eval(['bone_tr',num2str(loop),'_',num2str(b_t),'(number_tr',num2str(loop),'_',num2str(b_t),'+1,:)=data;'])
%                             eval(['bone_tr',num2str(loop),'_',num2str(b_t),'_c(number_tr',num2str(loop),'_',num2str(b_t),'_c+1,:)=data;'])
%                             eval(['number_tr',num2str(loop),'_',num2str(b_t),'=number_tr',num2str(loop),'_',num2str(b_t),'+1;'])
%                             eval(['number_tr',num2str(loop),'_',num2str(b_t),'_c=number_tr',num2str(loop),'_',num2str(b_t),'_c+1;'])
%                             eval(['number_tr',num2str(loop),'_',num2str(b_t),'_c1(',num2str(b_n),')=number_tr',num2str(loop),'_',num2str(b_t),'_c1(',num2str(b_n),')+1;'])
%                             eval(['no',num2str(b_t),'.ap(',num2str(b_n),',',num2str(s),')=0;'])
%                             eval(['no',num2str(b_t),'.tr(',num2str(nnn),',',num2str(s),')=0;'])
%                             continue
%                         end
%                     end
%                 end

                    
%                 eval(['aa=test',num2str(sample(s)),';'])
                eval(['aa=test',num2str(s),';'])
                temp_list=exclude_list(exclude_list(:,1)==b_t,:);
                temp_list=temp_list(temp_list(:,2)==str2num(nnn),:);
                temp_list=temp_list(temp_list(:,3)==s,:);

                if isempty(aa)==1 %| ~isempty(temp_list(temp_list(:,4)==1,:))
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
                        if isfield(a, 'intensity')
                            data=[a.surface.A_BS, a.surface.rA_BS, a.surface.gA_BS, a.surface.rgA_BS, a.surface.rA_BS+a.surface.rgA_BS, LA_BS, a.surface.A_only_BS, LA_BS/a.surface.A_BS*100, a.surface.A_only_BS/a.surface.A_BS*100, a.surface.rA_BS/a.surface.r_BS*100, a.surface.gA_BS/a.surface.g_BS*100, a.surface.rgA_BS/a.surface.rg_BS*100, (a.surface.rA_BS+a.surface.rgA_BS)/a.surface.r_BS*100, a.surface.AT_BS, a.surface.T_BS, LT_BS, a.surface.T_only_BS, (a.surface.rAT_BS+a.surface.rgAT_BS),...
                                a.volume.T_on_T, LT_BS/a.surface.T_BS*100, a.surface.T_only_BS/a.surface.T_BS*100, a.surface.GFP_BS, a.surface.rGFP_BS, a.surface.GFP_only_BS, a.roi.ratio, a.cortex.osteocytes.number_per_cortex_area*10^6, a.cortex.thickness, a.intensity.AC, a.intensity.Red_beads, a.intensity.AC_Red_beads, a.intensity.AC_BR, a.intensity.AC_BR_Red_beads, a.intensity.Calcein, a.intensity.Green_beads, a.intensity.Calcein_Green_beads, a.intensity.AP, a.intensity.AP_Background_mean, a.intensity.TRAP, a.intensity.TRAP_Background_mean];
                        else
                            data=[a.surface.A_BS, a.surface.rA_BS, a.surface.gA_BS, a.surface.rgA_BS, a.surface.rA_BS+a.surface.rgA_BS, LA_BS, a.surface.A_only_BS, LA_BS/a.surface.A_BS*100, a.surface.A_only_BS/a.surface.A_BS*100, a.surface.rA_BS/a.surface.r_BS*100, a.surface.gA_BS/a.surface.g_BS*100, a.surface.rgA_BS/a.surface.rg_BS*100, (a.surface.rA_BS+a.surface.rgA_BS)/a.surface.r_BS*100, a.surface.AT_BS, a.surface.T_BS, LT_BS, a.surface.T_only_BS, (a.surface.rAT_BS+a.surface.rgAT_BS),...
                                a.volume.T_on_T, LT_BS/a.surface.T_BS*100, a.surface.T_only_BS/a.surface.T_BS*100, a.surface.GFP_BS, a.surface.rGFP_BS, a.surface.GFP_only_BS, a.roi.ratio, a.cortex.osteocytes.number_per_cortex_area*10^6, a.cortex.thickness, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN];
                        end
                    elseif include_number_TRAP==1
                        if isfield(a, 'intensity')
                            data=[a.surface.A_BS, a.surface.rA_BS, a.surface.gA_BS, a.surface.rgA_BS, a.surface.rA_BS+a.surface.rgA_BS, LA_BS, a.surface.A_only_BS, LA_BS/a.surface.A_BS*100, a.surface.A_only_BS/a.surface.A_BS*100, a.surface.rA_BS/a.surface.r_BS*100, a.surface.gA_BS/a.surface.g_BS*100, a.surface.rgA_BS/a.surface.rg_BS*100, (a.surface.rA_BS+a.surface.rgA_BS)/a.surface.r_BS*100, a.surface.AT_BS, a.surface.T_BS, LT_BS, a.surface.T_only_BS, (a.surface.rAT_BS+a.surface.rgAT_BS),...
                                a.volume.T_on_T, LT_BS/a.surface.T_BS*100, a.surface.T_only_BS/a.surface.T_BS*100, a.surface.GFP_BS, a.surface.rGFP_BS, a.surface.GFP_only_BS, a.roi.ratio, a.cortex.osteocytes.number_per_cortex_area*10^6, a.cortex.thickness, a.intensity.AC, a.intensity.Red_beads, a.intensity.AC_Red_beads, a.intensity.AC_BR, a.intensity.AC_BR_Red_beads, a.intensity.Calcein, a.intensity.Green_beads, a.intensity.Calcein_Green_beads, a.intensity.AP, a.intensity.AP_Background_mean, a.intensity.TRAP, a.intensity.TRAP_Background_mean, a.volume_T_number];
                        else
                            data=[a.surface.A_BS, a.surface.rA_BS, a.surface.gA_BS, a.surface.rgA_BS, a.surface.rA_BS+a.surface.rgA_BS, LA_BS, a.surface.A_only_BS, LA_BS/a.surface.A_BS*100, a.surface.A_only_BS/a.surface.A_BS*100, a.surface.rA_BS/a.surface.r_BS*100, a.surface.gA_BS/a.surface.g_BS*100, a.surface.rgA_BS/a.surface.rg_BS*100, (a.surface.rA_BS+a.surface.rgA_BS)/a.surface.r_BS*100, a.surface.AT_BS, a.surface.T_BS, LT_BS, a.surface.T_only_BS, (a.surface.rAT_BS+a.surface.rgAT_BS),...
                                a.volume.T_on_T, LT_BS/a.surface.T_BS*100, a.surface.T_only_BS/a.surface.T_BS*100, a.surface.GFP_BS, a.surface.rGFP_BS, a.surface.GFP_only_BS, a.roi.ratio, a.cortex.osteocytes.number_per_cortex_area*10^6, a.cortex.thickness, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, a.volume_T_number];
                        end
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
                        if  ~isempty(temp_list(temp_list(:,5)==1,:))     % exclude AP
                            eval(['bone_tr',num2str(loop),'_',num2str(b_t),'(number_tr',num2str(loop),'_',num2str(b_t),',1:14)=ones(1,14)*NaN;'])
                            eval(['number_tr',num2str(loop),'_',num2str(b_t),'=number_tr',num2str(loop),'_',num2str(b_t),';'])
                            eval(['bone_tr',num2str(loop),'_',num2str(b_t),'(number_tr',num2str(loop),'_',num2str(b_t),',18)=ones(1,1)*NaN;'])
                            eval(['number_tr',num2str(loop),'_',num2str(b_t),'=number_tr',num2str(loop),'_',num2str(b_t),';'])
                        elseif ~isempty(temp_list(temp_list(:,6)==1,:)) % exclude TRAP
                            eval(['bone_tr',num2str(loop),'_',num2str(b_t),'(number_tr',num2str(loop),'_',num2str(b_t),',14:21)=ones(1,8)*NaN;'])
                            eval(['number_tr',num2str(loop),'_',num2str(b_t),'=number_tr',num2str(loop),'_',num2str(b_t),';'])
                        elseif ~isempty(temp_list(temp_list(:,4)==1,:)) % exclude All
                            eval(['bone_tr',num2str(loop),'_',num2str(b_t),'(number_tr',num2str(loop),'_',num2str(b_t),',:)=ones(1,no_column)*NaN;'])
                            eval(['bone_tr',num2str(loop),'_',num2str(b_t),'_c(number_tr',num2str(loop),'_',num2str(b_t),'_c,:)=ones(1,no_column)*NaN;'])
                            eval(['number_tr',num2str(loop),'_',num2str(b_t),'=number_tr',num2str(loop),'_',num2str(b_t),';'])
                            eval(['number_tr',num2str(loop),'_',num2str(b_t),'_c=number_tr',num2str(loop),'_',num2str(b_t),'_c;'])
                            eval(['number_tr',num2str(loop),'_',num2str(b_t),'_c1(',num2str(b_n),')=number_tr',num2str(loop),'_',num2str(b_t),'_c1(',num2str(b_n),');'])
                        end
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Remove 9th bone, if empty
%
%   Dec. 23 2020
%
        if isempty(bone_number{9,2}) | number_bones(1)==8
            bone_number=bone_number(1:8,:);
        end
%
%   Dec. 23 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
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
                    if exist('ind','var')
                        data_for_average=double(data_for_average(find(ind),:));
                    end
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
                    phras=[phras,'=mean(data_for_average1,1,''omitnan'');'];
                    eval(phras)
%                     phras=['std_tr',num2str(loop),'_',num2str(b_t),'(',num2str(b_n),',1:end_column)'];
%                     phras=[phras,'=nanstd(data_for_average1,0);'];
%                     eval(phras)

                    phras=['std_tr',num2str(loop),'_',num2str(b_t),'(',num2str(b_n),',1:end_column)'];
                    phras=[phras,'=std(data_for_average1,0,''omitnan'');'];
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
                    if exist('ind','var')
                        data_for_average=double(data_for_average(find(ind),:));
                    end
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
                    phras=[phras,'=mean(data_for_average1,1,''omitnan'');'];
                    eval(phras);
%                     phras=['std_tr',num2str(loop),'_',num2str(b_t),'(',num2str(b_n),',start_column:end_column)'];
%                     phras=[phras,'=nanstd(data_for_average1,0);'];
%                     eval(phras);

                    phras=['std_tr',num2str(loop),'_',num2str(b_t),'(',num2str(b_n),',start_column:end_column)'];
                    phras=[phras,'=std(data_for_average1,0,''omitnan'');'];
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
                    end_column=39;
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
                    if exist('ind','var')
                        data_for_average=double(data_for_average(find(ind),:));
                    end
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
                    phras=[phras,'=mean(data_for_average1,1,''omitnan'');'];
                    eval(phras);
%                     phras=['std_tr',num2str(loop),'_',num2str(b_t),'(',num2str(b_n),',start_column:end_column)'];
%                     phras=[phras,'=nanstd(data_for_average1,0);'];
%                     eval(phras);

                    phras=['std_tr',num2str(loop),'_',num2str(b_t),'(',num2str(b_n),',start_column:end_column)'];
                    phras=[phras,'=std(data_for_average1,0,''omitnan'');'];
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
        b_t;
        eval(['bone_number=bone_number',num2str(b_t),';']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Remove 9th bone, if empty
%
%   Dec. 23 2020
%
        if isempty(bone_number{9,2}) | number_bones(1)==8
            bone_number=bone_number(1:8,:);
        end
%
%   Dec. 23 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        count=0;
        minus_count(b_t,length(bone_number))=0;
        
        for b_n=1:size(bone_number,1)
            if b_n==3
                b_n;
            end
            nnn=bone_number{b_n};

            if ana_n==0
%                 phr1_1=[phrase1{b_t,:},num2str(bone_number(b_n)),'_1_analysis_tr',num2str(loop),'.mat'];
%                 phr1_2=[phrase1{b_t,:},num2str(bone_number(b_n)),'_2_analysis_tr',num2str(loop),'.mat'];
%                 phr1_3=[phrase1{b_t,:},num2str(bone_number(b_n)),'_3_analysis_tr',num2str(loop),'.mat'];
                phr1_1=[phrase1{b_t,:},'L1', '_s',bone_number{b_n},'_analysis_tr',num2str(loop),'.mat'];
                phr1_2=[phrase1{b_t,:},'L2', '_s',bone_number{b_n},'_analysis_tr',num2str(loop),'.mat'];
                phr1_3=[phrase1{b_t,:},'L3', '_s',bone_number{b_n},'_analysis_tr',num2str(loop),'.mat'];
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
                eval(['d=bone_tr',num2str(loop),'_',num2str(b_t),'(',num2str(count),'+',num2str(s),'-sum(minus_count(',num2str(b_t),',:)),:);'])                
                if loop==1
%                     if (s==1 & b_n==1) %| samples_per_bone==1
%                         ph=['fprintf(fid_g1,''%s\t%s\'];
%                         for col=1:no_column
%                             ph=[ph, 't%8.6f\'];
%                         end
%                         ph=[ph, 'n'', phrase2{b_t,:}, section_name, bone_tr1_',num2str(b_t),'(count+s-sum(minus_count(b_t,:)),:));'];
%                         eval(ph);
%                     else
%                         ph=['fprintf(fid_g1,''\t%s\'];
%                         for col=1:no_column
%                             ph=[ph, 't%8.6f\'];
%                         end
%                         ph=[ph, 'n'', section_name, bone_tr1_',num2str(b_t),'(count+s-sum(minus_count(b_t,:)),:));'];
%                         eval(ph);
%                     end
                    if (s==1 & b_n==1) %| samples_per_bone(b_t, b_n)==1
                        tem(row,1)={phrase2{b_t,:}};
                    end
                    tem(row,2)={section_name};
                    for var_loop=1:length(d)
                        tem(row,var_loop+2)={d(var_loop)};
                    end
                    row=row+1;
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
            
            eval(['d_mean=mean_tr',num2str(loop),'_',num2str(b_t),'(',num2str(b_n),',:);'])
            eval(['d_std=std_tr',num2str(loop),'_',num2str(b_t),'(',num2str(b_n),',:);'])
            eval(['d_RSE=RSE_tr',num2str(loop),'_',num2str(b_t),'(',num2str(b_n),',:);'])
%             ph=['fprintf(fid_g',num2str(loop),',''\t%s\'];
%             for col=1:no_column
%                 ph=[ph, 't%8.6f\'];
%             end
% %                ph=[ph, 'n'', ''mean'', mean_tr1_',num2str(b_t),'(b_n,:));'];
%             ph=[ph, 'n'', ''',section_name_mean,''', mean_tr1_',num2str(b_t),'(b_n,:));'];
%             eval(ph);
%             ph=['fprintf(fid_g',num2str(loop),',''\t%s\'];
%             for col=1:no_column
%                 ph=[ph, 't%8.6f\'];
%             end
% %                ph=[ph, 'n'', ''std'', std_tr1_',num2str(b_t),'(b_n,:));'];
%             ph=[ph, 'n'', ''',section_name_std,''', std_tr1_',num2str(b_t),'(b_n,:));'];
%             eval(ph);
%             ph=['fprintf(fid_g',num2str(loop),',''\t%s\'];
%             for col=1:no_column
%                 ph=[ph, 't%8.6f\'];
%             end
% %                ph=[ph, 'n'', ''RSE(%)'', RSE_tr1_',num2str(b_t),'(b_n,:));'];
%             ph=[ph, 'n'', ''',section_name_RSE,''', RSE_tr1_',num2str(b_t),'(b_n,:));'];
%             eval(ph);

            tem(row,2)={section_name_mean};
            for var_loop=1:length(d_mean)
                tem(row,var_loop+2)={d_mean(var_loop)};
            end
            row=row+1;
            tem(row,2)={section_name_std};
            for var_loop=1:length(d_std)
                tem(row,var_loop+2)={d_std(var_loop)};
            end
            row=row+1;                
            tem(row,2)={section_name_RSE};
            for var_loop=1:length(d_RSE)
                tem(row,var_loop+2)={d_RSE(var_loop)};
            end
            row=row+1;
            
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

        eval(['data_for_std=mean_tr',num2str(loop),'_',num2str(b_t),';'])
        for ppp=1:size(data_for_std,2)
            temp=(strfind(num2str((data_for_std(:,ppp))'),'NaN'))
            number_bone_sample(b_t, ppp)=size(data_for_std,1)-length(temp);
        end
%         average=sum(data_for_std,1,'omitnan')./number_bone_sample(b_t,:);
%         standard_deviation=sqrt(sum((data_for_std-ones(size(data_for_std,1),1,'omitnan')*average).^2)./(number_bone_sample(b_t,:)-1));
        average=mean(data_for_std,1,'omitnan');
        standard_deviation=std(data_for_std,1,'omitnan');
        RSE=standard_deviation./average/sqrt(number_bone_sample(b_t))*100;

%         eval(['fid_g=fid_g',num2str(loop),';'])
%         ph=['fprintf(fid_g,''\t%s\'];
%         for col=1:no_column
%             ph=[ph, 't%8.6f\'];
%         end
%         ph=[ph, 'n'', ''',section_name_mean,''', average);'];
%         eval(ph);
%         ph=['fprintf(fid_g,''\t%s\'];
%         for col=1:no_column
%             ph=[ph, 't%8.6f\'];
%         end
%         ph=[ph, 'n'', ''',section_name_std,''', standard_deviation);'];
%         eval(ph);
%         ph=['fprintf(fid_g,''\t%s\'];
%         for col=1:no_column
%             ph=[ph, 't%8.6f\'];
%         end
%         ph=[ph, 'n'', ''',section_name_RSE,''', RSE);'];
%         eval(ph);

        tem(row,2)={section_name_mean};
        for var_loop=1:length(average)
            tem(row,var_loop+2)={average(var_loop)};
        end
        row=row+1;
        tem(row,2)={section_name_std};
        for var_loop=1:length(standard_deviation)
            tem(row,var_loop+2)={standard_deviation(var_loop)};
        end
        row=row+1;                
        tem(row,2)={section_name_RSE};
        for var_loop=1:length(RSE)
            tem(row,var_loop+2)={RSE(var_loop)};
        end
        row=row+1;
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
%                 eval(['fid_g=fid_g',num2str(loop),';']);
%                 if loop==1
%                     ph=['fprintf(fid_g,''%s\t%s\'];
%                     for col=1:no_column
%                         ph=[ph, 't%4.2e\'];
%                     end
 %                   ph=[ph, 'n'', ''t-test'', ''p_value'', p_b_c_tr1);'];
%                     ph=[ph, 'n'', ''',p1,''', ''p_value'', pp);'];
%                     eval(ph);
%                     ph=['fprintf(fid_g,''%s\t\'];
%                     for col=1:no_column
%                         ph=[ph, 't%d\'];
%                     end
%                     ph=[ph, 'n'', ''',p1,''', hh);'];
%                     eval(ph);
                tem(row,1)={p1};
                tem(row,2)={'p_value'};
                for var_loop=1:length(pp)
                    tem(row,var_loop+2)={pp(var_loop)};
                end
                row=row+1; 

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

%fprintf(fid_g1,'\n');
row=row+1;
for loop=1:no_roi
    for b_t=1:size(phrase1,1)   
%         if bone_type=='V'
%             section_name_mean=[phrase1{b_t,:},bone_type,'_total mean'];
%         elseif bone_type=='F'
            section_name_mean=[phrase1{b_t,:},'_total mean'];
%         end
        eval(['data_for_std=mean_tr',num2str(loop),'_',num2str(b_t),';'])
        for ppp=1:size(data_for_std,2)
            temp=(strfind(num2str((data_for_std(:,ppp))'),'NaN'));
            number_bone_sample(b_t, ppp)=size(data_for_std,1)-length(temp);
        end
        average=sum(data_for_std,1,'omitnan')./number_bone_sample(b_t,:);

%         eval(['fid_g=fid_g',num2str(loop),';'])
%         ph=['fprintf(fid_g,''\t%s\'];
%         for col=1:no_column
%             ph=[ph, 't%8.6f\'];
%         end
%         ph=[ph, 'n'', ''',section_name_mean,''', average);'];
%         eval(ph);
        tem(row,2)={section_name_mean};
        for var_loop=1:length(average)
            tem(row,var_loop+2)={average(var_loop)};
        end
        row=row+1; 
    end
end
% eval(['xlswrite(''',out_file1,''',tem)'])
      

Bone_type=tem(2:end,1);
Bone_no=tem(2:end,2);
phase=['T=table(Bone_type, Bone_no,'];
for i=3:size(tem,2)
    eval(['name=cell2mat(col_name2(',num2str(i-2),',1));'])
    eval([name,'=cell2mat(tem(2:end-3,',num2str(i),'));'])
    eval([name,'(end+2:end+3)=cell2mat(tem(end-1:end,',num2str(i),'));'])
    phase=[phase,name,','];
end
phase=[phase(1:end-1),');'];
eval(phase);
eval(['writetable(T,''',out_file1,''')'])

% fclose all

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


% color={'red';
% 	   'blue';
%        'black';
%        'magenta';
%        'yellow';
%        'cyan';
%        'green'};
% 
% ma=0;
% for i=1:size(phrase1,1)
%     eval(['ma=max(ma,max(RSE_tr1_',num2str(i),'(:)));'])
% end
% ma=ma*1.2;
% 
% for loop=1:size(RSE_tr1_1,2)
%     plot(RSE_tr1_1(:,loop),color{1});
%     for i=2:size(phrase1,1)
%         hold on;
%         eval(['plot(RSE_tr1_',num2str(i),'(:,',num2str(loop),'), color{',num2str(i),'});'])
%     end
%     aaa=[];
%     for i=1:size(bone_number1,1)
%         aaa=[aaa,' ',bone_number1{i,1}];
%     end
%     ph=['legend([phrase2{1,:}, '' (',aaa,')'']'];
%     for j=2:size(phrase1,1)
%         eval(['a=bone_number',num2str(j),';'])
%         aaa=[];
%         for i=1:size(a,1)
%             aaa=[aaa,' ',a{i,1}];
%         end
%         ph=[ph,',[phrase2{',num2str(j),',:}, '' (',aaa,')'']'];
%     end
%     ph=[ph,');'];
%     eval(ph)
%     
%     xlabel('Bone number');
%     ylabel('RSE (%)');
%     title(['RSE of each bone for ',cell2mat(col_name1(loop))])
%     axis([1 size(RSE_tr1_1,1) 0 ma])
%     output_image_file=[direct,'RSE_cell_',cell2mat(col_name2(loop))];
%     eval(['print -djpeg100 ''',output_image_file,''''])
%     close all
% end

function read_excel_data_4(varargin)

%
%   read_excel_data_4(bone_excel_file, cell_excel_file, output_excel_file, bone_n, margin, use_mean_median)

%   bone_excel_file     : input excel file (bone histomorphometry) with possible outliers
%   cell_excel_file     : input excel file (cellular activity) with possible outliers
%   output_excel_file   : output excle file name
%   bone_n              : number of sections per bone type (b_t)  ex) [8,8]
%                           should be [8 8] % changed October 17, 2014
%   margin              : outlier boundary (default : 1.5 SD)

%   use_mean_median     : 0 --> use median (default)
%                         1 --> use mean

% clear
% [data, txtdata]= xlsread('Z:\Histomorphometry\His18\01_Submitted\Layers\Images\His18_3.6129F_analysis1_tr1_all.xls');
% 


bone_excel_file=varargin{1};
cell_excel_file=varargin{2};
output_excel_file=varargin{3};
% b_t=2;                                          % number of comparison
% bone_n=[8,8];                                   % number of bone section per bone group (b_t)
bone_n=varargin{4};
b_t=length(bone_n);

if nargin==4
    margin=1.5;
    use_mean_median=0;
elseif nargin==5
    margin=varargin{5};
    use_mean_median=0;
elseif nargin==6
    margin=varargin{5};
    use_mean_median=varargin{6};
else
    error('Error : Too many argments');
end

[data, txtdata]= xlsread(bone_excel_file);
[data_cell, txtdata_cell]=xlsread(cell_excel_file);

diff_size=size(txtdata,1)-size(data,1);
if diff_size>1
    data(end+1:end+diff_size-1,:)=NaN*ones(diff_size-1,size(data,2));
    data_cell(end+1:end+diff_size-1,:)=NaN*ones(diff_size-1,size(data_cell,2));
end
data=[data,data_cell];
txtdata=[txtdata,txtdata_cell(:,3:end)];
ttt=txtdata;
for i=1:size(data,1)
    for j=1:size(data,2)
%         if ~isnan(data(i,j))
            ttt(i+1,j+2)={data(i,j)};
%         end
    end
end
% ttt=ttt(1:end-(b_t+1),:);
ttt=ttt(1:end-(b_t),:);
% % p_pos=size(ttt,1)-1;
% count2=0;
% for i=1:size(ttt,1)
%     c=strfind(ttt{i,1},'t-test : ');
%     if ~isempty(c)
%         count2=count2+1;
%         ind_p(count2)=i;
%     end
% end

% margin=1.5;                                     % SD(standard deviation) search margin (2.0 : 95.5%, 1.5 : 86.6%, 1 : 68.3%), formula --> erf(z/sqrt(2))
% use_mean_median=1;                              % 0 : use median (default), 1: use mean

% major=[1, 2, 12, 14, 15, 18, 21, 24];                           % R/BS, G/BS, MAR, BV/TV, Tb.Th, AP/BS, GFP/BS, TRAP/BS
major=[1, 2, 15, 17, 18, 21, 35, 42, 45, 46, 47, 48, 53, 56, 57, 58, 59];         % R/BS, G/BS, MAR, BV/TV, Tb.Th, AP/BS, TRAP/BS, GFP/BS, Height/Width, Osteocytes_n, Cortex_width, AC_Intensity, Calcein_Intensity, AP_Intensity, AP_Background_Intensity, TRAP_Intensity, TRAP_Background_Intensity

% column=['CDEFGHIJKLMNOPQRSTUVWXYZ'];
col_name1={'R_BS';                  % 1
    'G_BS';                         % 2
    'R_G';                          % 3
    'Ronly_BS';                    % 4
    'Gonly_BS';                    % 5
    'sLS_BS';                       % 6
    'dLS_BS';                       % 7
    'LS_BS';                        % 8
    'R_LS';                         % 9
    'G_LS';                         % 10
    'MS_BS';                        % 11
    'sLS_LS';                       % 12
    'dLS_LS';                       % 13
    'dLS_sLS';                      % 14
    'MAR_um';                      % 15
    'BFR';                          % 16
    'BV_TV';                        % 17
    'Tb.Th';                        % 18
    'Tb.N';                         % 19
    'Tb.Sp';                        % 20
    'AP_BS';                        % 21
    'AP_R_BS';                      % 22
    'AP_G_BS';                      % 23
    'AP_RG_BS';                     % 24
    'AP_R_RG_BS';                   % 25
    'AP_L_BS';                      % 26
    'AP_NL_BS';                     % 27
    'AP_L_AP';                      % 28
    'AP_NL_AP';                     % 29
    'AP_R_R';                       % 30
    'AP_G_G';                       % 31
    'AP_RG_RG';                     % 32
    'AP_R_RG_R';                    % 33
    'AP_TRAP_BS';                   % 34
    'TRAP_BS';                      % 35
    'TRAP_L_BS';                    % 36
    'TRAP_NL_BS';                   % 37
    'AP_TRAP_R_RG_BS';              % 38
    'TRAP_on_TRAP';                 % 39
    'TRAP_L_TRAP_on';               % 40
    'TRAP_NL_TRAP_on';              % 41
    'GFP_BS';                       % 42
    'GFP_R_BS';                     % 43
    'GFPonly_BS';                   % 44
    'Height_Width';                 % 45
    'Osteocytes_Density';           % 46
    'Cortex_Width';                 % 47
    'AC_Intensity';                 % 48
    'RedBeads_Intensity';           % 49
    'AC_RedBeads';                  % 50
    'AC_BR_Intensity';              % 51
    'AC_BR_RedBeads';               % 52
    'Calcein_Intensity';            % 53
    'GreenBeads_Intensity';         % 54
    'Calcein_GreenBeads';           % 55
    'AP_Intensity';                 % 56
    'AP_Background_Intensity';      % 57
    'TRAP_Intensity';               % 58
    'TRAP_Background_Intensity';    % 59
    };

vari=length(col_name1);

for bt=1:b_t
    b_n=bone_n(bt);
    for bn=1:b_n
        st_n=(bn-1)*6+1;
        if bt>1
%             st_t=(bt-1)*(6*bone_n(bt-1)+3)+1;
            st_t=sum(bone_n(1:bt-1))*6+(bt-1)*3+1;
        else
            st_t=1;
        end
        st=st_t+st_n-1;
        eval(['bone_',num2str(bt),'_',num2str(bn),'(1:6,:)=data(',num2str(st),':',num2str(st+5),',:);'])
        eval(['section_',num2str(bt),'(',num2str(bn),',:)=sum(bone_',num2str(bt),'_',num2str(bn),'(1:3,:)>0);'])        
    end
end

for bt=1:b_t
    t_data=[];
    b_n=bone_n(bt);
    for bn=1:b_n
        eval(['temp_data=bone_',num2str(bt),'_',num2str(bn),';'])
        t_data=[t_data;temp_data(1:3,:)];
    end
%     if use_mean_median==1
%         max_range=nanmean(t_data)+nanstd(t_data)*margin;
%         min_range=nanmean(t_data)-nanstd(t_data)*margin;
%     elseif use_mean_median==0
%         max_range=nanmedian(t_data)+nanstd(t_data)*margin;
%         min_range=nanmedian(t_data)-nanstd(t_data)*margin;
%     end
    if use_mean_median==1
        for i=1:size(t_data,2)
            tttt=t_data(:,i);
            max_range(i)=mean(tttt(tttt>0),'omitnan')+std(tttt(tttt>0),'omitnan')*margin;
            min_range(i)=mean(tttt(tttt>0),'omitnan')-std(tttt(tttt>0),'omitnan')*margin;
        end
    elseif use_mean_median==0
        for i=1:size(t_data,2)
            tttt=t_data(:,i);
            max_range(i)=median(tttt(tttt>0),'omitnan')+std(tttt(tttt>0),'omitnan')*margin;
            min_range(i)=median(tttt(tttt>0),'omitnan')-std(tttt(tttt>0),'omitnan')*margin;
        end
    end
    b_n=bone_n(bt);
    for bn=1:b_n
        for i=1:3
            eval(['t_data=bone_',num2str(bt),'_',num2str(bn),';'])
            ma_ind=find((t_data(i,major)-max_range(1,major))>0);
            mi_ind=find((t_data(i,major)-min_range(1,major))<0);
            temp=sort(unique([ma_ind, mi_ind]));
            if isempty(temp)
                eval(['out_data_',num2str(bt),'_',num2str(bn),'(i,:)=t_data(i,:);'])
                eval(['outliers_',num2str(bt),'_',num2str(bn),'_',num2str(i),'={''''};'])
            else
                new_data=t_data(i,:);
                if ~isempty(temp)
                    for loop=1:length(temp)
                        temp1=temp(loop);                   
                        if temp1==1         % R/BS
%                             new_data([1,3,5:11,13, 19, 20, 22, 23, 25, 26])=NaN;
%                             new_data([1, 3, 5:11, 13, 19, 21:23, 25, 28:29, 33:34])=NaN;
                            new_data([1, 3:4, 6:14, 16, 22, 24:30, 32:33, 36:38, 40:41, 43:44])=NaN;
                        elseif temp1==2     % G/BS
%                             new_data([2,4:11,13])=NaN;
%                             new_data([2, 4:11,13, 20:22, 24:24])=NaN;
                            new_data([2:3, 5:14, 16, 23:29, 31:33, 36:38, 40:41, 44])=NaN;
                        elseif temp1==4     % BV/TV
%                             new_data([14,16:17])=NaN;
                            new_data([17,19:20])=NaN;
                        elseif temp1==5     % Tb.Th
%                             new_data(15:17)=NaN;
                            new_data(18:20)=NaN;
                        elseif temp1==3     % MAR
%                             new_data(12:13)=NaN;
                            new_data(15:16)=NaN;
                        elseif temp1==6     % AP/BS
%                             new_data(18:26)=NaN;
                            new_data([21:34, 38])=NaN;
                        elseif temp1==7     % TRAP/BS
%                             new_data(26:31)=NaN;
                            new_data(35:41)=NaN;
                        elseif temp1==8     % GFP/BS
%                             new_data(32:34)=NaN;
                            new_data(42:44)=NaN;
                        elseif temp1==9     % Height/Width
%                             new_data(32:34)=NaN;
                            new_data(45)=NaN;
                        elseif temp1==10     % Osteocytes_N/Cortex_area
%                             new_data(32:34)=NaN;
                            new_data(46)=NaN;
                        elseif temp1==11     % Cortex_Width
%                             new_data(32:34)=NaN;
                            new_data(47)=NaN;
                        elseif temp1==12     % AC_Intensity
%                             new_data(32:34)=NaN;
                            new_data([48,50:52])=NaN;
                        elseif temp1==13     % Calcein_Intensity
%                             new_data(32:34)=NaN;
                            new_data([53,55])=NaN;
                        elseif temp1==14     % AP_Intensity
%                             new_data(32:34)=NaN;
                            new_data([56])=NaN;
                        elseif temp1==15     % AP_Background_Intensity
%                             new_data(32:34)=NaN;
                            new_data([57])=NaN;
                        elseif temp1==16     % TRAP_Intensity
%                             new_data(32:34)=NaN;
                            new_data([58])=NaN;
                        elseif temp1==17     % TRAP_Background_Intensity
%                             new_data(32:34)=NaN;
                            new_data([59])=NaN;
                        end
                    end
                end
                eval(['out_data_',num2str(bt),'_',num2str(bn),'(i,:)=new_data;'])
                eval(['outliers_',num2str(bt),'_',num2str(bn),'_',num2str(i),'={major(temp)};'])
            end 
        end
    end
end

for bt=1:b_t
    eval(['no_',num2str(bt),'=zeros(1,',num2str(vari),');'])
    eval(['mice_no_',num2str(bt),'=zeros(1,',num2str(vari),');'])
    b_n=bone_n(bt);
    for bn=1:b_n
        eval(['temp=out_data_',num2str(bt),'_',num2str(bn),';'])
        temp1=temp; temp1(isnan(temp))=0;
        number=sum(~isnan(temp));
        eval(['no_',num2str(bt),'=no_',num2str(bt),'+number;'])
        temp(4,:)=mean(temp(1:3,:),'omitnan');
        mice_number=~isnan(temp(4,:));
        eval(['mice_no_',num2str(bt),'=mice_no_',num2str(bt),'+mice_number;'])        
        tt=std(temp(1:3,:),'omitnan'); tt(tt==0)=NaN;
        temp(5,:)=tt;
        temp(6,:)=temp(5,:)./temp(4,:)./sqrt(sum(temp1~=0))*100;
        eval(['out_data_',num2str(bt),'_',num2str(bn),'=temp;'])
        eval(['total_data_',num2str(bt),'(',num2str(bn),',:)=temp(4,:);'])
    end
end

for bt=1:b_t
    eval(['temp=total_data_',num2str(bt),';'])
    temp1=temp; temp1(isnan(temp))=0;
    eval(['total_',num2str(bt),'(1,:)=mean(temp,1,''omitnan'');'])
    eval(['total_',num2str(bt),'(2,:)=std(temp,0,''omitnan'');'])
    col=find(sum(total_data_1>0)==1);
    eval(['total_',num2str(bt),'(2,col)=NaN;'])
    eval(['total_',num2str(bt),'(3,:)=total_',num2str(bt),'(2,:)./mean(temp,1,''omitnan'')./sqrt(sum(temp1~=0))*100;'])
end

count=0;
for bt=1:b_t-1
    for i=bt+1:b_t
        eval(['test1=total_data_',num2str(bt),';']);
        eval(['test2=total_data_',num2str(i),';']);
        
        count=count+1;
        [h1, p1, ci1, stats1]=ttest2(test1,test2,[],[],'unequal');
%         h(count,:)=h1;
        p(count,:)=p1;
%         ci(count,:)=ci1;
%         stats(count,:)=stats1;        
    end
end

ttt1=ttt(1,:);
ttt1(2,1)=ttt(2,1);
for i=1:size(total_data_1,1)
    temp=cell2mat(ttt(1+(i-1)*6+4,2));
    ttt1(i+1,2)={temp(1:end-5)};
    for j=1:size(total_data_1,2)
        ttt1{i+1,2+j}=total_data_1(i,j);
    end
end

if exist('total_data_2')
    last=i+1;
    ttt1(i+2,1)=ttt(1+(i-1)*6+4+6,1);
    start=1+(i-1)*6+4+5;
    for ii=1:size(total_data_2,1)
        temp=cell2mat(ttt(start+(ii-1)*6+4,2));
        ttt1(ii+last,2)={temp(1:end-5)};
        for j=1:size(total_data_2,2)
            ttt1{ii+last,2+j}=total_data_2(ii,j);
        end
    end
end


a=txtdata(:,2);
total=[];
count3=1;
for i=1:length(a)
    c=strfind(a{i,:},'total');
    if ~isempty(c)
        total=[total i];
        name{numel(total)}=a{i}(1:c-1);
        d_total(length(total),:)=data(total(count3)-1,:);
        count3=count3+1;
    end
end
% tt(size(tt,1)+10,size(tt,2))={''};
tt=txtdata(1:end-b_t,:);
tt(1,size(tt,2)+1)={'Problematic columns'};
for i=2:size(tt,1)
    tt(i,62)={''};
end
% tt(size(tt,1)-1,:)=tt(size(tt,1)-2,:);
% tt(size(tt,1)-2,2)={'total mice (sections)'};
% tt(size(tt,1),2)={'total mice (sections)'};

count2=0;
for i=1:size(tt,1)
    c=strfind(tt{i,1},'t-test : ');
    if ~isempty(c)
        count2=count2+1;
        ind_p(count2)=i;
        eval(['temp',num2str(count2),'=tt{ind_p(',num2str(count2),'),1};'])
        eval(['p_string',num2str(count2),'=temp',num2str(count2),'(10:end);'])        
    end
end
% temp=tt{ind_p,1};
% % p_string=[temp(10),'vs',temp(strfind(temp,'vs')+3)];
% p_string=[temp(10:end)];

for c=1:length(total)-b_t
    name_total{c,:}=[tt{total(c),2},'_trim'];
end

% tt(total(1):end-3,:)=tt(total(1)+3:end,:);
% tt(total(4)-3:end-3,:)=tt(total(4):end,:);
% tt=tt(1:end-6,:);
% 
ori_pval=data(end-3,:);
% tt=tt(1:end-2,:);   % clear p-value

if b_t>1
    ttt(total(1):end-3,:)=ttt(total(1)+3:end,:);
    ttt(total(4)-3:end-3,:)=ttt(total(4):end,:);
    ttt=ttt(1:end-6,:);
    ttt=ttt(1:end-2,:);   % clear p-value
elseif b_t==1
    ttt=ttt(1:end-1,:);   % clear p-value
end


if b_t==1
    total_new=[total(1:3)-3,total(4)-4];
elseif b_t==2
    total_new=[total(1:3)-3,total(4:6)-6,total(7:8)-8];
end
count1=size(tt,1)+1;
% count1=size(tt,1);
for bt=1:b_t
    for i=1:3
        count1=count1+1;
        c=(bt-1)*3+i;
%         tt(size(tt,1)-(11-count1),2)=tt(total(c),2);
        tt(count1,2)=name_total(c,:);
    end
    count1=count1+1;
%     tt(size(tt,1)-(11-count1),2)={[name{end-(2-bt)},'total mice']};
    tt(count1,2)={[name{end-(b_t-bt)},'total mice_trim']};
    count1=count1+1;
%     tt(size(tt,1)-(11-count1),2)={[name{end-(2-bt)},'total sections']};
    tt(count1,2)={[name{end-(b_t-bt)},'total sections_trim']};
end
for i=1:b_t
    eval(['a',num2str(i),'=name{end-',num2str(i-1),'};'])
end

if bt==2
    eval(['ind=find(a1-a',num2str(b_t),');'])
elseif bt==3
    for i=1:min(length(a1),length(a3))
        if a1(i)~=a3(i)
            ind=i;
            break
        end
    end
elseif bt==4
    for i=1:min(length(a1),length(a4))
        if a1(i)~=a4(i)
            ind=i;
            break
        end
    end
end

count1=count1+1;
if b_t>1
    for i=1:factorial(b_t)/(factorial(b_t-2)*factorial(2))
        eval(['tt(count1,2)={[a1(1:ind-1),p_string',num2str(i),','' p_value_trim'']};'])
    %     eval(['tt(ind_p(',num2str(i),'),2)={[a1(1:ind-1),p_string',num2str(i),','' p_value'']};'])
        count1=count1+1;
    end
end

row=1;
for bt=1:b_t
    b_n=bone_n(bt);
    for bn=1:b_n
%         eval(['temp=out_data_',num2str(bt),'_',num2str(bn),';'])
%         for i=1:6
%             row=row+1;
%             for col=1:17
%                 tt(row,col+2)={temp(i,col)};            % new values
%             end
%         end

        eval(['temp=bone_',num2str(bt),'_',num2str(bn),';'])
        for i=1:3
            row=row+1;
            for col=1:vari
                tt(row,col+2)={temp(i,col)};            % old raw values
            end
            eval(['out=''outliers_',num2str(bt),'_',num2str(bn),'_',num2str(i),''';'])
            eval(['out1=outliers_',num2str(bt),'_',num2str(bn),'_',num2str(i),';'])
            if exist(out,'var')
                tem='';
                for k=1:length(out1{1})
                    if ~isempty(tem)
                        tem=[tem,', ',col_name1{out1{1}(k)}];
                    else
                        tem=[tem,col_name1{out1{1}(k)}];
                    end
                end
                tem={tem};
            else
                tem={''};                
            end
            tt(row,col+3)=tem;
        end
        eval(['temp=out_data_',num2str(bt),'_',num2str(bn),';'])
        for i=4:6
            row=row+1;
            for col=1:vari
                tt(row,col+2)={temp(i,col)};            % new values (mean, std, RSE)
            end
        end
    end
    
    eval(['temp=total_',num2str(bt),';'])
    for i=1:3
        row=row+1;
        for col=1:vari
            tt(row,col+2)={temp(i,col)};                % new total values
        end
    end
end

for c=1:count
    row=row+1;
    for col=1:vari
        tt(row,col+2)={p(c,col)};                       % p-value
    end
end

row=row+2;                                              % skip one line
% row=row+1;

tt(row,1)={'trimmed data summary'};
tt(row,2:vari+2)=txtdata(1,2:end);
% row=row+1;


% tt(total(1):end-3,:)=tt(total(1)+3:end,:);
% tt(total(4)-3:end-3,:)=tt(total(4):end,:);
% tt=tt(1:end-6,:);
% 
% tt=tt(1:end-2,:);   % clear p-value

for bt=1:b_t
    for i=1:3
        row=row+1;
        eval(['temp=total_',num2str(bt),';'])
        for col=1:vari
            tt(row,col+2)={temp(i,col)};                    % new total mean at the bottom
        end
    end
    row=row+1;
    for col=1:vari
        eval(['temp=mice_no_',num2str(bt),'(1,col);'])
%         number=[num2str(temp)];
%         tt(row,col+2)={number};                         % new total mice number
        tt(row,col+2)={temp};                         % new total mice number
    end
    row=row+1;
    for col=1:vari
        eval(['temp1=no_',num2str(bt),'(1,col);'])
%         number=[num2str(temp1)];
%         tt(row,col+2)={number};                         % new total section number
        tt(row,col+2)={temp1};                         % new total section number
    end
end

for c=1:count
    row=row+1;
    for col=1:vari
        tt(row,col+2)={p(c,col)};                       % p-value
    end
end

% tt(total(1):end-3,:)=tt(total(1)+3:end,:);
% tt(total(4)-3:end-3,:)=tt(total(4):end,:);
% tt=tt(1:end-6,:);


if b_t>1
    new_summary=tt(end-(5*b_t+factorial(b_t)/(factorial(b_t-2)*factorial(2))):end,:);
else
    new_summary=tt(end-(5*b_t):end,:);
end



count1=size(ttt,1);
if b_t>1
    ttt(size(ttt,1)+11,size(ttt,2))={''};
elseif b_t==1
end
for bt=1:b_t
    eval(['total_mice',num2str(bt),'=sum(section_',num2str(bt),'>0);'])
    eval(['total_section',num2str(bt),'=sum(section_',num2str(bt),');'])
%     total_mice1=sum(section_1>0);
%     total_mice2=sum(section_2>0);
%     total_section1=sum(section_1);
%     total_section2=sum(section_2);
end


% count1=count1+1;
% ttt(count1,1)={'raw data summary'};
% ttt(count1,2:vari+2)=txtdata(1,2:end);
for bt=1:b_t
    eval(['tm=total_mice',num2str(bt),';'])
    eval(['ts=total_section',num2str(bt),';'])
    for i=1:3
        count1=count1+1;
        c=(bt-1)*3+i;

%         ttt(count1,2)=name_total(c,:);
        switch i
            case 1
                ttt(count1,2)={[name{c},'total mean']};
            case 2
                ttt(count1,2)={[name{c},'total std']};
            case 3
                ttt(count1,2)={[name{c},'total RSE(%)']};
        end        
        for k=3:vari+2
%             ttt(count1,k)=ttt(total(c),k);
            ttt{count1,k}=d_total(c,k-2);
        end
    end
    count1=count1+1;
    ttt(count1,2)={[name{end-(b_t-bt)},'total mice']};
    for k=3:vari+2
        ttt(count1,k)={tm(k-2)};
    end
    count1=count1+1;
    ttt(count1,2)={[name{end-(b_t-bt)},'total section']};
    for k=3:vari+2
        ttt(count1,k)={ts(k-2)};
    end
end

% ttt(total(1):end-3,:)=ttt(total(1)+3:end,:);
% ttt(total(4)-3:end-3,:)=ttt(total(4):end,:);
% ttt=ttt(1:end-6,:);

if bt==2
    a1=name{end};
    a2=name{end-1};
    ind=find(a1-a2);
    % count1=count1+1;
    % ttt(count1,2)={[a1(1:ind-1),p_string,' p_value']};
    % ttt(ind_p,2)={[a1(1:ind-1),p_string,' p_value']};
    % % count1=count1+1;
    % % ttt(count1,2)={'p_value'};
elseif bt==3
    for i=1:min(length(a1),length(a3))
        if a1(i)~=a3(i)
            ind=i;
            break
        end
    end
end

count1=count1+1;
if b_t>1
    for i=1:factorial(b_t)/(factorial(b_t-2)*factorial(2))
        eval(['ttt(count1,2)={[a1(1:ind-1),p_string',num2str(i),','' p_value'']};'])
    %     eval(['ttt(ind_p(',num2str(i),'),2)={[a1(1:ind-1),p_string',num2str(i),','' p_value'']};'])
        for k=3:vari+2
    %         ttt(count1,k)=ttt(p_pos,k);
            ttt{count1,k}=ori_pval(1,k-2);
        end
        count1=count1+1;
    end
end

% for i=1:length(ind_p)
%     for k=3:27
% %         ttt(count1,k)=ttt(p_pos,k);
%         ttt(count1,k)=ttt(ind_p(i),k);
%     end
% end

count1=size(ttt,1)+1;
for i=1:size(new_summary,1)%12
    count1=count1+1;
    for j=1:vari+2
        ttt(count1,j)=new_summary(i,j);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Added for raw_mean beneath trimmed_mean
%   February 03, 2015
%
ChangeRow_tt=[];
for i=1:size(tt,1)
    phr=tt{i,2};
    if strfind(phr, '_mean')
        ChangeRow_tt=[ChangeRow_tt, i];
    end
end

Row_ttt=[];
for i=1:size(ttt,1)
    phr=ttt{i,2};
    if strfind(phr, '_mean')
        Row_ttt=[Row_ttt, i];
    end
end

for i=length(ChangeRow_tt):-1:1
    r=ChangeRow_tt(i);
    phr=tt{r,2};
    ind=strfind(phr, '_mean');
    tt{r,2}=[phr(1:ind-1),'_mean_trim'];
    tt{r+1,2}=[phr(1:ind-1),'_std_trim'];
    tt{r+2,2}=[phr(1:ind-1),'_RSE(%)_trim'];
    [tt{r+2:end+1,:}]=tt{r+1:end,:};
    [tt{r+1,1:end-1}]=ttt{Row_ttt(i),:};
end
%
%   Added for raw_mean beneath trimmed_mean
%   February 03, 2015
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Added for raw_mean beneath trimmed_mean
%   February 11, 2015
%
for i=1:size(ttt,1)
    phr=ttt{i,2};
    if strfind(phr, 'Bone_no')
       new_row=i;
       break
    end
end
tt{end+2,1}='raw data summary';
[tt{end,2:end-1}]=ttt{new_row,2:end};
for i=1:size(ttt,1)
    phr=ttt{i,2};
    if strfind(phr, 'total mean')
       st_row=i;
       break
    end
end
for i=1:size(ttt,1)
    phr=ttt{i,2};
    if strfind(phr, 'p_value')
       ed_row=i;
       break
    end
end
if exist('ed_row','var')==0
    ed_row=size(ttt,1);
end    
len=ed_row-st_row+1;
[tt{end+1:end+1+len-1,1:end-1}]=ttt{st_row:ed_row,:};
%
%   Added for raw_mean beneath trimmed_mean
%   February 11, 2015
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Added for delete total mean,std RSE female and total mean,std RSE male in trimmed tab
%   March 05, 2015
%
count=0;
for i=1:size(tt,1)
    phr=tt{i,2};
    if strfind(phr, '_total mean') & strfind(tt{i+1,2}, '_total std') & strfind(tt{i+2,2}, '_total RSE(%)')
        count=count+1;
        if count<b_t
            [tt{i:end-3,:}]=tt{i+3:end,:};
            tt=tt(1:end-3,:);
        elseif count==b_t
            [tt{i:end-6,:}]=tt{i+6:end,:};
            tt=tt(1:end-6,:);
            break
        end
    end
end
tt{i,1}='trimmed data summary';
%
%   Added for delete total mean,std RSE female and total mean,std RSE male in trimmed tab
%   March 05, 2015
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Added for 'raw data summary' at the 1st comumn of raw data summary
%   March 05, 2015
%
for i=1:size(ttt,1)
    phr=ttt{i,2};
    if strfind(phr, '_total mean')
        ttt{i,1}='raw data summary';
        break
    end
end
%
%   Added for 'raw data summary' at the 1st comumn of raw data summary
%   March 05, 2015
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

col_name1(18)={'Tb_Th'};
col_name1(19)={'Tb_N'};
col_name1(20)={'Tb_SP'};

% if use_mean_median==1
%     eval(['xlswrite(''',output_excel_file,num2str(margin),'SD_mean.xls'',ttt,1)'])
%     eval(['xlswrite(''',output_excel_file,''',ttt,''raw data'')'])
% else
%     eval(['xlswrite(''',output_excel_file,num2str(margin),'SD_median.xls'',ttt,1)'])
% end

% Bone_type=ttt(2:end-13,1);
Bone_no=ttt(2:end-13,2);
phrase=['T=table(Bone_no,'];
for i=3:size(ttt,2)
    eval(['name=cell2mat(col_name1(',num2str(i-2),',1));'])
    eval([name,'=cell2mat(ttt(2:end-13,',num2str(i),'));'])
%     eval([name,'(end+2:end+13)=cell2mat(ttt(end-11:end,',num2str(i),'));'])
    phrase=[phrase,name,','];
end
phrase=[phrase(1:end-1),');'];
eval(phrase);
eval(['writetable(T,''',output_excel_file,''',''sheet''',',''raw data'',''range'',''B1:BI108'');'])

% Bone_type=ttt(end-10:end,1);
Bone_no=ttt(end-10:end,2);
phrase=['T=table(Bone_no,'];
for i=3:size(ttt,2)
    eval(['name=cell2mat(col_name1(',num2str(i-2),',1));'])
    eval([name,'=cell2mat(ttt(end-10:end,',num2str(i),'));'])
%     eval([name,'(end+2:end+13)=cell2mat(ttt(end-11:end,',num2str(i),'));'])
    phrase=[phrase,name,','];
end
phrase=[phrase(1:end-1),');'];
eval(phrase);
eval(['writetable(T,''',output_excel_file,''',''sheet''',',''raw data'',''range'',''B110:BI121'');'])




% if use_mean_median==1
%     eval(['xlswrite(''',output_excel_file,num2str(margin),'SD_mean.xls'',tt,2)'])
%     eval(['xlswrite(''',output_excel_file,''',tt,''trimmed data'')'])
%     eval(['writematrix(''',output_excel_file,''',tt,''trimmed data'')'])

Bone_no=tt(2:end-13,2);
phrase=['T=table(Bone_no,'];
for i=3:size(tt,2)-1
    eval(['name=cell2mat(col_name1(',num2str(i-2),',1));'])
    eval([name,'=cell2mat(tt(2:end-13,',num2str(i),'));'])
%     eval([name,'(end+2:end+13)=cell2mat(ttt(end-11:end,',num2str(i),'));'])
    phrase=[phrase,name,','];
end

phrase=[phrase(1:end-1),');'];
eval(phrase);
eval(['writetable(T,''',output_excel_file,''',''sheet''',',''trimmed data'',''range'',''B1:BI124'');'])

% Bone_type=ttt(end-10:end,1);
Bone_no=tt(end-10:end,2);
phrase=['T=table(Bone_no,'];
for i=3:size(tt,2)-1
    eval(['name=cell2mat(col_name1(',num2str(i-2),',1));'])
    eval([name,'=cell2mat(tt(end-10:end,',num2str(i),'));'])
%     eval([name,'(end+2:end+13)=cell2mat(ttt(end-11:end,',num2str(i),'));'])
    phrase=[phrase,name,','];
end
phrase=[phrase(1:end-1),');'];
eval(phrase);
eval(['writetable(T,''',output_excel_file,''',''sheet''',',''trimmed data'',''range'',''B126:BI137'');'])

Problematic_columns=tt(2:109,end);
T=table(Problematic_columns);
writetable(T,output_excel_file,'sheet','trimmed data','range','BJ1:BJ109');


% else
%     eval(['xlswrite(''',output_excel_file,num2str(margin),'SD_median.xls'',tt,2)'])
% end
%     eval(['xlswrite(''',output_excel_file,''',ttt1,''trimmed average data'')'])
%     eval(['writematrix(''',output_excel_file,''',ttt1,''trimmed average data'')'])

Bone_no=ttt1(2:end,2);
phrase=['T=table(Bone_no,'];
for i=3:size(ttt1,2)
    eval(['name=cell2mat(col_name1(',num2str(i-2),',1));'])
    eval([name,'=cell2mat(ttt1(2:end,',num2str(i),'));'])
%     eval([name,'(end+2:end+13)=cell2mat(ttt(end-11:end,',num2str(i),'));'])
    phrase=[phrase,name,','];
end

phrase=[phrase(1:end-1),');'];
eval(phrase);
eval(['writetable(T,''',output_excel_file,''',''sheet''',',''trimmed average data'',''range'',''B1:BI17'');'])
