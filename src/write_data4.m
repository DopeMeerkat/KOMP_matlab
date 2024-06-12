function write_data4(varargin)

% function write_data4(direct, phrase1, phrase2, bone_type, delimeter)

% direct='Z:\KOMP\Bzw_E01V\01_Submitted\Layers\Images\data\';
% phrase1={'Bzw2_E01_hV_F';'Bzw2_E01_hV_M'};
% phrase2={'Female';                              % Exp name for output excel file
%       'Male'};
% bone_type='V';                                % Femur, Vertebra
% delimeter='\';

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
    ana_n=1;               % ana_n=0 : first analysis, ana_n=1 : second analysis
% end


no_roi=1;

tic

if strcmp(direct(end),'\')~=1
    direct=[direct,'\'];
end

fclose('all')

aaa=phrase1{1,:};

eval(['out_file1=''',direct,aaa(1:end-1),'analysis1_tr1.txt'';'])
% eval(['fid_g1 = fopen(''',out_file1,''', ''w'');'])
fid_g1 = fopen(out_file1,'w');
fprintf(fid_g1,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n', ...
        'Bone type','Bone #','R_BS', 'G_BS', 'R_G', 'Ronly_BS', 'Gonly_BS', 'sLS_BS', 'dLS_BS', 'LS_BS', 'R_LS', 'G_LS', 'MS_BS', 'sLS_LS', 'dLS_LS', 'dLS_sLS', 'MAR_um', 'BFR', 'BV_TV', 'Tb.Th', 'Tb.N', 'Tb.Sp');
if no_roi==2
    eval(['out_file2=''',direct,phrase1{1,:},'analysis1_tr2.txt'';'])
    fid_g2 = fopen(out_file2, 'w');
    fprintf(fid_g2,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n', ...
        'Bone type','Bone #','R_BS', 'G_BS', 'R_G', 'Ronly_BS', 'Gonly_BS', 'sLS_BS', 'dLS_BS', 'LS_BS', 'R_LS', 'G_LS', 'MS_BS', 'sLS_LS', 'dLS_LS', 'dLS_sLS', 'MAR_um', 'BFR', 'BV_TV', 'Tb.Th', 'Tb.N', 'Tb.Sp');
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
end

NaN_flag=0;             % if Nan_flag==0 --> there is no NaN, if 1 --> there is NaN

new_bone=0;



for loop=1:no_roi
    for b_t=1:size(phrase1,1)
        b_t
        eval(['bone_number=bone_number',num2str(b_t),';']);
        for b_n=1:size(bone_number,1)

            bn=str2num(bone_number{b_n,1})


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
                if isempty(aa)==1
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
                phras=[phras,'=nanmean(data,1);'];
                eval(phras)
                
                if size(data,1)>1
                    phras=['std_tr',num2str(loop),'_',num2str(b_t),'(',num2str(b_n),',:)'];
                    phras=[phras,'=nanstd(data);'];
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




for loop=1:no_roi
    for b_t=1:size(phrase1,1)       
        eval(['bone_number=bone_number',num2str(b_t),';']);
        count=0;
        minus_count(b_t,length(bone_number))=0;
        
        for b_n=1:size(bone_number,1)
            if b_n==3
                b_n;
            end
            
%             nnn=num2str(bone_number(b_n));
            nnn=bone_number{b_n,1};

            if ana_n==0
                phr1_1=[phrase1{b_t,:},num2str(bone_number(b_n)),'_1_analysis_tr',num2str(loop),'.mat'];
                phr1_2=[phrase1{b_t,:},num2str(bone_number(b_n)),'_2_analysis_tr',num2str(loop),'.mat'];
                phr1_3=[phrase1{b_t,:},num2str(bone_number(b_n)),'_3_analysis_tr',num2str(loop),'.mat'];
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
                if (s==1 & b_n==1) %| samples_per_bone(b_t, b_n)==1
                    fprintf(fid_g1,'%s\t%s\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\n',...
                        phrase2{b_t,:}, section_name, d);
                else
                    fprintf(fid_g1,'\t%s\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\n',...
                        section_name, d);
                end
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
            eval(['fid_g=fid_g',num2str(loop),';'])

            fprintf(fid_g,'\t%s\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\n',...
                section_name_mean, d_mean);
            fprintf(fid_g,'\t%s\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\n',...
                section_name_std, d_std);
            fprintf(fid_g,'\t%s\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\n',...
                section_name_RSE, d_RSE);

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
        eval(['fid_g=fid_g',num2str(loop),';'])
        eval(['len=(size(mean_tr',num2str(loop),'_',num2str(b_t),',1)-minus_count(',num2str(b_n),'));'])
%         eval(['me=nansum(mean_tr',num2str(loop),'_',num2str(b_t),',1)/len;'])
%         eval(['stdev=sqrt(nansum((mean_tr',num2str(loop),'_',num2str(b_t),'-ones(size(mean_tr',num2str(loop),'_',num2str(b_t),',1),1)*me).^2)/(len-1));'])
        eval(['me=nanmean(mean_tr',num2str(loop),'_',num2str(b_t),',1);'])
        eval(['stdev=nanstd(mean_tr',num2str(loop),'_',num2str(b_t),',0);'])
        fprintf(fid_g1,'\t%s\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\n',...
            section_name_mean, me);
        fprintf(fid_g1,'\t%s\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\n',...
            section_name_std, stdev);
        fprintf(fid_g1,'\t%s\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\n',...
            section_name_RSE, stdev./me/sqrt(len)*100);

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
            eval(['fid_g=fid_g',num2str(loop),';']);
            for loop=1:no_roi
                fprintf(fid_g,'%s\t%s\t%6.4e\t%6.4e\t%6.4e\t%6.4e\t%6.4e\t%6.4e\t%6.4e\t%6.4e\t%6.4e\t%6.4e\t%6.4e\t%6.4e\t%6.4e\t%6.4e\t%6.4e\t%6.4e\t%6.4e\t%6.4e\t%6.4e\t%6.4e\n',...
                    p1,'p_value', pp);
            end

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
        eval(['fid_g=fid_g',num2str(loop),';'])
        eval(['me=nanmean(mean_tr',num2str(loop),'_',num2str(b_t),',1);'])
        fprintf(fid_g1,'\t%s\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\n',...
            section_name_mean, me);
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
    output_image_file=[direct,'RSE_Histo_',cell2mat(col_name2(loop))];
    eval(['print -djpeg100 ''',output_image_file,''''])
    close all
end
