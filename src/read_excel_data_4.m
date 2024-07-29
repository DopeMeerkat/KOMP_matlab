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
major=[1, 2, 15, 17, 18, 21, 35, 42, 45, 46, 47, 48, 53];         % R/BS, G/BS, MAR, BV/TV, Tb.Th, AP/BS, TRAP/BS, GFP/BS, Height/Width, Osteocytes_n, Cortex_width, AC_Intensity, Calcein_Intensity 

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
            max_range(i)=nanmean(tttt(tttt>0))+nanstd(tttt(tttt>0))*margin;
            min_range(i)=nanmean(tttt(tttt>0))-nanstd(tttt(tttt>0))*margin;
        end
    elseif use_mean_median==0
        for i=1:size(t_data,2)
            tttt=t_data(:,i);
            max_range(i)=nanmedian(tttt(tttt>0))+nanstd(tttt(tttt>0))*margin;
            min_range(i)=nanmedian(tttt(tttt>0))-nanstd(tttt(tttt>0))*margin;
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
        temp(4,:)=nanmean(temp(1:3,:));
        mice_number=~isnan(temp(4,:));
        eval(['mice_no_',num2str(bt),'=mice_no_',num2str(bt),'+mice_number;'])        
        tt=nanstd(temp(1:3,:)); tt(tt==0)=NaN;
        temp(5,:)=tt;
        temp(6,:)=temp(5,:)./temp(4,:)./sqrt(sum(temp1~=0))*100;
        eval(['out_data_',num2str(bt),'_',num2str(bn),'=temp;'])
        eval(['total_data_',num2str(bt),'(',num2str(bn),',:)=temp(4,:);'])
    end
end

for bt=1:b_t
    eval(['temp=total_data_',num2str(bt),';'])
    temp1=temp; temp1(isnan(temp))=0;
    eval(['total_',num2str(bt),'(1,:)=nanmean(temp,1);'])
    eval(['total_',num2str(bt),'(2,:)=nanstd(temp,0);'])
    col=find(sum(total_data_1>0)==1);
    eval(['total_',num2str(bt),'(2,col)=NaN;'])
    eval(['total_',num2str(bt),'(3,:)=total_',num2str(bt),'(2,:)./nanmean(temp,1)./sqrt(sum(temp1~=0))*100;'])
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
        number=[num2str(temp)];
        tt(row,col+2)={number};                         % new total mice number
    end
    row=row+1;
    for col=1:vari
        eval(['temp1=no_',num2str(bt),'(1,col);'])
        number=[num2str(temp1)];
        tt(row,col+2)={number};                         % new total section number
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
    if strfind(phr, 'Bone #')
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


% if use_mean_median==1
%     eval(['xlswrite(''',output_excel_file,num2str(margin),'SD_mean.xls'',ttt,1)'])
    eval(['xlswrite(''',output_excel_file,''',ttt,''raw data'')'])
% else
%     eval(['xlswrite(''',output_excel_file,num2str(margin),'SD_median.xls'',ttt,1)'])
% end

% if use_mean_median==1
%     eval(['xlswrite(''',output_excel_file,num2str(margin),'SD_mean.xls'',tt,2)'])
    eval(['xlswrite(''',output_excel_file,''',tt,''trimmed data'')'])
% else
%     eval(['xlswrite(''',output_excel_file,num2str(margin),'SD_median.xls'',tt,2)'])
% end
    eval(['xlswrite(''',output_excel_file,''',ttt1,''trimmed average data'')'])
