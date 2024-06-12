function analize_KOMP_for_BioHead_6(root_di, direct, out_dir, phr, delimeter, bt, bn, section, GFP_flag, RG_label, mouse, rat_mouse)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   ver 2.7 : using move_3 : (output : minimum distance from the label, GFP)
%   ver 2.8 : using move_4 : (decision of green label and GFP is done by comparison of distances
%                             distance between the red label and the boundary, and
%                             distance between the green signal and the boundary)
%                             if dr>ds --> (signal = GFP)
%                             if dr<ds --> (signal = green label)
%   ver 2.9 : using find_min : find the minimum of two green(red) trimmed labels (1/7/08)
%                              if there are more than 2 labels are overlapped on the same trimmed region
%   ver 2.10 : separate single into red only and green only (1/11/08)
%   ver 2.11 : add dLS/sLS (1/17/08)
%               substitute find with logical indexing such as
%               'stats(II(find(tt==s_tt(1))))'-->'stats(II(logical([1 0 0 0])))'
%   ver 2.12 : using move_5 : save alternate trim candidate (1/23/08)
%              save alternate trim candidate
%              calculate MAR(Mineral Apposition Rate), then recalculate MAR within sigma*n
%               if MAR of a certain label is greater than sigma*n, then treat that green label as a single label
%   ver 2.13 : using move_6 : recalculate with MAR with trimmed red_label and trimmed green_label(1/23/08)
%   ver 2.14 : if distance between red and green is negative,
%               it is GFP and put those green area with average distance of green in that specific TB(2/4/08)
%              using move1_2(2/4/08)
%   ver 2.15 : 10x magnification with resizing the image to 1/2
%   ver 2.16 : Find the broken cortical bones and put them into one cortical bone (3/5/08)
%   ver 2.17 : Find left, right cortical bones and put them into one cortical bone (3/6/08)
%
%   ver 9    : find the mid-points of labels instead of leading edges (12/01/10)
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    



% root_di='C:\Users\shhong\Desktop\Histomorphometry\GS03   (sent 07-22-13)\01_Submitted\Layers';
% direct='C:\Users\shhong\Desktop\Histomorphometry\GS03   (sent 07-22-13)\01_Submitted\Layers\Images\';
% out_dir='C:\Users\shhong\Desktop\Histomorphometry\GS03   (sent 07-22-13)\02_Analyzed\';          % directory to write the analized image
% phr={'GS03_CF';
%       'GS03_CM';
%       'GS03_TF';
%       'GS03_TM'};
% auto_exp=0;
% mouse = 1; % --> mouse(default)
            % 0; % --> rat
% rat_mouse=5028/1913; % ratio of size rat over mouse

if strcmp(root_di(end),delimeter)==0
    root_di=[root_di,delimeter];
end
thresh_dir=[root_di,'Threshold',delimeter];

if nargin==10
    mouse=1;
    rat_mouse=1;
elseif nargin==11
    if mouse==1
        rat_mouse=1;
    elseif mouse==0
        rat_mouse=5028/1913;
    end
end

phr1=phr;
% for i=1:length(phr1)
%     phr1(i,:)={[phr1{i,:},'_']};
% end

  
  
  
% segment tetraycline with 2*std_dev 
close all
%clear
matlab_ver=7.1;
warning('off')

exp_type='F';
%exp_type='';
trap_flag=0;
ap_flag=0;
bubble_threshold=0.5;                       % Bubble intensity threshold
segment_threshold=0.0002;                   % Threshold for segmenting DIC and labels
start_threshold_ratio=-60;                  % start_search of multiplication number for segmenting DIC and labels (start_threshold_ratio/10)

dis= [10, 1;
      10, 0;
      6, 0];
%   [green_horizontal, green_vertical;
%   [red_horizontal,   red_vertical;
%   [yellow_horizontal, yellow_vertical;]

%
%dis_x=6;
%dis_y=1;

if ~strcmp(root_di(end),delimeter)
    root_di=[root_di,delimeter];
end

if ~strcmp(direct(end),delimeter)
    direct=[direct,delimeter];
end

direc_ana=[direct,'data',delimeter];
if isdir(direc_ana)==0
    mkdir(direc_ana)
end


folder_name=[phr(1,:)];
%out_direct1='C:\Users\shhong\Desktop\jpeg_HisBL6_8wk_V';
out_direct1=direct;
aa=dir(out_direct1);
d=[out_direct1,cell2mat(folder_name)];
uchc_d=[out_dir];%,phr];
% if isdir(d)==0
%    mkdir(d)
% end

med_bubble=0;                                           % 0 : shift to match beads
                                                        % 1 : shift + median filter
                                                        % 2 : shift + median filter + remove bubbles
if med_bubble==0
    mb_com='_shift3.jpg';
elseif med_bubble==1
    mb_com='_med.jpg';
elseif med_bubble==2
    mb_com='_med_bubble.jpg';
end
                                                        
red_first_green_last=0;                                 % if 1 : 1st injection --> red
                                                        %        2nd injection --> green
                                                        % if 0 : 1st injection --> green

% for i=1:size(phr,1)
%     eval(['load ''',root_di,'threshold_values',num2str(i),'''']);
% end

% eval(['load ''',root_di,'threshold_values2''']);
% eval(['load ''',root_di,'threshold_values2''']);

eval(['load ''',root_di,'info''']);
% end_growth_plate=growth_plate1(:,:,1);
% end_growth_plate(:,:,2)=growth_plate2(:,:,1);
% end_points=growth_plate1(:,:,2);
% end_points(:,:,2)=growth_plate2(:,:,2);

angle=0;%-90;
%[section bone]=find(threshold(:,:,1)');                                                      %        2nd injection --> red


                                                        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       Read individual images and Save combined image
%
%read_write_images(direct, dis, color_order, bubble_threshold, phr1, samples_per_bone, bone_number1);
%read_write_images(direc, dis, color_order, bubble_threshold, phr, sample_per_bone, bone_number1, bone_number2, bone_number3, bone_number4);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


tic
%std_no_green=2;
%std_no_red=1.5;
down_margin=50;
down_margin=150;
std_no_green=0.5;
std_no_red=-0.2;
no_analysis=13;
area_roi_tb=800;

% if mouse==1
%     ratio=2;
    ratio=1;
% else
%     ratio=1;
% end


std_const=2;
days2_1=4;                         % days between injection 2 and injection 1 (for labels)
%dis_5x=312;                    % 312 pixels in 5x --> 400 um
% dis_5x=round(312*rat_mouse);                    % 312 pixels in 5x --> 400 um
% dis_5x=round(624*rat_mouse);                    % 312 pixels in 5x --> 400 um
dis_5x=round((144*4)*rat_mouse);                    % 576 pixels in 5x --> 400 um
% dis_10x=round(dis_5x*2/ratio);              % 312*2 pixels in 10x --> 400 um (based on the image 'C:\Users\shhong\Desktop\CJake\scale.jpg')

dis_10x=round(dis_5x*2/ratio);              % 312*2 pixels in 10x --> 400 um (based on the image 'C:\Users\shhong\Desktop\CJake\scale.jpg')

dis_10x_for_400micron_per_312pixels=400/dis_10x;    % distance for 10x image as ratio of 400 um per 312 pixels
%dis_5x=dis_5x*2;
margin_5x=200*rat_mouse;                  % 230 pixel
%margin_5x=margin_5x*2;                  % 230 pixel

start_trabecula=dis_10x;            % 400 um
%start_trabecula=150;            % 150 um
%start_trabecula=76;            % 150 um
sp_tb=dis_10x;                  %*start_trabecula/400;      % start point of trabecula from the growth plate (400 um)
%lr_sp_tb=round(dis_10x*start_trabecula/4/400);      % start point of trabecula from the endosteum (250 um)
lr_sp_tb=round(dis_10x/2);      % start point of trabecula from the endosteum (200 um)
margi=margin_5x;
margi_lr=156;
%margi_lr=margi_lr*2;
%start_point=550/ratio;
%end_point=4700/ratio;


% start_point=round(1/ratio);
% 
% %start_point=start_point+300+dis_10x;
% %start_point=start_point+1300+dis_10x;
% 
% end_point=1500/ratio;

% red_flag=1;
% GFP_flag=0;
% TRAP_flag=1;
% AP_flag=1;
exact_boundary=1;
loop_count=20;
%rot_angle=-15;
rot_angle=0;%-90;
thresh=0.0001;
%thresh=0.00005;
write=1;                % 0 : test without writing
                        % 1 : normal analizing with writing data and images
with_rg_label=1;        % 1 : using red label and green label to segment DIC
                        % 0 : not using red label and green label to segment DIC

save parameters


max_dis=max(dis);
%if isempty(bone_number3)==1
%    if isempty(bone_number2)==1
%        no_bone_type=1;
%    else
%        no_bone_type=2;
%    end
%else
%    no_bone_type=3;
%end
no_bone_type=size(phr1,1);

             
for b_t=bt      % 1:no_bone_type
    
%     eval(['threshold=threshold_values',num2str(b_t),';'])
%     [x y]=find(threshold(:,:,1));
%     bone_number=unique(x);
    eval(['ifm=info',num2str(b_t),';'])
    eval(['inform=ifm.im;'])
    fn=fieldnames(ifm);
    for i=1:length(fn)
        if strcmp(fn{i},'ap')==0
            tap_flag=0;
            continue
        else
            eval(['inform_ap=info',num2str(b_t),'.ap;'])
            tap_flag=1;
            break
        end
    end
    for i=1:length(fn)
        if strcmp(fn{i},'tr')==0
            ttr_flag=0;
            continue
        else
            eval(['inform_tr=info',num2str(b_t),'.tr;'])
            ttr_flag=1;
            break
        end
    end
    
%     bone_number=find(sum(inform,2));
    temp=inform;
%     temp1=sum(cell2mat(temp(:,2:end)),2);
    bone_number=temp(:,1);
    
    eval(['bone_number',num2str(b_t),'=bone_number;']);
    eval(['load ''',root_di,'growth_plate',num2str(b_t),'''']);
    eval(['load ''',root_di,'rot_ang''']);

%     eval(['load ''',root_di,'growth_plate1''']);
%     eval(['load ''',root_di,'growth_plate2''']);
    
%     eval(['end_growth_plate=growth_plate',num2str(b_t),'(:,:,1);'])
%     eval(['end_points=growth_plate',num2str(b_t),'(:,:,2);'])
    eval(['growth=growth_plate',num2str(b_t),';'])
    eval(['rot_an=rot_ang',num2str(b_t),';'])
    
    sample=size(growth,3);
    for sss=1:sample
%         end_growth_plate(:,sss)=[growth{:,2,sss}]';
%         end_points(:,sss)=[growth{:,3,sss}]';
        for i=1:size(growth,1)
            temp=growth{i,2,sss};
            temp1=growth{i,3,sss};
            if isempty(temp)
                end_growth_plate(i,sss)=0;
                end_points(i,sss)=0;
            else
                end_growth_plate(i,sss)=temp;
                end_points(i,sss)=temp1;
            end
        end
    end
    %eval(['bone_number=bone_number',num2str(b_t),';']);
    for b_n=bn      % 1:length(bone_number)
%         samples_per_bone=length(find(threshold(bone_number(b_n),:,1)));
        if size(end_growth_plate,1)<b_n
            continue
        end
        sample=zeros(1,length(end_growth_plate(b_n,:)));
        x=find(end_growth_plate(b_n,:));
        sample(x)=x;
        samples_per_bone=length(sample);
    
        
        for ss=section      % 3:samples_per_bone
            tic
%             archive_d=[uchc_d,phr1{b_t,:},bone_number{b_n}];
            archive_d=[uchc_d];
%             archive_d=[uchc_d,phr1{b_t,:},'0',bone_number{b_n}];
            if isdir(archive_d)==0
                mkdir(archive_d)
            end
            
            cc=0;
            for lo=1:size(rot_an,1)
                if ~isempty(strfind(rot_an{lo,1},bone_number{b_n}))
                    rot_a=rot_an{lo,section+1};
                    if isempty(rot_a)
                        rot_a=0;
                    else
                        cc=1;
                    end
                else
                    continue
                end
            end
            if cc==0 | rot_a==0
                ang=0;
            else
                if rot_a>0;
                    ang=(90-rot_a);
                elseif rot_a<0
                    ang=-(90+rot_a);
                end
            end
            
            if isempty(find(sample==ss))
                continue
            end
            
            s=sample(ss);
            ['Start ', phr1{b_t,:},'L',num2str(s),'_s',bone_number{b_n}, ' --> ', datestr(clock)]
%             ['Start ', phr1{b_t,:},'0',bone_number{b_n}, '_',num2str(s), ' --> ', datestr(clock)]
            
            
           %if (b_n==9 | b_n==14) & s==3
            %    break
            %end

%             load parameters
%             if size(threshold,3)>3
%                 if threshold(bone_number(b_n),s,4)~=0
%                     trap_flag=1;
%                 else
%                 trap_flag=0;
%                 end
%             else
%                 trap_flag=0;
%             end
% 
%             if size(threshold,3)>4
%                 if threshold(bone_number(b_n),s,6)~=0
%                     ap_flag=1;
%                 else
%                     ap_flag=0;
%                 end
%             else
%                 ap_flag=0;
%             end
            
            if tap_flag==1
                if inform_ap{b_n,s+1}==1
                    ap_flag=1;
                else
                    ap_flag=0;
                end
            end
            if ttr_flag==1
                if inform_tr{b_n,s+1}==1
                    trap_flag=1;
                else
                    trap_flag=0;
                end
            end
            
            
            %if samples_per_bone==1
            %    phr=[phr1(b_t,:), bone_number{b_n}, 'FL.jpg_Files\'];
            %else
%                 phr=[phr1(b_t,:), bone_number{b_n}, 'FL_',num2str(s),'.jpg_Files\'];
%                 phr_trap=[phr1(b_t,:), bone_number{b_n}, 'FL_',num2str(s),'(trap).jpg_Files\'];
%                 phr_ap=[phr1(b_t,:), bone_number{b_n}, 'FL_',num2str(s),'(DAPI-AP).jpg_Files\'];
            %end
%             direc=[direct, phr];
%             direc_trap=[direct, phr_trap];
%             direc_ap=[direct, phr_ap];
            direc=direct;
            nnn=bone_number{b_n};
%             nnn=['0',bone_number{b_n}];
            %if s==1 | s==2 | s==3

            c=phr1{b_t,:};
%            c1=[c,nnn,'_c'];
            %if samples_per_bone==1
            %    c1=[c,nnn,'FL_c'];
            %else
%                c1=[c,nnn,'FL_',num2str(s),'_c'];
                if strcmp(exp_type,'V')
                    c1=[c,nnn,'_h',exp_type,'_s',num2str(s)];
                elseif strcmp(exp_type,'F')
                    c1=[c,'L',num2str(s),'_s',nnn];
                end
%                 c1_trap=[c,nnn,'FL_',num2str(s),'(trap)_c'];
%                 c1_ap=[c,nnn,'FL_',num2str(s),'(DAPI-AP)_c'];
            %end
%            c1=[c,'\',c,'_',nnn,'FL\',c,'_',nnn,'FL.jpg_Files\',c,'_',nnn,'FL_c'];
            if ap_flag==1
                file_name={[thresh_dir,c1,'_2',mb_com];
                        [thresh_dir,c1,'_1',mb_com];
                        [thresh_dir,c1,'_0',mb_com];
                        [thresh_dir,c1,'_5',mb_com]};
%                 thre=[threshold(bone_number(b_n),s,3),threshold(bone_number(b_n),s,2),threshold(bone_number(b_n),s,1),threshold(bone_number(b_n),s,6)];
            else
                file_name={[thresh_dir,c1,'_2',mb_com];
                        [thresh_dir,c1,'_1',mb_com];
                        [thresh_dir,c1,'_0',mb_com]};
%                 thre=[threshold(bone_number(b_n),s,3),threshold(bone_number(b_n),s,2),threshold(bone_number(b_n),s,1)];
            end
            
%             [r, g, AP]=clear_beads(file_name, thre); %, mouse);
%              eval(['g=imread(''',direc,c1,'_0',mb_com,''',''jpg'');']);
%              eval(['r=imread(''',direc,c1,'_1',mb_com,''',''jpg'');']);

             eval(['a=dir(''',thresh_dir,c1,'_2',mb_com,''');']);
             if isempty(a) 
                 ['No file : ',thresh_dir,c1,'_2',mb_com]
                 return
             end
             
            'ImageToolbox License test start'
            [l m]=license('checkout','image_Toolbox');
            while (~l)
               pause(1);
            [l m]=license('checkout','image_Toolbox');
            end
            'ImageToolbox License test end'
            
             eval(['a=imread(''',thresh_dir,c1,'_2',mb_com,''',''jpg'');']);
             a=a(:,:,1);
             a=a>127; ori_size=size(a,1);
             a=imrotate(a,ang);
%              a=(a>graythresh(a)*max(a(:)));
             eval(['r=imread(''',thresh_dir,c1,'_1',mb_com,''',''jpg'');']);
             r=r>127;
             r=imrotate(r,ang);
%              r=(r>graythresh(r)*max(r(:)));
             eval(['g=imread(''',thresh_dir,c1,'_0',mb_com,''',''jpg'');']);
             g=g>127;
             g=imrotate(g,ang);
%              g=(g>graythresh(g)*max(g(:)));
%%%%%%%%%%%%%%%%%%%%%
%
%	June, 7th 2018
%
			if ~RG_label % 1st label : red,  2nd label : green --> Reversed labels
				temp=g;
				g=r;
				r=temp;
			end
%
%	June, 7th 2018
%
%%%%%%%%%%%%%%%%%%%%%
			
            r=xor(r,g) & (r & ~g);  % logical(r-g) 

            if ap_flag==1
                check_file=dir([thresh_dir,c1,'_5',mb_com]);
                if ~isempty(check_file)
                    eval(['AP=imread(''',thresh_dir,c1,'_5',mb_com,''',''jpg'');']);
                    if size(AP,3)==3
                        AP=AP(:,:,1);
                    end
                    AP=(AP>127);
                    AP=imrotate(AP,ang);
                else
                    AP=false(size(a));
                end

                check_file=dir([thresh_dir,c1,'_6',mb_com]);                              %   March 24, 2016
                if ~isempty(check_file)                                                   %   March 24, 2016
                    eval(['DAPI=imread(''',thresh_dir,c1,'_6',mb_com,''',''jpg'');']);    %   March 24, 2016
                    if size(DAPI,3)==3                                                    %   March 24, 2016
                        DAPI=DAPI(:,:,3);                                                 %   March 24, 2016
                    end                                                                   %   March 24, 2016
                    DAPI=(DAPI>127);                                                      %   March 24, 2016
                    DAPI=imrotate(DAPI,ang);                                              %   March 24, 2016
                else                                                                      %   March 24, 2016
                    DAPI=false(size(a));                                                  %   March 24, 2016
                end                                                                       %   March 24, 2016
            else
                AP=false(size(r));
                DAPI=false(size(r));
            end
            
            if trap_flag==1
                check_file=dir([thresh_dir,c1,'_3',mb_com]);
                if ~isempty(check_file)
                    eval(['TRAP=imread(''',thresh_dir,c1,'_3',mb_com,''',''jpg'');']);
                    if size(TRAP,3)==3
                        TRAP=TRAP(:,:,1);
                    end
                    TRAP=(TRAP>127);
                    TRAP=imrotate(TRAP,ang);
                else
                    TRAP=false(size(a));
                end
            else
                TRAP=false(size(r));
            end
           
%              a=a(:,:,1);
%              a=(a>threshold(bone_number(b_n),s,3));
            if ap_flag==1
                [r, g, AP]=clear_beads_1(file_name, a,r,g,AP);
            else
                [r, g, AP]=clear_beads_1(file_name, a,r,g);
            end

            no_flag=0;
%         elseif s==4
%             c=phr(b_t,:);
%             c1=[c,'\',c,'_',nnn,'FL\',c,'_',nnn,'FL.jpg_Files\',c,'_',nnn,'FL_c'];
%             eval(['a=imread(''',direc,c1,'0',mb_com,''',''jpg'');']);
%             eval(['g=imread(''',direc,c1,'1',mb_com,''',''jpg'');']);
%             %eval(['y=imread(''',direc,c1,'2.jpg'',''jpg'');']);
%             eval(['r=imread(''',direc,c1,'2',mb_com,''',''jpg'');']);
%             y=uint8(zeros(size(r)));
% 
%             no_flag=0;
%             end
        %end
%        a=imrotate(a(max_dis(2)+1:end,max_dis(1)+1:end),rot_angle);
%        %a=a(dis_y+1:end,dis_x+1:end);
%        g=imrotate(g(max_dis(2)-dis(1,2)+1:end-(dis(1,2)),max_dis(1)-dis(1,1)+1:end-(dis(1,1))),rot_angle);
%        %y=imrotate(y(max_dis(2)-dis(2,2)+1:end-(dis(2,2)),max_dis(1)-dis(2,1)+1:end-(dis(2,1))),rot_angle);
%        r=imrotate(r(max_dis(2)-dis(2,2)+1:end-(dis(2,2)),max_dis(1)-dis(2,1)+1:end-(dis(2,1))),rot_angle); % red beads
%         a=imrotate(a,rot_angle);
        %a=a(dis_y+1:end,dis_x+1:end);
%         g=imrotate(g,rot_angle);
        if GFP_flag==1                                                        %   March 24, 2016
%             [g, G]=separate_GFP(g,a);       % when GFP is green channel     %   March 24, 2016%
            eval(['G=imread(''',thresh_dir,c1,'_8',mb_com,''',''jpg'');']);   %   March 24, 2016
            G=G>127;                                                          %   March 24, 2016
            G=imrotate(G,ang);                                                %   March 24, 2016
            G=xor(G, g|r) & (G & ~(g|r));                                     %   March 24, 2016
        else                                                                  %   March 24, 2016
            G=false(size(g));                                                 %   March 24, 2016
        end                                                                   %   March 24, 2016
        %y=imrotate(y(max_dis(2)-dis(2,2)+1:end-(dis(2,2)),max_dis(1)-dis(2,1)+1:end-(dis(2,1))),rot_angle);
%         r=imrotate(r,rot_angle); % red beads
        %y=imrotate(y,rot_angle); % red beads
        %r=imsubtract(y,g);      % tetracycline 



%if b_n==1 | b_n==2 | b_n==4 |b_n==5 | b_n==6 | b_n==8 | b_n==9 | b_n==11 | b_n==12 | b_n==13 | b_n==14
%    start_point=start_point+1300/2+dis_10x;
%elseif b_n==3 | b_n==10
%    start_point=start_point+850/2+dis_10x;
%elseif b_n==7
%    if s==1
%        start_point=start_point+850/2+dis_10x;
%    else
%        start_point=start_point+1300/2+dis_10x;
%    end
%end

%eval(['start_point=start_point+bottom_growth_plate_',num2str(nnn),'(',num2str(s),')/ratio+dis_10x;'])
%eval(['start_point=start_point+1950/2;'])%bottom_growth_plate_',num2str(nnn),'(',num2str(s),')/ratio;'])
% start_point=round(end_growth_plate(b_n,s)/ratio);
[start_point]=find_points_after_rotate(ori_size, size(a,1), end_growth_plate(b_n,s), ang);
start_point=round(start_point/ratio);
% end_point=round(end_points(b_n,s)/ratio);             % 4/16/2008
[end_point]=find_points_after_rotate(ori_size, size(a,1), end_points(b_n,s), ang);
end_point=round(end_point/ratio);             % 4/16/2008


% if auto_exp==1return
%     a=(a>(graythresh(a)*255-std(double(a(:)))));
% else
%     a=(a>threshold(bone_number(b_n),s,3));
% end
a=imresize(a,1/ratio);
a=a(start_point:end_point,:,:);

AP=imresize(AP,1/ratio);
AP=AP(start_point:end_point,:,:);

DAPI=imresize(DAPI,1/ratio);                                                  %   March 24, 2016
DAPI=DAPI(start_point:end_point,:,:);                                         %   March 24, 2016

TRAP=imresize(TRAP,1/ratio);
TRAP=TRAP(start_point:end_point,:,:);

G=imresize(G,1/ratio);                                                        %   March 24, 2016
G=G(start_point:end_point,:,:);                                               %   March 24, 2016






%G=G/max(G(:));

%a=a(start_point:end,:,:);
% figure;imshow(a)
% a=double(sum(a,3));
% a=a/max(a(:));


%eval(['a=a.^(1/1-beta_',num2str(nnn),'(',num2str(s),'));'])
% a=a/max(a(:));

%tic
%a=a-medfilt2(a, [100 100]);
%a=a-min(a(:));
%a=a/max(a(:));
%toc

%a=a(1630:3500,500:2550);


%am=medfilt2(a,[5,5]);
%%am=medfilt2(a,[11,11]);
%figure;imshow(am)


%figure;imshow(a)

%r=imread('C:\Users\shhong\Desktop\Jake New\11473-1\IP-ConvertImage-06\IP-ConvertImage-06_c3_1.JPG','jpg');
%   r=imread('C:\Users\shhong\Desktop\Jake_122607\Jake_5.jpg_Files\Jake_5_c1_1.jpg','jpg');
%r=imread('E:\Jake New\11494\IP-ConvertImage-01\IP-ConvertImage-01_c3_1.JPG','jpg');
%r=imread('E:\Jake New\11473-1\IP-ConvertImage-06\IP-ConvertImage-06_c3_1.JPG','jpg');

%if rot_angle~=0
%    r=imrotate(r,rot_angle);
%end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Search best threshold
%   January 7, 2009
%
%r_temp=r(:,round(size(r,2)/2)-20:round(size(r,2)/2)+20);
%
%%%%%%%%%%%%%
%%   AAA
%%%%%%%%%%%%%
%
%rr=(r_temp>mean(r_temp(:)));
%r1=r_temp.*double(rr);
%[I J V]=find(r1);
%%if b_t==1 & (b_n==1 | b_n==4 | b_n==5)
%%    r=(r>graythresh(V)+std(V(:))*3);
%%elseif b_t==1 & (b_n==2)
%%    r=(r>graythresh(V)+std(V(:))*(-2));
%%elseif b_t==2 & (b_n==1 | b_n==3 | b_n==4 | b_n==5 | b_n==6 | b_n==7)
%%    r=(r>graythresh(V)+std(V(:))*(-1.5));
%%else
%%    r=(r>graythresh(V)+std(V(:))*0);
%%end    
%%r=(r>graythresh(r_temp));
%%r=(r>graythresh(V)+std(V(:))*-1.5);


% r=(r>threshold(bone_number(b_n),s,2));
r=imresize(r,1/ratio);
r=r(start_point:end_point,:,:);
%%r=r(start_point:end,:,:);
%figure;imshow(r)
%r=r(:,:,1);
% r=double(r);
% r=r/max(r(:));

a=a|r;
number_r=length(find(r>0));
% a=double(a>0);

[L n_a]=bwlabel(a);
stats_a=regionprops(L,'Area','Image','BoundingBox'); clear L
a=false(size(a));
for i=1:n_a
%     temp=zeros(size(a));
    start_x=round(stats_a(i).BoundingBox(2));
    start_y=round(stats_a(i).BoundingBox(1));
    end_x=start_x+stats_a(i).BoundingBox(4)-1;
    end_y=start_y+stats_a(i).BoundingBox(3)-1;
    temp=stats_a(i).Image;
    area_temp=stats_a(i).Area;
    temp_r=temp & r(start_x:end_x,start_y:end_y);
    if area_temp==length(find(temp_r(:)))
        r(start_x:end_x,start_y:end_y)=xor(r(start_x:end_x,start_y:end_y),temp_r) & (r(start_x:end_x,start_y:end_y) & ~temp_r);
        continue
    else
        a(start_x:end_x,start_y:end_y)=a(start_x:end_x,start_y:end_y) | stats_a(i).Image;
    end
end
% r=double(r>0);
% a=double(a>0);

% for i=start_threshold_ratio:15
%     r1=(r>graythresh(r)+std(r(:))*(i/10));
%     threshold_ratio(i-start_threshold_ratio+1)=sum(r1(:))/(size(r1,1)*size(r1,2));
% end
% double_differentiation=abs(diff(diff(threshold_ratio)));
% start_min_point=min(find(double_differentiation==max(double_differentiation)));
% index=min(start_threshold_ratio+start_min_point+(min(find(double_differentiation(start_min_point:end)<segment_threshold & ...
%     double_differentiation(start_min_point:end)>0.00005))-1));
% if isempty(index)==1
%     index=min(find(double_differentiation(start_min_point:end)==min(double_differentiation(start_min_point:end))))+start_min_point+2;
% end
% r1=(r>graythresh(r)+std(r(:))*(index/10));
% r=r1; clear r1
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%r=r(1630:3500,500:2550);
%g=imread('C:\Users\shhong\Desktop\Jake New\11473-1\IP-ConvertImage-06\IP-ConvertImage-06_c2_1.JPG','jpg');
%   g=imread('C:\Users\shhong\Desktop\Jake_122607\Jake_5.jpg_Files\Jake_5_c0_2.jpg','jpg');
%g=imread('E:\Jake New\11494\IP-ConvertImage-01\IP-ConvertImage-01_c2_1.JPG','jpg');
%g=imread('E:\Jake New\11473-1\IP-ConvertImage-06\IP-ConvertImage-06_c2_1.JPG','jpg');

%if rot_angle~=0
%    g=imrotate(g,rot_angle);
%end

%g=(g>threshold(bone_number(b_n),s,2));

g=imresize(g,1/ratio);
g=g(start_point:end_point,:,:);
%g=g(start_point:end,:,:);
%g=g(:,:,2);
% figure;imshow(g)
% g=double(g);
% g=g/max(g(:));

a=a|g;
% at=double(a>0);
at=a;

[L n_a]=bwlabel(at);
stats_a=regionprops(L,'Area','Image','BoundingBox'); clear L
a=false(size(at));
for i=1:n_a
%     temp=zeros(size(at));
    start_x=round(stats_a(i).BoundingBox(2));
    start_y=round(stats_a(i).BoundingBox(1));
    end_x=start_x+stats_a(i).BoundingBox(4)-1;
    end_y=start_y+stats_a(i).BoundingBox(3)-1;
    temp=stats_a(i).Image;
    area_temp=stats_a(i).Area;
    temp_g=temp & g(start_x:end_x,start_y:end_y);
    if area_temp==length(find(temp_g(:)))
        g(start_x:end_x,start_y:end_y)=xor(g(start_x:end_x,start_y:end_y), temp_g) & (g(start_x:end_x,start_y:end_y) & ~temp_g);
        continue
    else
        a(start_x:end_x,start_y:end_y)=a(start_x:end_x,start_y:end_y) | stats_a(i).Image;
    end
end
% g=double(g>0);

% a=remove_DAPI_from_DIC(a,DAPI);                   % added on October 24, 2011
a=remove_unlabelled_DIC(a,r|g);          % added on October 24, 2011

% a=double(a>0);
%g=g-medfilt2(g, [30 30]);

number_g=length(find(g>0));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Search best threshold
%   January 7, 2009
%

%g_temp=g(:,round(size(g,2)/2)-20:round(size(g,2)/2)+20);

%gg=(g_temp>mean(g_temp(:)));
%g1=g_temp.*double(gg);
%[I J V]=find(g1);

% for i=start_threshold_ratio:15
%     g1=(g>graythresh(g)+std(g(:))*(i/10));
%     threshold_ratio(i-start_threshold_ratio+1)=sum(g1(:))/(size(g1,1)*size(g1,2));
% end
% double_differentiation=abs(diff(diff(threshold_ratio)));
% start_min_point=min(find(double_differentiation==max(double_differentiation)));
% index=min(start_threshold_ratio+start_min_point+(min(find(double_differentiation(start_min_point:end)<segment_threshold & ...
%     double_differentiation(start_min_point:end)>0.00005))-1));
% if isempty(index)==1
%     index=min(find(double_differentiation(start_min_point:end)==min(double_differentiation(start_min_point:end))))+start_min_point+2;
% end
% g1=(g>graythresh(g)+std(g(:))*(index/10));
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %g=(g>graythresh(V)+std(V(:))*-1.5);
% 
% g=g1;
% %g1=(g>graythresh(g)-std(g(:))*std_no_green);

%if length(size(threshold,3))>3





%aa=a;
%
%           red = TRITC + mCherry
%
%           October 9 2009
%
%
% y=imresize(y,1/ratio);
% y=y(start_point:end_point,:,:);
% %r=r(start_point:end,:,:);
% figure;imshow(y)
% y=y(:,:,1);
% y=double(y);
% y=y/max(y(:));
% aa=a+y;
% 
% %r=r-medfilt2(r, [30 30]);
% 
% number_y=length(find(y>0.5));
% %r=(r>graythresh(r)-std(r(:))*std_no_red);
% for i=start_threshold_ratio:15
%     y1=(y>graythresh(y)+std(y(:))*(i/10));
%     threshold_ratio(i-start_threshold_ratio+1)=sum(y1(:))/(size(y1,1)*size(y1,2));
% end
% double_differentiation=abs(diff(diff(threshold_ratio)));
% start_min_point=min(find(double_differentiation==max(double_differentiation)));
% index=min(start_threshold_ratio+start_min_point+(min(find(double_differentiation(start_min_point:end)<segment_threshold & ...
%     double_differentiation(start_min_point:end)>0.00005))-1));
% if isempty(index)==1
%     index=min(find(double_differentiation(start_min_point:end)==min(double_differentiation(start_min_point:end))))+start_min_point+2;
% end
% y1=(y>graythresh(y)+std(y(:))*(index/10));
% y=y1; clear y1
% 
% y=y-g;
% g=g-y;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%
%   For tetracycline
%   Dec. 8 2008

%r=r-(r&g);

number_r=length(find(r>0));
%
%

% if b_t==1
%     temp_color=r;
%     number_temp=number_r;
%     %r=g;
% %    g=y;
%     %number_r=number_g;
% %    number_g=number_y;
% elseif b_t==2
%     temp_color=r;
%     number_temp=number_r;
% %    r=y;
% %    number_r=number_y;
% elseif b_t==4
%     temp_color=r;
%     number_temp=number_r;
%     r=g;
%     g=temp_color;
%     number_r=number_g;
%     number_g=number_temp;   
% end
g1=g;

if (number_r/number_g)>300
    g=false(size(g));
    g1=g;
    no_green=1;
else
    no_green=0;
end

if (number_g/number_r)>100
    r=false(size(r));
    r1=r;
    no_red=1;
else
    no_red=0;
end

%g=g(1630:3500,500:2550);






%[a, exact]=get_rid_of_bubbles(a, bubble_threshold);

%%ar=ar-ar.*imdilate(exact,strel('disk',5));
%r=r-exact;
%g=g-exact;






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Search best threshold
%   January 7, 2009
%

%a1=ar.^0.5;
%%a1=a.^1.5;
%at=double(a1>graythresh(a1)*max(a1(:)));

% for i=start_threshold_ratio:15
%     a1=(aa>graythresh(aa)+std(aa(:))*(i/10));
%     threshold_ratio(i-start_threshold_ratio+1)=sum(a1(:))/(size(a1,1)*size(a1,2));
% end
% double_differentiation=abs(diff(diff(threshold_ratio)));
% start_min_point=min(find(double_differentiation==max(double_differentiation)));
% index=min(start_threshold_ratio+start_min_point+(min(find(double_differentiation(start_min_point:end)<segment_threshold & ...
%     double_differentiation(start_min_point:end)>0.00005))-1));
% if isempty(index)==1
%     %index=min(find(double_differentiation(start_min_point:end)==min(double_differentiation(start_min_point:end))))+start_min_point+2;
%     if isempty(max(find(threshold_ratio>=0.06)))==1
%         index=start_threshold_ratio;
%     else
%         index=start_threshold_ratio+max(find(threshold_ratio>=0.06));        
%     end
% end


%at=(aa>graythresh(aa)+std(aa(:))*(index/10))+g+r;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % at=(aa>graythresh(aa))+g+r;
% arg1=a;%+r+double(g);
% % arg2=double(arg1>=1); clear arg1
% arg2=double(arg1>0); clear arg1
% %arg=arg1.*arg2+arg1.*(1-arg2);
% arg=arg2/max(arg2(:)); clear arg2
arg=a;

%arg=(a.^1.5)*2+r+g;
%arg=arg/max(arg(:));
% figure;imshow(arg)
% argm=imerode(imdilate(arg, strel('disk',5)),strel('disk',5)); clear arg
argm=imerode(imdilate(arg, strel('disk',3)),strel('disk',3)); clear arg
%argm=medfilt2(arg,[5,5]); %clear a ar arg
% argm=argm/max(argm(:));
% % figure;imshow(argm)
% argm=double(argm>0);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       For test Jan. 8 2009

%%ar=ar+r;
%ar=a+r;
%ar=ar/max(ar(:));
%figure;imshow(ar)
%%arm=medfilt2(ar,[5,5]); %clear a ar arg
%arm=imerode(imdilate(ar,strel('disk',5)),strel('disk',5));
%arm=arm/max(arm(:));
%figure;imshow(arm)
%
%
%%
%%   3/13
%%
%%
%if with_rg_label==1                 % using labels
%    argm=(argm>graythresh(argm));
%    argm=((argm|g1|r)-(exact>0));
%    %argm=imerode(imdilate(argm, strel('disk',5)),strel('disk',5));
%elseif with_rg_label==0             % not using labels
%    argm=(am>graythresh(am));
%end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%argm=medfilt2(argm, [9 9]);
% argm=medfilt2(argm, [13 13]);%+medfilt2(g, [15 15])+medfilt2(r, [15 15]);

[L n_argm]=bwlabel(argm);
stats_argm=regionprops(L,'Area','Image','FilledImage','BoundingBox','Orientation'); clear L
armt=false(size(argm));
for i=1:n_argm
    start_x=round(stats_argm(i).BoundingBox(2));
    start_y=round(stats_argm(i).BoundingBox(1));
    end_x=start_x+stats_argm(i).BoundingBox(4)-1;
    end_y=start_y+stats_argm(i).BoundingBox(3)-1;
    if length(find(stats_argm(i).FilledImage))-length(find(stats_argm(i).Image))<100
        armt(start_x:end_x,start_y:end_y)=armt(start_x:end_x,start_y:end_y)|stats_argm(i).FilledImage;
    else
        armt(start_x:end_x,start_y:end_y)=armt(start_x:end_x,start_y:end_y)|stats_argm(i).Image;
    end
end

% G=(1-armt)&g1;

g=xor(g1, G) & (g1 & ~G);


%
%   Added on 3/13
%
% G=zeros(size(g));
%g=g-r;             % June 4, 2008
% gg=(g>0);
% g=g.*gg;%-min(g(:));
% if max(g(:))~=0
%     g=g/max(g(:));
% end
%   Added on 3/13



%%g2=imread('C:\Users\shhong\Desktop\Jake New\11473-1\IP-ConvertImage-06\IP-ConvertImage-06_c2_2.JPG','jpg');
%g2=imread('C:\Users\shhong\Desktop\Jake_121807\Jake_10\Jake_10_c1_2.jpg','jpg');
%%g2=imread('E:\Jake New\11494\IP-ConvertImage-01\IP-ConvertImage-01_c2_2.JPG','jpg');
%%g2=imread('E:\Jake New\11473-1\IP-ConvertImage-06\IP-ConvertImage-06_c2_2.JPG','jpg');
%g2=g2(587:end,:,:);
%g2=g2(:,:,2);
%g2=double(g2);
%g2=g2/max(g2(:));
%if rot_angle~=0
%    g2=imrotate(g2,rot_angle);
%end
%%g2=g2(1630:3500,500:2550);
%g2t=(g2>graythresh(g2));
%fclose all

%argm1=(argm>graythresh(argm));

%[L n_argm]=bwlabel(argm1);
%stats_argm=regionprops(L,'Area','Image','FilledImage','BoundingBox','Orientation'); clear L
%argmt=zeros(size(argm));
%for i=1:n_argm
%    start_x=round(stats_argm(i).BoundingBox(2));
%    start_y=round(stats_argm(i).BoundingBox(1));
%    end_x=start_x+stats_argm(i).BoundingBox(4)-1;
%    end_y=start_y+stats_argm(i).BoundingBox(3)-1;
%    argmt(start_x:end_x,start_y:end_y)=argmt(start_x:end_x,start_y:end_y)+stats_argm(i).FilledImage;
%end

%argmtg2t=(1-argmt)&g2t;
%figure;imshow(argmtg2t)

%g=g2t-argmtg2t;


%argm_argmtg2t=argm+argmtg2t;
%argm_argmtg2t=argmt+argmtg2t;
%argm_argmtg2t=argm_argmtg2t/max(argm_argmtg2t(:));

%figure;imshow(argm_argmtg2t)

%figure;imshow(g2)


% armt=armt+r+g;
% 
% armt=imerode(imdilate(armt,strel('disk',4)), strel('disk',4));


t=uint8(armt)*255; clear am
t(:,:,2)=uint8((g+G))*255;
t(:,:,3)=uint8(0);
% figure;imshow(t)


%t(:,:,1)=uint8(armt*255);
%figure;imshow(t); clear t
%axis on
%grid

%argmt=(argm>graythresh(argm));

[L n]=bwlabel(armt);
stats=regionprops(L,'Area','Image','FilledImage','BoundingBox','Orientation'); clear L

result=false(size(armt));

area=[stats.Area];
s_area=sort(area);
[I J]=find(area/sum(area)>=thresh);
for i=1:length(J)
    start_x=round(stats(J(i)).BoundingBox(2));
    start_y=round(stats(J(i)).BoundingBox(1));
    end_x=start_x+stats(J(i)).BoundingBox(4)-1;
    end_y=start_y+stats(J(i)).BoundingBox(3)-1;
    result(start_x:end_x,start_y:end_y)=result(start_x:end_x,start_y:end_y)|stats(J(i)).Image;
end

[I J]=find(area/sum(area)<thresh);
for i=1:length(J)
    start_x=round(stats(J(i)).BoundingBox(2));
    start_y=round(stats(J(i)).BoundingBox(1));
    end_x=start_x+stats(J(i)).BoundingBox(4)-1;
    end_y=start_y+stats(J(i)).BoundingBox(3)-1;
    if no_green==0 & no_red==0
        tmp=stats(J(i)).Image & (r(start_x:end_x,start_y:end_y) | g(start_x:end_x,start_y:end_y));
    elseif no_green==0 & no_red==1
        tmp=stats(J(i)).Image & (g(start_x:end_x,start_y:end_y));
    elseif no_green==1 & no_red==0
        tmp=stats(J(i)).Image & (r(start_x:end_x,start_y:end_y));
    elseif no_green==1 & no_red==1
        error('No green label or red label');
    end
%    if ~isempty(find(tmp==1))       % Matlab 6.5
    if ~isempty(find(tmp,1))       % Matlab 7.1
        result(start_x:end_x,start_y:end_y)=result(start_x:end_x,start_y:end_y)|stats(J(i)).Image;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Find broken cortical bones
%
%   March 5, 2008
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%II(1)=find(area==s_area(length(s_area)));
%II(2)=find(area==s_area(length(s_area)-1));
%
%if stats(II(1)).BoundingBox(1)<stats(II(1)).BoundingBox(1)
%    left_cortical=stats(II(1)).FilledImage;
%    right_cortical=stats(II(2)).FilledImage;
%else
%    left_cortical=stats(II(2)).FilledImage;
%    right_cortical=stats(II(1)).FilledImage;
%end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%armt1=armt;
%%armt(150-2:150+2,:)=zeros(5,size(armt,2));
%armt=armt(dis_10x+1:end,:);

cortical=false(size(armt));
periosteum=false(size(armt));
endosteum=false(size(armt));

%tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Find left cortical bone (3/6/08)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%third=round(size(armt,2)/3);
% third=round(size(armt,2)/2);
[third]=find_middle_point_of_bone(a);
armt1=imclose(armt, strel('disk',15));
[L n]=bwlabel(armt1(1:end,1:third));
stats=regionprops(L,'Area','Image','FilledImage','BoundingBox','Orientation'); clear L

area=[stats.Area];
s_area=sort(area);

h_left=0;
left_overlap=zeros(size(armt));
clear II
for i=1:length(s_area)
    temp=find(area==s_area(length(s_area)-(i-1)));
    for j=1:length(temp)
        II(i)=temp(j);
        i=i+1;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   find cortical that is 10% of total area and leftmost
%
%   April 07, 2011
%
% i=1;
temp=find(area/sum(area)*100>10);
for i=1:length(temp)
    st(i)=round(stats(temp(i)).BoundingBox(1));
end

% i=temp(find(st==min(st)));
% i=(find(st==min(st)));
% i=i(1);
% sum_x=0;
old=0;
ii=1;
for i=1:length(temp)
%     st(i)=round(stats(II(temp(i))).BoundingBox(1));
    st(i)=round(stats(temp(i)).BoundingBox(1));
    
    start_x(i)=round(stats(temp(i)).BoundingBox(2));
    start_y(i)=round(stats(temp(i)).BoundingBox(1));
    end_x(i)=start_x(i)+stats(temp(i)).BoundingBox(4)-1;
    end_y(i)=start_y(i)+stats(temp(i)).BoundingBox(3)-1;
    [X Y]=find(stats(temp(i)).Image);
    test_center(i)=median(Y)+start_y(i);
    if stats(temp(i)).BoundingBox(4)>old
        main_center=test_center(i);
        old=test_center(i);
        ii=i;
    end
%     sum_x=sum_x+stats(temp(i)).BoundingBox(4);
end


% i=temp(find(st==max(st)));
% i=(find(st==min(st)));
% i=i(1);
i=ii;
cortical(start_x(i):end_x(i),start_y(i):end_y(i))=cortical(start_x(i):end_x(i),start_y(i):end_y(i))|stats(temp(i)).Image;
sum_x=stats(temp(i)).BoundingBox(4);

for t_i=1:length(temp)
    if i==t_i
%         sum_x=sum_x+stats(temp(t_i)).BoundingBox(4);
        continue
    end
    if test_center(t_i)-main_center<350
        if start_x(t_i)>=start_x(i) & end_x(t_i)<=end_x(i)
            continue
        end
        cortical(start_x(t_i):end_x(t_i),start_y(t_i):end_y(t_i))=cortical(start_x(t_i):end_x(t_i),start_y(t_i):end_y(t_i))|stats(temp(t_i)).Image;
        sum_x=sum_x+stats(temp(t_i)).BoundingBox(4);
    end
end


% start_x(i)=round(stats(II(i)).BoundingBox(2));
% start_y(i)=round(stats(II(i)).BoundingBox(1));
% end_x(i)=start_x(i)+stats(II(i)).BoundingBox(4)-1;
% end_y(i)=start_y(i)+stats(II(i)).BoundingBox(3)-1;
% [X Y]=find(stats(II(1)).Image);

% test_center(i)=mean(Y)+start_y(i);
% main_center=test_center(1);
% 
% %cortical(start_x:end_x,start_y:end_y)=cortical(start_x:end_x,start_y:end_y)+stats(II(i)).FilledImage;
% cortical(start_x:end_x,start_y:end_y)=cortical(start_x:end_x,start_y:end_y)+stats(II(i)).Image;
%
% start_x(i)=round(stats(i).BoundingBox(2));
% start_y(i)=round(stats(i).BoundingBox(1));
% end_x(i)=start_x(i)+stats(i).BoundingBox(4)-1;
% end_y(i)=start_y(i)+stats(i).BoundingBox(3)-1;
% [X Y]=find(stats(i).Image);
% test_center=mean(Y)+start_y(i);
% main_center=test_center;
% 
% %cortical(start_x:end_x,start_y:end_y)=cortical(start_x:end_x,start_y:end_y)+stats(II(i)).FilledImage;
% cortical(start_x(i):end_x(i),start_y(i):end_y(i))=cortical(start_x(i):end_x(i),start_y(i):end_y(i))+stats(i).Image;
% %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% if sum_x/size(cortical,1)<0.95
%     h_left=h_left+(end_x(i)-start_x(i)+1);
%     left_overlap(start_x(i):end_x(i),start_y(i):end_y(i))=left_overlap(start_x(i):end_x(i),start_y(i):end_y(i)) + ones(end_x(i)-start_x(i)+1,end_y(i)-start_y(i)+1);
% 
%     t_st_x=start_x(1)+100;
%     t_ed_x=end_x(1)-100;
%     t_st_y=max(1,start_y(1)-100);
%     t_ed_y=min(end_y(1)+100,third);
%     temp_im=armt(1:end,1:third);
%     t_im=zeros(size(temp_im));
%     t_im(:,t_st_y:t_ed_y)=ones(size(t_im,1),t_ed_y-t_st_y+1);
%     t_im(t_st_x:t_ed_x,t_st_y:t_ed_y)=zeros(t_ed_x-t_st_x+1,t_ed_y-t_st_y+1);
%     temp_im=temp_im&t_im;
%     temp_im=temp_im-(temp_im&cortical(:,1:third));
%     [L n]=bwlabel(temp_im);
%     stats=regionprops(L,'Area','Image','FilledImage','BoundingBox','Orientation'); clear L
%     clear start_x, end_x, start_y, end_y
% 
%     area=[stats.Area];
%     s_area=sort(area);
%     nn=find(s_area/sum(s_area)>0.005);
%     clear II
%     for i=1:length(nn)
%         temp=find(area==s_area(length(s_area)-(i-1)));
%         for j=1:length(temp)
%             II(i)=temp(j);
%             i=i+1;
%         end
%     end
% 
% 
%     for i=1:length(nn)
%         start_x(i)=round(stats(II(i)).BoundingBox(2));
%         start_y(i)=round(stats(II(i)).BoundingBox(1));
%         end_x(i)=start_x(i)+stats(II(i)).BoundingBox(4)-1;
%         end_y(i)=start_y(i)+stats(II(i)).BoundingBox(3)-1;
%         [X Y]=find(stats(II(i)).Image);
%     %    test_center(i)=round((end_y(i)-start_y(i))/2)+start_y(i);
%         test_center(i)=mean(Y)+start_y(i);
%     end
%     center=abs(test_center-main_center);
%     s_center=sort(center(1:end-1));
%     j=0;
%     while h_left<size(armt,1)*1.1
%         j=j+1;
%         if j>numel(s_center)
%             break
%         end
%         ii=find(center==s_center(j));
%     %    if ~(start_x(1)<start_x(ii) & end_x(1)>end_y(ii))
%         if length(ii)>1
%             ii=ii(1);
%         end
%         if ~(start_x(1)<start_x(ii) & end_x(1)>end_x(ii))
%     %        cortical(start_x(ii):end_x(ii),start_y(ii):end_y(ii))=cortical(start_x(ii):end_x(ii),start_y(ii):end_y(ii))+stats(II(ii)).FilledImage;
%             cortical(start_x(ii):end_x(ii),start_y(ii):end_y(ii))=cortical(start_x(ii):end_x(ii),start_y(ii):end_y(ii))+stats(II(ii)).Image;
%             h_left=h_left+(end_x(ii)-start_x(ii)+1);
%             left_overlap(start_x(ii):end_x(ii),start_y(ii):end_y(ii))=left_overlap(start_x(ii):end_x(ii),start_y(ii):end_y(ii))...
%                 + ones(end_x(ii)-start_x(ii)+1,end_y(ii)-start_y(ii)+1);
%         end
%     end
% end

% [L n]=bwlabel(cortical(1:end,1:third));
[L n]=bwlabel(cortical);
stats=regionprops(L,'Area','Image','FilledImage','BoundingBox','Orientation'); clear L
for i=1:n
   t_order(i)=stats(i).BoundingBox(2);
end
s_order=sort(t_order);
p=0;
for i=1:n
    temp=find(t_order==s_order(i));
    if length(temp)>1
        p=p+1;
        order(i)=temp(p);
        if p>=length(temp)
            p=0;
        end
    else
        p=0;
        order(i)=temp;
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   remove small bones which are overlapped with cortical bone in vertical direction
%   4/16/08
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

temp=round([stats.BoundingBox]);
temp=reshape(temp',4,length(temp)/4);
%temp=[temp(2,:);temp(2,:)+temp(4,:)];
temp=[temp(2,:);temp(4,:)];
temp(2,:)=temp(1,:)+temp(2,:)-1;


%[I V]=sort(temp(2,:),2);
s_temp=sort(temp,2);
%for i=1:length(V)
%   s_temp(:,i)=temp(:,V(i));
%end


p=1; j=1;
i=1;
clear s_temp1
while i<size(s_temp,2)+1
    if j>1
        i=i+j-1;
        j=1;
    end
    if i>size(s_temp,2)
        break
    end
    x=find(temp(1,:)==s_temp(1,i));
    for j=1:length(x)
        s_temp1(p)=temp(2,x(j));
        p=p+1;
    end
    i=i+1;
end
s_temp(2,:)=s_temp1;
x=100;
size_s_temp=size(s_temp,2);
while isempty(x)==0
    d_temp=diff(s_temp(2,:),1,2);
    x=find(d_temp<=0);
    if isempty(x)==1
        break
    end
    s_temp(:,x(1)+1:size_s_temp-1)=s_temp(:,x(1)+2:size_s_temp);
    size_s_temp=size_s_temp-1;
    s_temp=s_temp(:,1:size_s_temp);
end

clear order
for i=1:size(s_temp,2)
    x=find(sum(temp-s_temp(:,i)*ones(1,size(temp,2)),1)==0);
    %order(i)=x;    % Dec. 9.2008
    order(i)=min(x);
end
    
    
    
    
%for i=1:n
%   t_order(i)=stats(i).BoundingBox(2);
%end
%s_order=sort(t_order);
%for i=1:n
%    order(i)=find(t_order==s_order(i));
%end

if length(order)>1
    %%%%%%%%%%%%%%%%%%%%%%%%
    %
    %   left broken bone (3/5/08)
    %
    %%%%%%%%%%%%%%%%%%%%%%%%
    for i=1:length(order)-1
        outline1=false(size(cortical)); outline2=outline1;
        for j=1:2
            temp_im=stats(order(i-1+j)).Image;
            start_x(j)=round(stats(order(i-1+j)).BoundingBox(2));
            start_y(j)=round(stats(order(i-1+j)).BoundingBox(1));
            end_x(j)=start_x(j)+stats(order(i-1+j)).BoundingBox(4)-1;
            end_y(j)=start_y(j)+stats(order(i-1+j)).BoundingBox(3)-1;
            eval(['outline',num2str(j),'(start_x(j):end_x(j),start_y(j):end_y(j))=outline',num2str(j),'(start_x(j):end_x(j),start_y(j):end_y(j))|temp_im;'])
            eval(['temp=imerode(outline',num2str(j),',strel(''disk'',1));'])
            eval(['outline',num2str(j),'=xor(outline',num2str(j),',temp) & (outline',num2str(j),' & ~temp);'])
        end
        [x y d]=find_closest_points(outline1, outline2);
        
        x1=x(1);
        y1=y(1);
        x2=x(2);
        y2=y(2);
        a=(x2-x1)/(y2-y1);
        
        if y1==y2 & y1>1
            cortical(min(x1,x2):max(x1,x2),y1-1:y1+1)=true(abs(x2-x1)+1,3);
        elseif x1==x2 & x1>1
            cortical(x1-1:x1+1,min(y1,y2):max(y1,y2))=true(3,abs(y2-y1)+1);
        elseif abs(a)<1
            b=x1-a*y1;
            if y1<y2
                for y=y1:y2
                    cortical(max(floor(a*y+b)-1,1):floor(a*y+b)+1,y-1:y+1)=true(floor(a*y+b)+1-max(floor(a*y+b)-1,1)+1,3);
                    cortical(max(ceil(a*y+b)-1,1):ceil(a*y+b)+1,y-1:y+1)=true(ceil(a*y+b)+1-max(ceil(a*y+b)-1,1)+1,3);
                end
            else
                for y=y1:-1:y2
                    cortical(max(floor(a*y+b)-1,1):floor(a*y+b)+1,y-1:y+1)=true(floor(a*y+b)+1-max(floor(a*y+b)-1,1)+1,3);
                    cortical(max(ceil(a*y+b)-1,1):ceil(a*y+b)+1,y-1:y+1)=true(ceil(a*y+b)+1-max(ceil(a*y+b)-1,1)+1,3);
                end
            end
        elseif abs(a)>=1
            a=(y2-y1)/(x2-x1);
            b=y1-a*x1;
            if x1<x2
                for x=x1:x2
                    cortical(x-1:x+1, max(floor(a*x+b)-1,1):floor(a*x+b)+1)=true(3,floor(a*x+b)+1-max(floor(a*x+b)-1,1)+1);
                    cortical(x-1:x+1, max(ceil(a*x+b)-1,1):ceil(a*x+b)+1)=true(3,ceil(a*x+b)+1-max(ceil(a*x+b)-1,1)+1);
                end
            else
                for x=x1:-1:x2
                    cortical(x-1:x+1, max(floor(a*x+b)-1,1):floor(a*x+b)+1)=true(3,floor(a*x+b)+1-max(floor(a*x+b)-1,1)+1);
                    cortical(x-1:x+1, max(ceil(a*x+b)-1,1):ceil(a*x+b)+1)=true(3,ceil(a*x+b)+1-max(ceil(a*x+b)-1,1)+1);
                end
            end
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Find right cortical bone (3/6/08)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h_right=0;
right_overlap=false(size(armt));
% r_third=size(armt,2)-third;
r_third=third;
armt1=armt1&(1-cortical);

[L n]=bwlabel(armt1(1:end,r_third+1:end));
stats=regionprops(L,'Area','Image','FilledImage','BoundingBox','Orientation'); clear L

area=[stats.Area];
s_area=sort(area);
clear II
for i=1:length(s_area)
    temp=find(area==s_area(length(s_area)-(i-1)));
    for j=1:length(temp)
        II(i)=temp(j);
        i=i+1;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   find cortical that is 10% of total area and rightmost
%
%   April 07, 2011
%
% i=1;
temp=find(area/sum(area)*100>10);
clear st
sum_x=0;
old=0;
ii=1;
for i=1:length(temp)
%     st(i)=round(stats(II(temp(i))).BoundingBox(1));
    st(i)=round(stats(temp(i)).BoundingBox(1));
    
    start_x(i)=round(stats(temp(i)).BoundingBox(2));
    start_y(i)=round(stats(temp(i)).BoundingBox(1));
    end_x(i)=start_x(i)+stats(temp(i)).BoundingBox(4)-1;
    end_y(i)=start_y(i)+stats(temp(i)).BoundingBox(3)-1;
    [X Y]=find(stats(temp(i)).Image);
    test_center(i)=median(Y)+start_y(i)+r_third;
%     if stats(temp(i)).BoundingBox(4)>old
    if test_center(i)>old
        main_center=test_center(i);
        old=test_center(i);
        ii=i;
    end
%     sum_x=sum_x+stats(temp(i)).BoundingBox(4);
end


% i=temp(find(st==max(st)));
% i=(find(st==max(st)));
% i=i(1);
i=ii;
r_cortical=false(size(cortical));
r_cortical(start_x(i):end_x(i),start_y(i)+r_third:end_y(i)+r_third)=r_cortical(start_x(i):end_x(i),start_y(i)+r_third:end_y(i)+r_third)|stats(temp(i)).Image;

sum_x=stats(temp(i)).BoundingBox(4);

for t_i=1:length(temp)
    if i==t_i
%         sum_x=sum_x+stats(temp(t_i)).BoundingBox(4);
        continue
    end
    if test_center(t_i)-main_center<350
        if start_x(t_i)>=start_x(i) & end_x(t_i)<=end_x(i)
            continue
        end
        r_cortical(start_x(t_i):end_x(t_i),start_y(t_i)+r_third:end_y(t_i)+r_third)=r_cortical(start_x(t_i):end_x(t_i),start_y(t_i)+r_third:end_y(t_i)+r_third)|stats(temp(t_i)).Image;
        sum_x=sum_x+stats(temp(t_i)).BoundingBox(4);
    end
end
%
% start_x(i)=round(stats(II(i)).BoundingBox(2));
% start_y(i)=round(stats(II(i)).BoundingBox(1));
% end_x(i)=start_x(i)+stats(II(i)).BoundingBox(4)-1;
% end_y(i)=start_y(i)+stats(II(i)).BoundingBox(3)-1;
% [X Y]=find(stats(II(1)).Image);
% test_center(i)=median(Y)+start_y(i)+r_third;
% main_center=test_center(1);
% 
% %cortical(start_x(i):end_x(i),start_y(i)+r_third:end_y(i)+r_third)=cortical(start_x(i):end_x(i),start_y(i)+r_third:end_y(i)+r_third)+stats(II(i)).FilledImage;
% cortical(start_x(i):end_x(i),start_y(i)+r_third:end_y(i)+r_third)=cortical(start_x(i):end_x(i),start_y(i)+r_third:end_y(i)+r_third)+stats(II(i)).Image;

% start_x(i)=round(stats(i).BoundingBox(2));
% start_y(i)=round(stats(i).BoundingBox(1));
% end_x(i)=start_x(i)+stats(i).BoundingBox(4)-1;
% end_y(i)=start_y(i)+stats(i).BoundingBox(3)-1;
% [X Y]=find(stats(i).Image);
% test_center=median(Y)+start_y(i)+r_third;
% main_center=test_center;

%cortical(start_x(i):end_x(i),start_y(i)+r_third:end_y(i)+r_third)=cortical(start_x(i):end_x(i),start_y(i)+r_third:end_y(i)+r_third)+stats(II(i)).FilledImage;
% cortical(start_x(i):end_x(i),start_y(i)+r_third:end_y(i)+r_third)=cortical(start_x(i):end_x(i),start_y(i)+r_third:end_y(i)+r_third)+stats(i).Image;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





% if sum_x/size(cortical,1)<0.95
%     h_right=h_right+(end_x(i)-start_x(i)+1);
%     right_overlap(start_x(i):end_x(i),start_y(i)+r_third:end_y(i)+r_third)=right_overlap(start_x(i):end_x(i),start_y(i)+r_third:end_y(i)+r_third) + ones(end_x(i)-start_x(i)+1,end_y(i)-start_y(i)+1);
% 
%     t_st_x=start_x(1)+100;
%     t_ed_x=end_x(1)-100;
%     t_st_y=max(r_third,start_y(1)-100+r_third);
%     t_ed_y=min(end_y(1)+100+r_third,size(armt,2));
%     temp_im=armt(1:end,1:end);

%     t_im=zeros(size(temp_im));
%     t_im(:,t_st_y:t_ed_y)=ones(size(t_im,1),t_ed_y-t_st_y+1);
%     t_im(t_st_x:t_ed_x,t_st_y:t_ed_y)=zeros(t_ed_x-t_st_x+1,t_ed_y-t_st_y+1);
%     temp_im=temp_im&t_im;
%     [L n]=bwlabel(temp_im);
%     stats=regionprops(L,'Area','Image','FilledImage','BoundingBox','Orientation'); clear L
% 
%     area=[stats.Area];
%     s_area=sort(area);
%     nn=find(s_area/sum(s_area)>=0.005);
%     clear II
%     for i=1:length(nn)
%         temp=find(area==s_area(length(s_area)-(i-1)));
%         for j=1:length(temp)
%             II(i)=temp(j);
%             i=i+1;
%         end
%     end
% 
%     %j=0;
%     %for i=1:length(nn)
%     %    temp=find(area==s_area(length(s_area)-(i-1)));
%     %    for j=1:length(temp)
%     %        j=j+j;
%     %        II(j+1)=temp(j);
%     %    end
%     %end
% 
%     %clear test_center
%     for i=1:length(nn)
%         start_x(i)=round(stats(II(i)).BoundingBox(2));
%         start_y(i)=round(stats(II(i)).BoundingBox(1));
%         end_x(i)=start_x(i)+stats(II(i)).BoundingBox(4)-1;
%         end_y(i)=start_y(i)+stats(II(i)).BoundingBox(3)-1;
%         [X Y]=find(stats(II(i)).Image);
%     %    test_center(i)=round((end_y(i)-start_y(i))/2)+start_y(i);
%         test_center(i)=median(Y)+start_y(i);
%     end
%     center=abs(test_center-main_center);
%     s_center=sort(center(1:end-1));
%     j=0;
%     while h_right<size(armt,1)*1.1
%         j=j+1;
%         if j>numel(s_center)
%             break
%         end
%         ii=find(center==s_center(j));
%         if length(ii)>1
%             ii=ii(1);
%         end
%     %    if ~(start_x(1)<start_x(ii) & end_x(1)>end_y(ii))
%         if ~(start_x(1)<start_x(ii) & end_x(1)>end_x(ii))
%             %cortical(start_x(ii):end_x(ii),start_y(ii):end_y(ii))=cortical(start_x(ii):end_x(ii),start_y(ii):end_y(ii))+stats(II(ii)).FilledImage;
%             cortical(start_x(ii):end_x(ii),start_y(ii):end_y(ii))=cortical(start_x(ii):end_x(ii),start_y(ii):end_y(ii))+stats(II(ii)).Image;
%             h_right=h_right+(end_x(ii)-start_x(ii)+1);
%             right_overlap(start_x(ii):end_x(ii),start_y(ii):end_y(ii))=right_overlap(start_x(ii):end_x(ii),start_y(ii):end_y(ii))...
%                 + ones(end_x(ii)-start_x(ii)+1,end_y(ii)-start_y(ii)+1);
%         end
%     end
% end






%for i=1:length(nn)
%    start_x(i)=round(stats(II(i)).BoundingBox(2));
%    start_y(i)=round(stats(II(i)).BoundingBox(1));
%    end_x(i)=start_x(i)+stats(II(i)).BoundingBox(4)-1;
%    end_y(i)=start_y(i)+stats(II(i)).BoundingBox(3)-1;
%    [X Y]=find(stats(II(i)).Image);
%%    test_center(i)=round((end_y(i)-start_y(i))/2)+start_y(i);
%    test_center(i)=mean(Y)+start_y(i);
%    if i==1
%        main_center=test_center(1)+start_y(i);
%        cortical(start_x:end_x,r_third+start_y:r_third+end_y)=cortical(start_x:end_x,r_third+start_y:r_third+end_y)+stats(II(i)).FilledImage;
%        h_right=h_right+(end_x-start_x+1);
%        right_overlap(start_x:end_x,r_third+start_y:r_third+end_y)=right_overlap(start_x:end_x,r_third+start_y:r_third+end_y)...
%            + ones(end_x-start_x+1,end_y-start_y+1);
%    elseif h_right>=size(armt,1)
%        break;
%    end
%end
%center=abs(test_center-main_center);
%s_center=sort(center(2:end-1));
%j=0;
%while h_right<size(armt,1)
%    j=j+1;
%    ii=find(center==s_center(j));
%    if ~(start_x(1)<start_x(ii) & end_x(1)>end_y(ii))
%        cortical(start_x(ii):end_x(ii),r_third+start_y(ii):r_third+end_y(ii))=cortical(start_x(ii):end_x(ii),r_third+start_y(ii):r_third+end_y(ii))+stats(II(ii)).FilledImage;
%        h_right=h_right+(end_x(ii)-start_x(ii)+1);
%        right_overlap(start_x(ii):end_x(ii),r_third+start_y(ii):r_third+end_y(ii))=right_overlap(start_x(ii):end_x(ii),r_third+start_y(ii):r_third+end_y(ii))...
%            + ones(end_x(ii)-start_x(ii)+1,end_y(ii)-start_y(ii)+1);
%    end
%end

[L n]=bwlabel(r_cortical(1:end,r_third+1:end));
stats=regionprops(L,'Area','Image','FilledImage','BoundingBox','Orientation'); clear L
clear t_order
for i=1:n
   t_order(i)=stats(i).BoundingBox(2);
end
s_order=sort(t_order);
p=0;
for i=1:n
    temp=find(t_order==s_order(i));
    if length(temp)>1
        p=p+1;
        order(i)=temp(p);
        if p>=length(temp)
            p=0;
        end
    else
        p=0;
        order(i)=temp;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   remove small bones which are overlapped with cortical bone in vertical direction
%   4/16/08
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
temp=round([stats.BoundingBox]);
temp=reshape(temp',4,length(temp)/4);
temp=[temp(2,:);temp(2,:)+temp(4,:)];
temp(2,:)=temp(1,:)+temp(2,:)-1;

s_temp=sort(temp,2);
p=1; j=1;
i=1;
clear s_temp1
while i<size(s_temp,2)+1
    if j>1
        i=i+j-1;
        j=1;
    end
    if i>size(s_temp,2)
        break
    end
    x=find(temp(1,:)==s_temp(1,i));
    for j=1:length(x)
        s_temp1(p)=temp(2,x(j));
        p=p+1;
    end
    i=i+1;
end
s_temp(2,:)=s_temp1;
x=100;
size_s_temp=size(s_temp,2);
while isempty(x)==0
    d_temp=diff(s_temp(2,:),1,2);
    x=find(d_temp<=0);
    if isempty(x)==1
        break
    end
    s_temp(:,x(1)+1:size_s_temp-1)=s_temp(:,x(1)+2:size_s_temp);
    size_s_temp=size_s_temp-1;
    s_temp=s_temp(:,1:size_s_temp);
end

clear order
for i=1:size(s_temp,2)
    x=find(sum(temp-s_temp(:,i)*ones(1,size(temp,2)),1)==0);
    order(i)=min(x);
end

cortical=cortical|r_cortical;
if length(order)>1
    %%%%%%%%%%%%%%%%%%%%%%%%
    %
    %   right broken bone (3/5/08)
    %
    %%%%%%%%%%%%%%%%%%%%%%%%
    for i=1:length(order)-1
        outline1=zeros(size(cortical)); outline2=outline1;
        for j=1:2
            temp_im=stats(order(i-1+j)).Image;
            start_x(j)=round(stats(order(i-1+j)).BoundingBox(2));
            start_y(j)=round(stats(order(i-1+j)).BoundingBox(1));
            end_x(j)=start_x(j)+stats(order(i-1+j)).BoundingBox(4)-1;
            end_y(j)=start_y(j)+stats(order(i-1+j)).BoundingBox(3)-1;
            eval(['outline',num2str(j),'(start_x(j):end_x(j),start_y(j)+r_third:end_y(j)+r_third)=outline',num2str(j),'(start_x(j):end_x(j),start_y(j)+r_third:end_y(j)+r_third)|temp_im;'])
            eval(['temp=imerode(outline',num2str(j),',strel(''disk'',1));'])
            eval(['outline',num2str(j),'=xor(outline',num2str(j),',temp) & (outline',num2str(j),' & ~temp);'])
        end
        [x y d]=find_closest_points(outline1, outline2);
        
        x1=x(1);
        y1=y(1);
        x2=x(2);
        y2=y(2);
        a=(x2-x1)/(y2-y1);
        
        
        if y1==y2
            cortical(min(x1,x2):max(x1,x2),y1-1:y1+1)=true(abs(x2-x1)+1,3);
        elseif x1==x2 & x1>1
            cortical(x1-1:x1+1,min(y1,y2):max(y1,y2))=true(3,abs(y2-y1)+1);
        elseif x1==x2 & x1==1
            cortical(x1:x1+1,min(y1,y2):max(y1,y2))=true(2,abs(y2-y1)+1);
        elseif abs(a)<1
            b=x1-a*y1;
            if y1<y2
                for y=y1:y2
                    if floor(a*y+b)-1<1 | ceil(a*y+b)+1>size(cortical,1)
                        break
                    end
                    cortical(floor(a*y+b)-1:floor(a*y+b)+1,y-1:y+1)=true(3,3);
                    cortical(ceil(a*y+b)-1:ceil(a*y+b)+1,y-1:y+1)=true(3,3);
                end
            else
                for y=y1:-1:y2
                    if floor(a*y+b)-1<1 | ceil(a*y+b)+1>size(cortical,1)
                        break
                    end
                    cortical(floor(a*y+b)-1:floor(a*y+b)+1,y-1:y+1)=true(3,3);
                    cortical(ceil(a*y+b)-1:ceil(a*y+b)+1,y-1:y+1)=true(3,3);
                end
            end
        elseif abs(a)>=1
            a=(y2-y1)/(x2-x1);
            b=y1-a*x1;
            if x1<x2
                for x=x1:x2
                    if floor(a*x+b)-1<1 | ceil(a*x+b)+1>size(cortical,1)
                        break
                    end
                    cortical(x-1:x+1, floor(a*x+b)-1:floor(a*x+b)+1)=true(3,3);
                    cortical(x-1:x+1, ceil(a*x+b)-1:ceil(a*x+b)+1)=true(3,3);
                end
            else
                for x=x1:-1:x2
                    if floor(a*x+b)-1<1 | ceil(a*x+b)+1>size(cortical,1)
                        break
                    end
                    cortical(x-1:x+1, floor(a*x+b)-1:floor(a*x+b)+1)=true(3,3);
                    cortical(x-1:x+1, ceil(a*x+b)-1:ceil(a*x+b)+1)=true(3,3);
                end
            end
        end
    end
end
%for i=1:50
%    cortical=imerode(imdilate(cortical,strel('disk',3)),strel('disk',3));
%    [L n]=bwlabel(cortical);
%    if n==2
%        cortical=imerode(imdilate(cortical,strel('disk',3)),strel('disk',3));
%        [L n]=bwlabel(cortical);
%        break
%    end
%end

%outline=cortical-imerode(cortical,strel('disk',1));
%[L n]=bwlabel(outline);

%cortical=cortical(dis_10x:end,:);
cortical=cortical(start_trabecula:end,:);

cortical(:,third-2:third+2)=false(size(cortical,1),5);
cortical(:,r_third-2:r_third+2)=false(size(cortical,1),5);

[L n]=bwlabel(cortical);
stats=regionprops(L,'Area','Image','FilledImage','BoundingBox','Orientation'); clear L
area=[stats.Area];
s_area=sort(area);

for i=1:2
    II(i)=find(area==s_area(length(s_area)-(i-1)));
end

left_flag_top=0;
right_flag_top=0;
left_flag_bottom=0;
right_flag_bottom=0;
cortical=zeros(size(cortical));
test_left_cortical=cortical;
test_right_cortical=cortical;
clear left_cortical right_cortical

if stats(II(1)).BoundingBox(1)<stats(II(2)).BoundingBox(1)
    st_y=round(stats(II(1)).BoundingBox(1));
    st_x=round(stats(II(1)).BoundingBox(2));
    ed_y=st_y+round(stats(II(1)).BoundingBox(3))-1;
    ed_x=st_x+round(stats(II(1)).BoundingBox(4))-1;
    if st_x~=1
        temp=stats(II(1)).FilledImage;
        left_cortical(1:st_x-1,1:ed_y-st_y+1)=logical(ones(st_x-1,1)*temp(1,:));
        left_flag_top=1;
    end
    if ed_x~=size(cortical,1)
        temp=stats(II(1)).FilledImage;
        left_cortical(ed_x+1:size(cortical,1),1:ed_y-st_y+1)=logical(ones(size(cortical,1)-ed_x,1)*temp(end,:));
        left_flag_bottom=1;
    end
    %left_cortical=stats(II(1)).FilledImage;
    %cortical(st_x:ed_x,st_y:ed_y)=cortical(st_x:ed_x,st_y:ed_y)+left_cortical;
    if left_flag_top==0 & left_flag_bottom==0
        left_cortical=stats(II(1)).Image;
    elseif left_flag_top==1 & left_flag_bottom==0
        left_cortical(st_x:size(cortical,1),1:ed_y-st_y+1)=stats(II(1)).Image;
    elseif left_flag_top==0 & left_flag_bottom==1
       left_cortical(st_x:ed_x,1:ed_y-st_y+1)=stats(II(1)).Image;
    elseif left_flag_top==1 & left_flag_bottom==1
       left_cortical(st_x:ed_x,1:ed_y-st_y+1)=stats(II(1)).Image;
    end
    cortical(:,st_y:ed_y)=cortical(:,st_y:ed_y)|left_cortical;
    test_left_cortical(:,st_y:ed_y)=left_cortical;
    start_left=st_y;
    %end_left=ed_y;   
    
    st_y=round(stats(II(2)).BoundingBox(1));
    st_x=round(stats(II(2)).BoundingBox(2));
    ed_y=st_y+round(stats(II(2)).BoundingBox(3))-1;
    ed_x=st_x+round(stats(II(2)).BoundingBox(4))-1;
    if st_x~=1
        temp=stats(II(2)).FilledImage;
        right_cortical(1:st_x-1,1:ed_y-st_y+1)=logical(ones(st_x-1,1)*temp(1,:));
        right_flag_top=1;
    end
    if ed_x~=size(cortical,1)
        temp=stats(II(2)).FilledImage;
        right_cortical(ed_x+1:size(cortical,1),1:ed_y-st_y+1)=logical(ones(size(cortical,1)-ed_x,1)*temp(end,:));
        right_flag_bottom=1;
    end
    %right_cortical=stats(II(2)).FilledImage;
    %cortical(st_x:ed_x,st_y:ed_y)=cortical(st_x:ed_x,st_y:ed_y)+right_cortical;
    if right_flag_top==0 & right_flag_bottom==0
        right_cortical=stats(II(2)).Image;
    elseif right_flag_top==1 & right_flag_bottom==0
        right_cortical(st_x:size(cortical,1),1:ed_y-st_y+1)=stats(II(2)).Image;
    elseif right_flag_top==0 & right_flag_bottom==1
       right_cortical(st_x:ed_x,1:ed_y-st_y+1)=stats(II(2)).Image;
    elseif right_flag_top==1 & right_flag_bottom==1
       right_cortical(st_x:ed_x,1:ed_y-st_y+1)=stats(II(2)).Image;
    end
    cortical(:,st_y:ed_y)=cortical(:,st_y:ed_y)|right_cortical;
    test_right_cortical(:,st_y:ed_y)=right_cortical;
    %start_right=st_y;
    end_right=ed_y;
else
    st_y=round(stats(II(2)).BoundingBox(1));
    st_x=round(stats(II(2)).BoundingBox(2));
    ed_y=st_y+round(stats(II(2)).BoundingBox(3))-1;
    ed_x=st_x+round(stats(II(2)).BoundingBox(4))-1;
    if st_x~=1
        temp=stats(II(2)).FilledImage;
        left_cortical(1:st_x-1,1:ed_y-st_y+1)=logical(ones(st_x-1,1)*temp(1,:));
        left_flag_top=1;
    end
    if ed_x~=size(cortical,1)
        temp=stats(II(2)).FilledImage;
        left_cortical(ed_x+1:size(cortical,1),1:ed_y-st_y+1)=logical(ones(size(cortical,1)-ed_x,1)*temp(end,:));
        left_flag_bottom=1;
    end
    
    %left_cortical=stats(II(2)).FilledImage;
    %cortical(st_x:ed_x,st_y:ed_y)=cortical(st_x:ed_x,st_y:ed_y)+left_cortical;
    if left_flag_top==0 & left_flag_bottom==0
        left_cortical=stats(II(2)).Image;
    elseif left_flag_top==1 & left_flag_bottom==0
        left_cortical(st_x:size(cortical,1),1:ed_y-st_y+1)=stats(II(2)).Image;
    elseif left_flag_top==0 & left_flag_bottom==1
       left_cortical(st_x:ed_x,1:ed_y-st_y+1)=stats(II(2)).Image;
    elseif left_flag_top==1 & left_flag_bottom==1
       left_cortical(st_x:ed_x,1:ed_y-st_y+1)=stats(II(2)).Image;
    end
    %if left_flag==0
    %    left_cortical(:,1:ed_y-st_y+1)=stats(II(2)).FilledImage;
    %else
    %    left_cortical(:,1:ed_y-st_y+1)=left_cortical(:,1:ed_y-st_y+1)+stats(II(2)).FilledImage;
    %end
    cortical(:,st_y:ed_y)=cortical(:,st_y:ed_y)|left_cortical;
    test_left_cortical(:,st_y:ed_y)=left_cortical;
    start_left=st_y;
    %end_left=ed_y;   
    
    st_y=round(stats(II(1)).BoundingBox(1));
    st_x=round(stats(II(1)).BoundingBox(2));
    ed_y=st_y+round(stats(II(1)).BoundingBox(3))-1;
    ed_x=st_x+round(stats(II(1)).BoundingBox(4))-1;
    if st_x~=1
        temp=stats(II(1)).FilledImage;
        right_cortical(1:st_x-1,1:ed_y-st_y+1)=logical(ones(st_x-1,1)*temp(1,:));
        right_flag_top=1;
    end
    if ed_x~=size(cortical,1)
        temp=stats(II(1)).FilledImage;
        right_cortical(ed_x+1:size(cortical,1),1:ed_y-st_y+1)=logical(ones(size(cortical,1)-ed_x,1)*temp(end,:));
        right_flag_bottom=1;
    end

    %right_cortical=stats(II(1)).FilledImage;
    %cortical(st_x:ed_x,st_y:ed_y)=cortical(st_x:ed_x,st_y:ed_y)+right_cortical;
    if right_flag_top==0 & right_flag_bottom==0
        right_cortical=stats(II(1)).Image;
    elseif right_flag_top==1 & right_flag_bottom==0
        right_cortical(st_x:size(cortical,1),1:ed_y-st_y+1)=stats(II(1)).Image;
    elseif right_flag_top==0 & right_flag_bottom==1
       right_cortical(st_x:ed_x,1:ed_y-st_y+1)=stats(II(1)).Image;
    elseif right_flag_top==1 & right_flag_bottom==1
       right_cortical(st_x:ed_x,1:ed_y-st_y+1)=stats(II(1)).Image;
    end
    %if right_flag==0
    %    right_cortical(:,1:ed_y-st_y+1)=stats(II(1)).FilledImage;
    %else
    %    right_cortical(:,1:ed_y-st_y+1)=right_cortical(:,1:ed_y-st_y+1)+stats(II(1)).FilledImage;
    %end
    cortical(:,st_y:ed_y)=cortical(:,st_y:ed_y)|right_cortical;
    test_right_cortical(:,st_y:ed_y)=right_cortical;
    %start_right=st_y;
    end_right=ed_y;
end


%e=logical(left_cortical)-imerode(left_cortical,strel('disk',1));

%[L n]=bwlabel(cortical);
%stats=regionprops(L,'Area','Image','FilledImage','BoundingBox','Orientation'); clear L
%
%cortical=zeros(size(cortical));
%
%for i=1:n
%    start_x=round(stats(i).BoundingBox(2));
%    start_y=round(stats(i).BoundingBox(1));
%    end_x=start_x+stats(i).BoundingBox(4)-1;
%    end_y=start_y+stats(i).BoundingBox(3)-1;
%    cortical(start_x:end_x,start_y:end_y)=cortical(start_x:end_x,start_y:end_y)+stats(i).FilledImage;
%end



cortical=imdilate(cortical,strel('disk',3)); %cortical=(cortical>0);
cortical(:,1)=false(size(cortical,1),1);
cortical(:,end)=false(size(cortical,1),1);
osteum=xor(cortical,imerode(cortical,strel('disk',1))) & (cortical & ~(imerode(cortical,strel('disk',1))));

[L n]=bwlabel(osteum);
stats=regionprops(L,'Area','Image','FilledImage','BoundingBox','Orientation'); clear L

for i=1:n
    height(i)=stats(i).BoundingBox(4);
end

s_height=sort(height);
II=find(height==s_height(length(s_height)));

%II(1)=find(area==s_area(length(s_area)));
%II(2)=find(area==s_area(length(s_area)-1));
%II(3)=find(area==s_area(length(s_area)-2));
%II(4)=find(area==s_area(length(s_area)-3));

tt=zeros(1,length(II));
for i=1:length(II)
    tt(i)=stats(II(i)).BoundingBox(1);
end
s_tt=sort(tt);

%l_p=logical([1 0 0 0]);
%r_p=logical([0 0 0 1]);
%l_e=logical([0 1 0 0]);
%r_e=logical([0 0 1 0]);

%r_e=logical([1 0 0 0]);
%l_e=logical([0 0 0 1]);
%l_p=logical([0 1 0 0]);
%r_p=logical([0 0 1 0]);

r_e=find(tt==s_tt(length(s_tt)-1));
l_e=find(tt==s_tt(2));
r_p=find(tt==s_tt(length(s_tt)));
l_p=find(tt==s_tt(1));

left_periosteum=stats(II(l_p)).Image;
right_periosteum=stats(II(r_p)).Image;
end_left=start_left+size(left_periosteum,2)-1;
start_right=end_right-size(right_periosteum,2)+1;





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   replaced with routine to find new endosteum
%
%   Jan. 8, 2009
%
%
%left_endosteum=stats(II(l_e)).Image;
%right_endosteum=stats(II(r_e)).Image;
%%left_periosteum=stats(II(l_p)).Image;
%%right_periosteum=stats(II(r_p)).Image;
%%left_endosteum=stats(II(l_e)).Image;
%%right_endosteum=stats(II(r_e)).Image;
%
%left_endosteum_x=round(stats(II(l_e)).BoundingBox(2));
%end_left_endosteum_x=left_endosteum_x+round(stats(II(l_e)).BoundingBox(4))-1;
%left_endosteum_y=round(stats(II(l_e)).BoundingBox(1));
%end_left_endosteum_y=left_endosteum_y+round(stats(II(l_e)).BoundingBox(3))-1;
%endosteum(left_endosteum_x:end_left_endosteum_x,left_endosteum_y:end_left_endosteum_y)...
%    =endosteum(left_endosteum_x:end_left_endosteum_x,left_endosteum_y:end_left_endosteum_y)+stats(II(l_e)).Image;
%
%right_endosteum_x=round(stats(II(r_e)).BoundingBox(2));
%end_right_endosteum_x=right_endosteum_x+round(stats(II(r_e)).BoundingBox(4))-1;
%right_endosteum_y=round(stats(II(r_e)).BoundingBox(1));
%end_right_endosteum_y=right_endosteum_y+round(stats(II(r_e)).BoundingBox(3))-1;
%endosteum(right_endosteum_x:end_right_endosteum_x,right_endosteum_y:end_right_endosteum_y)...
%    =endosteum(right_endosteum_x:end_right_endosteum_x,right_endosteum_y:end_right_endosteum_y)+stats(II(r_e)).Image;
%
%
%
test_count_right_origin=length(find(test_right_cortical(:,start_right:end_right) & right_periosteum));
new_right_cortical=false(size(cortical));
new_right_cortical(:,start_right:end_right)=right_periosteum;
start_right1=start_right;
end_right1=end_right;
ratio_count_right=1;
count_right=1;
cc_r=0;

thresh_cut_protrude_from_cortical_right=0.7;
thresh_cut_protrude_from_cortical_left=0.7;

while ratio_count_right(count_right)>thresh_cut_protrude_from_cortical_right
    count_right=count_right+1;
    start_right1=start_right1-1;
    end_right1=end_right1-1;
    cc_r=cc_r+1;
    
    test_count_right(count_right)=length(find(test_right_cortical(:,start_right1:end_right1) & right_periosteum));
%     ratio_count_right(count_right)=test_count_right(count_right)/test_cou
%     nt_right_origin; %May 19 2013
    ratio_count_right(count_right)=test_count_right(count_right)/max(test_count_right);
    new_right_cortical(:,start_right1:end_right1)=new_right_cortical(:,start_right1:end_right1)|(test_right_cortical(:,start_right1:end_right1)&right_periosteum);

    if cc_r>500*rat_mouse
        break;
    end
end
new_right_cortical(:,1)=false(size(new_right_cortical,1),1);
new_right_cortical(:,end)=false(size(new_right_cortical,1),1);
new_right_cortical_boundary=xor(new_right_cortical,imerode(new_right_cortical,strel('disk',1))) & (new_right_cortical & ~imerode(new_right_cortical,strel('disk',1)));

[L n]=bwlabel(new_right_cortical_boundary);
stats_new_r=regionprops(L,'Area','Image','FilledImage','BoundingBox','Orientation'); clear L
if cc_r>500*rat_mouse
    right_endosteum=stats(II(r_e)).Image;
    
    right_endosteum_x=round(stats(II(r_e)).BoundingBox(2));
    end_right_endosteum_x=right_endosteum_x+round(stats(II(r_e)).BoundingBox(4))-1;
    right_endosteum_y=round(stats(II(r_e)).BoundingBox(1));
    end_right_endosteum_y=right_endosteum_y+round(stats(II(r_e)).BoundingBox(3))-1;
    right_endosteum=stats(II(r_e)).Image;
    endosteum(right_endosteum_x:end_right_endosteum_x,right_endosteum_y:end_right_endosteum_y)...
        =endosteum(right_endosteum_x:end_right_endosteum_x,right_endosteum_y:end_right_endosteum_y)|stats(II(r_e)).Image;
else
    for length_n_r=1:n
        if stats_new_r(length_n_r).BoundingBox(4)>size(new_right_cortical_boundary,1)*0.8
            right_endosteum_x=round(stats_new_r(length_n_r).BoundingBox(2));
            end_right_endosteum_x=right_endosteum_x+round(stats_new_r(length_n_r).BoundingBox(4))-1;
            right_endosteum_y=round(stats_new_r(length_n_r).BoundingBox(1));
            end_right_endosteum_y=right_endosteum_y+round(stats_new_r(length_n_r).BoundingBox(3))-1;
            right_endosteum=stats_new_r(length_n_r).Image;
            endosteum(right_endosteum_x:end_right_endosteum_x,right_endosteum_y:end_right_endosteum_y)...
                =endosteum(right_endosteum_x:end_right_endosteum_x,right_endosteum_y:end_right_endosteum_y)|stats_new_r(length_n_r).Image;
            break;
        end
    end
end

test_count_left_origin=length(find(test_left_cortical(:,start_left:end_left)&left_periosteum));
new_left_cortical=false(size(cortical));
new_left_cortical(:,start_left:end_left)=left_periosteum;
start_left1=start_left;
end_left1=end_left;
ratio_count_left=1;
count_left=1;
cc=0;
while ratio_count_left(count_left)>thresh_cut_protrude_from_cortical_left
    count_left=count_left+1;
    start_left1=start_left1+1;
    end_left1=end_left1+1;
    cc=cc+1;
    
    test_count_left(count_left)=length(find(test_left_cortical(:,start_left1:end_left1)&left_periosteum));
%     ratio_count_left(count_left)=test_count_left(count_left)/test_count_left_origin; %May 19 2013
    ratio_count_left(count_left)=test_count_left(count_left)/max(test_count_left);
    new_left_cortical(:,start_left1:end_left1)=new_left_cortical(:,start_left1:end_left1)|(test_left_cortical(:,start_left1:end_left1)&left_periosteum);
    
    if cc>500*rat_mouse
        break;
    end
end
new_left_cortical(:,1)=false(size(new_left_cortical,1),1);
new_left_cortical(:,end)=false(size(new_left_cortical,1),1);
new_left_cortical_boundary=xor(new_left_cortical,imerode(new_left_cortical,strel('disk',1))) & (new_left_cortical & ~imerode(new_left_cortical,strel('disk',1)));

[L n]=bwlabel(new_left_cortical_boundary);
stats_new=regionprops(L,'Area','Image','FilledImage','BoundingBox','Orientation'); clear L

if cc>500*rat_mouse
    left_endosteum=stats(II(l_e)).Image;
    
    left_endosteum_x=round(stats(II(l_e)).BoundingBox(2));
    end_left_endosteum_x=left_endosteum_x+round(stats(II(l_e)).BoundingBox(4))-1;
    left_endosteum_y=round(stats(II(l_e)).BoundingBox(1));
    end_left_endosteum_y=left_endosteum_y+round(stats(II(l_e)).BoundingBox(3))-1;
    endosteum(left_endosteum_x:end_left_endosteum_x,left_endosteum_y:end_left_endosteum_y)...
        =endosteum(left_endosteum_x:end_left_endosteum_x,left_endosteum_y:end_left_endosteum_y)|stats(II(l_e)).Image;
else
    for length_n=n:-1:1
        if stats_new(length_n).BoundingBox(4)>size(new_right_cortical_boundary,1)*0.8
            left_endosteum_x=round(stats_new(length_n).BoundingBox(2));
            end_left_endosteum_x=left_endosteum_x+round(stats_new(length_n).BoundingBox(4))-1;
            left_endosteum_y=round(stats_new(length_n).BoundingBox(1));
            end_left_endosteum_y=left_endosteum_y+round(stats_new(length_n).BoundingBox(3))-1;
            left_endosteum=stats_new(length_n).Image;
            endosteum(left_endosteum_x:end_left_endosteum_x,left_endosteum_y:end_left_endosteum_y)...
                =endosteum(left_endosteum_x:end_left_endosteum_x,left_endosteum_y:end_left_endosteum_y)|stats_new(length_n).Image;
            break;
        end
    end
end
%
%   Jan. 8, 2009
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





[x y]=find(left_endosteum'==1);

mm=zeros(1,max(y));
m=zeros(1,length(I));
for i=1:max(y)
    I=find(y==i);
    for j=1:length(I)
        m(j)=x(I(j));
    end
    %mm(i)=mean(m(1:j));
    mm(i)=max(m(1:j));
end

%st=zeros(1,floor((size(left_endosteum,1)-dis_10x)/lr_sp_tb));
%for i=1:floor((size(left_endosteum,1)-dis_10x)/lr_sp_tb)
%%    st(i)=mean(mm(lr_sp_tb*(i-1)+1:lr_sp_tb*i))+left_endosteum_y+lr_sp_tb;    
%    st(i)=max(mm(lr_sp_tb*(i-1)+dis_10x+1:lr_sp_tb*i+dis_10x))+left_endosteum_y+lr_sp_tb;    
%end
%%st(i+1)=mean(mm(lr_sp_tb*(i)+1:end))+left_endosteum_y+lr_sp_tb;    
%st(i+1)=max(mm(lr_sp_tb*(i)+1+dis_10x:end))+left_endosteum_y+lr_sp_tb;    
%%st=round(s_tt(2))+stats(II(l_e)).BoundingBox(3)+lr_sp_tb;
%
%[x y]=find(right_endosteum'==1);
%
%mm=zeros(1,max(y));
%m=zeros(1,length(I));
%for i=1:max(y)
%    I=find(y==i);
%    for j=1:length(I)
%        m(j)=x(I(j));
%    end
%    %mm(i)=mean(m(1:j));
%    mm(i)=min(m(1:j));
%end
%
%et=zeros(1,floor((size(left_endosteum,1)-dis_10x)/lr_sp_tb));
%for i=1:floor((size(left_endosteum,1)-dis_10x)/lr_sp_tb)
%    et(i)=end_right_endosteum_y-(stats(II(r_e)).BoundingBox(3)-mean(mm(lr_sp_tb*(i-1)+dis_10x+1:lr_sp_tb*i+dis_10x)))-lr_sp_tb;
%end
%et(i+1)=end_right_endosteum_y-(stats(II(r_e)).BoundingBox(3)-mean(mm(lr_sp_tb*(i)+dis_10x+1:end)))-lr_sp_tb;
%%et=round(s_tt(3))-lr_sp_tb;%end_right_endosteum_x-size(left_endosteum,2);

lr_sp_tb1=dis_10x/2;

st=zeros(1,floor((size(left_endosteum,1))/lr_sp_tb));
for i=1:floor((size(left_endosteum,1))/lr_sp_tb)
%    st(i)=mean(mm(lr_sp_tb*(i-1)+1:lr_sp_tb*i))+left_endosteum_y+lr_sp_tb;    
    st(i)=max(mm(lr_sp_tb*(i-1)+1:lr_sp_tb*i))+left_endosteum_y+lr_sp_tb1;    
end
%st(i+1)=mean(mm(lr_sp_tb*(i)+1:end))+left_endosteum_y+lr_sp_tb;    

if floor((size(left_endosteum,1))/lr_sp_tb)<(size(left_endosteum,1))/lr_sp_tb
    st(i+1)=max(mm(lr_sp_tb*(i)+1:end))+left_endosteum_y+lr_sp_tb1;    
end
%st=round(s_tt(2))+stats(II(l_e)).BoundingBox(3)+lr_sp_tb;

[x y]=find(right_endosteum'==1);

mm=zeros(1,max(y));
m=zeros(1,length(I));
for i=1:max(y)
    I=find(y==i);
    for j=1:length(I)
        m(j)=x(I(j));
    end
    %mm(i)=mean(m(1:j));
    mm(i)=min(m(1:j));
end

et=zeros(1,floor((size(left_endosteum,1))/lr_sp_tb));
et_count=floor((size(left_endosteum,1))/lr_sp_tb);
for i=1:et_count
%    et(i)=end_right_endosteum_y-(stats(II(r_e)).BoundingBox(3)-mean(mm(lr_sp_tb*(i-1)+1:lr_sp_tb*i)))-lr_sp_tb;
    if cc_r>500*rat_mouse
        et(i)=end_right_endosteum_y-(stats(II(r_e)).BoundingBox(3)-min(mm(lr_sp_tb*(i-1)+1:lr_sp_tb*i)))-lr_sp_tb1;
    else
        et(i)=end_right_endosteum_y-(round(stats_new_r(length_n_r).BoundingBox(3))-min(mm(lr_sp_tb*(i-1)+1:lr_sp_tb*i)))-lr_sp_tb1;
    end
end
%et(i+1)=end_right_endosteum_y-(stats(II(r_e)).BoundingBox(3)-mean(mm(lr_sp_tb*(i)+1:end)))-lr_sp_tb;



if cc_r>500*rat_mouse & et_count<(size(left_endosteum,1))/lr_sp_tb
    et(i+1)=end_right_endosteum_y-(stats(II(r_e)).BoundingBox(3)-min(mm(lr_sp_tb*(i)+1:end)))-lr_sp_tb1;
elseif et_count<(size(left_endosteum,1))/lr_sp_tb
    et(i+1)=end_right_endosteum_y-(round(stats_new_r(length_n_r).BoundingBox(3))-min(mm(lr_sp_tb*(i)+1:end)))-lr_sp_tb1;
end
%et=round(s_tt(3))-lr_sp_tb;%end_right_endosteum_x-size(left_endosteum,2);

%trabecula=zeros(size(armt));
trabecula=false(size(result));

% 11/30/07 for i=2:floor(size(right_endosteum,1)/margi)
% 12/19/07 for i=2:length(find(cumsum(et-st)-margi*8<200))

%h=250*lr_sp_tb/dis_10x;
% if mouse==1
    band=lr_sp_tb;
% else
%     band=lr_sp_tb/2;
% end

% for i=1:length(et)
for i=1:min(11,length(et))          % limit the depth to 11xlr_sp_tb, July 14, 2015
%for i=1:length(find((cumsum(et-st)*lr_sp_tb)*(100/dis_10x*2)^2-2.115*10^6<100))
    if band*i>size(armt,1)
        break
    end
    if i>2 & size(armt,1)-(start_trabecula+band*(i-1)+1)<band
        trabecula(start_trabecula+band*(i-1)+1:end,st(i):et(i))=true(size(armt,1)-(start_trabecula+band*(i-1)),et(i)-st(i)+1);
        break
    end
    trabecula(start_trabecula+band*(i-1)+1:start_trabecula+band*i,st(i):et(i))=true(band,et(i)-st(i)+1);
end
% 11/30/07 trabecula(margi*(i)+1:end,st(i+1):et(i+1))=ones(size(right_endosteum,1)-margi*i,et(i+1)-st(i+1)+1);
%tissue_volume=sum(start_trabecula*(et-st));

tissue_volume=length(find(trabecula))*dis_10x_for_400micron_per_312pixels^2;%*(400/dis_10x)^2;  % 3/3
%tr_frame=uint8(imrotate(trabecula-imerode(trabecula,strel('disk',3)),-rot_angle));

tr=false(size(trabecula));
total=0;
row=1;

% if mouse
    roi_size=round(2.1125*10^6*(rat_mouse^2)/(dis_10x_for_400micron_per_312pixels^2));
% else
%     roi_size=2.1125*10^6/4;
% end

while (total<roi_size & row<size(trabecula,1)+1)
%while (total<1125600/((area_roi_tb/(dis_10x*2))^2) & row<size(trabecula,1)+1)
%     total=total+((400/dis_10x*rat_mouse)^2)*length(find(trabecula(row,:)));
    total=total+length(find(trabecula(row,:)));
    tr(row,:)=trabecula(row,:);
    row=row+1;
end
% figure;imshow(tr)
% if row<size(trabecula,1)
%     trabecula1=tr;
%     trabecula2=trabecula-trabecula1;
%     tr_frame=uint8(imrotate(trabecula1-imerode(trabecula1,strel('disk',3)),-rot_angle));
%     tr_frame2=uint8(imrotate(trabecula2-imerode(trabecula2,strel('disk',3)),-rot_angle));
%     no_tr_region=2;
%     tissue_volume1=length(find(trabecula1));         % Dec. 08, 2009
%     tissue_volume2=length(find(trabecula2));         % Dec. 08, 2009
% else
%     tr_frame=uint8(imrotate(trabecula-imerode(trabecula,strel('disk',3)),-rot_angle));
%     tissue_volume=length(find(trabecula));
%     tissue_volume1=tissue_volume;
%     no_tr_region=1;
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   March 15, 2016
% 
%   Remove thick trabeculae attached to growth plate
%
ima_temp=tr&result;
[L n]=bwlabel(ima_temp);
stat_ROI=regionprops(L, 'BoundingBox', 'Image', 'Area'); clear L
ima=false(size(ima_temp));
for i=1:n
    if round(stat_ROI(i).BoundingBox(2))~=313
        continue
    else
        ima_temp=stat_ROI(i).Image;
        ima_temp1=ima_temp(1,:);
        ima_temp1(1)=0;         ima_temp1(end)=0;
        cand=find(diff(ima_temp1,1, 2));
        if ~isempty(cand)
%             if ima_temp(1,1)==1
%                 cand=[0, cand];
%             end
            cand=diff(cand);
            max_d=max(cand(find(mod([1:length(cand)],2))));
        else
            max_d=stat_ROI(i).BoundingBox(4);
        end
        if max_d>100
            st_x=round(stat_ROI(i).BoundingBox(2));
            st_y=round(stat_ROI(i).BoundingBox(1));
            ed_x=st_x+stat_ROI(i).BoundingBox(4)-1;
            ed_y=st_y+stat_ROI(i).BoundingBox(3)-1;
            ima(st_x:ed_x,st_y:ed_y)=ima(st_x:ed_x,st_y:ed_y)|stat_ROI(i).Image;
        end
    end
end
dif=diff(ima,1, 2);
for i=1:size(dif,1)
    temp=find(dif(i,:));
    if ~isempty(temp)
        if dif(i,1)==1
            temp=[0 temp];
        end
        cand=diff(temp);
        max_width(i)=max(cand(find(mod([1:length(cand)],2))));
    else
        max_width(i)=0;
    end
end
% cutoff_ROI_row=find(max_width>200,1,'last')-1;   % from row 1 to row cutoff_ROI_row
temp=find(max_width>200);
if isempty(temp)
    cutoff_ROI_row=0;
else
    ind=find(diff(temp)<100,1,'last');
    if isempty(ind)
        cutoff_ROI_row=find(max_width>200,1,'last');
    else
        cutoff_ROI_row=temp(ind+1);   % from row 1 to row cutoff_ROI_row
    end
end
%
%   March 15, 2016
% 
%   Remove thick trabeculae attached to growth plate
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if row<size(trabecula,1)
    trabecula1=tr;
    
    trabecula_for_cortex=tr;        %   March 15, 2016
    
    trabecula2=xor(trabecula,trabecula1)&(trabecula & ~trabecula1);
    temp_frame=imerode(trabecula1,strel('disk',3));
    tr_frame=xor(trabecula1,temp_frame) & (trabecula1 & ~temp_frame);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   March 15, 2016
% 
%   Remove thick trabeculae attached to growth plate
%    
    if cutoff_ROI_row>0
        trabecula1(1:cutoff_ROI_row,:)=false(cutoff_ROI_row,size(trabecula,2));
    end
    temp_frame=imerode(trabecula1,strel('disk',3));
    tr_frame=tr_frame | xor(trabecula1,temp_frame) & (trabecula1 & ~temp_frame);
%
%   March 15, 2016
% 
%   Remove thick trabeculae attached to growth plate
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    temp_frame=imerode(trabecula2,strel('disk',3));
    tr_frame2=xor(trabecula2,temp_frame) & (trabecula2 & ~temp_frame);
%     no_tr_region=2;
    no_tr_region=1;
    tissue_volume1=length(find(trabecula1))*dis_10x_for_400micron_per_312pixels^2;         % Dec. 08, 2009
    tissue_volume2=length(find(trabecula2))*dis_10x_for_400micron_per_312pixels^2;         % Dec. 08, 2009
else
    trabecula1=tr;
    
    trabecula_for_cortex=tr;        %   March 15, 2016
    
    temp_frame=imerode(trabecula1,strel('disk',3));
    tr_frame=xor(trabecula1,temp_frame) & (trabecula1 & ~temp_frame);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   March 15, 2016
% 
%   Remove thick trabeculae attached to growth plate
%   
    if cutoff_ROI_row>0
        trabecula1(1:cutoff_ROI_row,:)=false(cutoff_ROI_row,size(trabecula,2));
    end
    temp_frame=imerode(trabecula1,strel('disk',3));
    tr_frame=tr_frame | xor(trabecula1,temp_frame) & (trabecula1 & ~temp_frame);
%
%   March 15, 2016
% 
%   Remove thick trabeculae attached to growth plate
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tissue_volume=length(find(trabecula1))*dis_10x_for_400micron_per_312pixels^2;
    tissue_volume1=tissue_volume;
    no_tr_region=1;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Analysis of Signal Intensity
%
%   November 4, 2014
%
c=phr1{b_t,:};
c1=[c,'L',num2str(s),'_s',nnn];
eval(['red=imrotate(imread(''',direc,c1,'_1_shift3.jpg'',''jpg''),rot_angle);']);
red=red(:,:,1);
red_beads=find_circles(red>180, 30, 400)>0;
[x, y, I]=find(immultiply(uint8(red_beads>0),uint8(red)));
I=double(I);
if ~isempty(I)
    analysis.intensity.Red_beads=median(I);
    analysis.intensity.maxRed_beads=max(I);
else
    analysis.intensity.Red_beads=NaN;
    analysis.intensity.maxRed_beads=NaN;
end
eval(['check_file=dir(''',direc,c1,'_1_shift3_BR.jpg'');']);
if isempty(check_file)
    red1=red;
else
    eval(['red1=imrotate(imread(''',direc,c1,'_1_shift3_BR.jpg'',''jpg''),rot_angle);']);
    if size(red1)~=size(red)
        red1=red;
    end
end
red=imrotate(red,ang);
red=imresize(red,1/ratio);
red=red(start_point:end_point,:,:);

[x, y, I]=find(immultiply(uint8(r&trabecula1),uint8(red)));
I=double(I);
if ~isempty(I)
    analysis.intensity.AC=median(I);
    analysis.intensity.maxAC=max(I);
else
    analysis.intensity.AC=NaN;
    analysis.intensity.maxAC=NaN;
end

if ~isnan(analysis.intensity.AC) & ~isnan(analysis.intensity.Red_beads)
    analysis.intensity.AC_Red_beads=analysis.intensity.AC/analysis.intensity.Red_beads;
else
    analysis.intensity.AC_Red_beads=NaN;
end

% eval(['red1=imrotate(imread(''',direc,c1,'_1_shift2_BR.jpg'',''jpg''),rot_angle);']);
red1=red1-medfilt2(red1, [51 51]);   % image resize scale is 0.5
red1=imrotate(red1,ang);
red1=imresize(red1,1/ratio);
red1=red1(start_point:end_point,:,:);

[x, y, I]=find(immultiply(uint8(r&trabecula1),uint8(red1)));
I=double(I);
if ~isempty(I)
    analysis.intensity.AC_BR=median(I);
    analysis.intensity.maxAC_BR=max(I);
else
    analysis.intensity.AC_BR=NaN;
    analysis.intensity.maxAC_BR=NaN;
end

if ~isnan(analysis.intensity.AC_BR) & ~isnan(analysis.intensity.Red_beads)
    analysis.intensity.AC_BR_Red_beads=analysis.intensity.AC_BR/analysis.intensity.Red_beads;
else
    analysis.intensity.AC_BR_Red_beads=NaN;
end
clear red red1

eval(['green=imrotate(imread(''',direc,c1,'_0_shift3.jpg'',''jpg''),rot_angle);']);

if size(green,3)==3
    green=green(:,:,2);
end
green_beads=find_circles(green>180, 30, 400)>0;
[x, y, I]=find(immultiply(uint8(green_beads>0),uint8(green)));
I=double(I);
if ~isempty(I)
    analysis.intensity.Green_beads=median(I);
    analysis.intensity.maxGreen_beads=max(I);
else
    analysis.intensity.Green_beads=NaN;
    analysis.intensity.maxGreen_beads=NaN;
end
green=imrotate(green,ang);
green=imresize(green,1/ratio);
green=green(start_point:end_point,:,:);
[x, y, I]=find(immultiply(uint8(g&trabecula1),uint8(green)));
I=double(I);
if ~isempty(I)
    analysis.intensity.Calcein=median(I);
    analysis.intensity.maxCalcein=max(I);
else
    analysis.intensity.Calcein=NaN;
    analysis.intensity.maxCalcein=NaN;
end

if ~isnan(analysis.intensity.Calcein) & ~isnan(analysis.intensity.Green_beads)
    analysis.intensity.Calcein_Green_beads=analysis.intensity.Calcein/analysis.intensity.Green_beads;
else
    analysis.intensity.Calcein_Green_beads=NaN;
end
clear green
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   AP and TRAP Intensity Analysis
%
%   March 24, 2016
%
%   Should be moved to bottom of AC, Calcein intensity analysis
%
% load([direct_for_old,c1,'_trabecula_bone1_tr1.mat'])
% load([direct_for_old,'data\',c1,'_analysis_n_tr1.mat'])
% load([direct_for_old,'data\',c1,'_analysis_tr1.mat'])

eval(['AP_exist_test=dir(''',direc,c1,'_5_shift3.jpg'');']);
if isempty(AP_exist_test)
    analysis.intensity.AP=NaN;
    analysis.intensity.maxAP=NaN;
    analysis_n.intensity.AP=NaN;
    analysis_n.intensity.maxAP=NaN;
    analysis.intensity.AP_Background_mean=NaN;
    analysis.intensity.AP_Background_median=NaN;
    analysis.intensity.AP_Background_max=NaN;
    analysis.intensity.AP_Background_min=NaN;
	analysis_n.intensity.AP_Background_mean=NaN;
    analysis_n.intensity.AP_Background_median=NaN;
    analysis_n.intensity.AP_Background_max=NaN;
    analysis_n.intensity.AP_Background_min=NaN;
else
    eval(['AP_intensity=imrotate(imread(''',direc,c1,'_5_shift3.jpg'',''jpg''),rot_angle);']);
    AP_intensity=AP_intensity(:,:,1);
    AP_intensity=imrotate(AP_intensity,ang);
    AP_intensity=imresize(AP_intensity,1/ratio);
    AP_intensity=AP_intensity(start_point:end_point,:,:);
    [x, y, I]=find(immultiply(uint8(AP&trabecula1),uint8(AP_intensity)));
    I=double(I);
    if ~isempty(I)
        analysis.intensity.AP=median(I);
        analysis.intensity.maxAP=max(I);
        analysis_n.intensity.AP=median(I);
        analysis_n.intensity.maxAP=max(I);
    else
        analysis.intensity.AP=NaN;
        analysis.intensity.maxAP=NaN;
        analysis_n.intensity.AP=NaN;
        analysis_n.intensity.maxAP=NaN;
    end
    [x, y, I]=find(immultiply(uint8((1-AP)&trabecula1),uint8(AP_intensity)));
    I=double(I);
    if ~isempty(I)
		analysis.intensity.AP_Background_mean=mean(I);
		analysis.intensity.AP_Background_median=median(I);
		analysis.intensity.AP_Background_max=max(I);
		analysis.intensity.AP_Background_min=min(I);
		analysis_n.intensity.AP_Background_mean=mean(I);
		analysis_n.intensity.AP_Background_median=median(I);
		analysis_n.intensity.AP_Background_max=max(I);
		analysis_n.intensity.AP_Background_min=min(I);
    else
		analysis.intensity.AP_Background_mean=NaN;
		analysis.intensity.AP_Background_median=NaN;
		analysis.intensity.AP_Background_max=NaN;
		analysis.intensity.AP_Background_min=NaN;
		analysis_n.intensity.AP_Background_mean=NaN;
		analysis_n.intensity.AP_Background_median=NaN;
		analysis_n.intensity.AP_Background_max=NaN;
		analysis_n.intensity.AP_Background_min=NaN;
    end
end

eval(['TRAP_exist_test=dir(''',direc,c1,'_3_shift3.jpg'');']);
if isempty(TRAP_exist_test)
    analysis.intensity.TRAP=NaN;
    analysis.intensity.maxTRAP=NaN;
    analysis_n.intensity.TRAP=NaN;
    analysis_n.intensity.maxTRAP=NaN;
    analysis.intensity.TRAP_Background_mean=NaN;
    analysis.intensity.TRAP_Background_median=NaN;
    analysis.intensity.TRAP_Background_max=NaN;
    analysis.intensity.TRAP_Background_min=NaN;
	analysis_n.intensity.TRAP_Background_mean=NaN;
    analysis_n.intensity.TRAP_Background_median=NaN;
    analysis_n.intensity.TRAP_Background_max=NaN;
    analysis_n.intensity.TRAP_Background_min=NaN;
else
    eval(['TRAP_intensity=imrotate(imread(''',direc,c1,'_3_shift3.jpg'',''jpg''),rot_angle);']);
    TRAP_intensity=TRAP_intensity(:,:,1);
    TRAP_intensity=imrotate(TRAP_intensity,ang);
    TRAP_intensity=imresize(TRAP_intensity,1/ratio);
    TRAP_intensity=TRAP_intensity(start_point:end_point,:,:);
    [x, y, I]=find(immultiply(uint8(TRAP&trabecula1),uint8(TRAP_intensity)));
    I=double(I);
    if ~isempty(I)
        analysis.intensity.TRAP=median(I);
        analysis.intensity.maxTRAP=max(I);
        analysis_n.intensity.TRAP=median(I);
        analysis_n.intensity.maxTRAP=max(I);
    else
        analysis.intensity.TRAP=NaN;
        analysis.intensity.maxTRAP=NaN;
        analysis_n.intensity.TRAP=NaN;
        analysis_n.intensity.maxTRAP=NaN;
    end
    [x, y, I]=find(immultiply(uint8((1-TRAP)&trabecula1),uint8(TRAP_intensity)));
    I=double(I);
    if ~isempty(I)
		analysis.intensity.TRAP_Background_mean=mean(I);
		analysis.intensity.TRAP_Background_median=median(I);
		analysis.intensity.TRAP_Background_max=max(I);
		analysis.intensity.TRAP_Background_min=min(I);
		analysis_n.intensity.TRAP_Background_mean=mean(I);
		analysis_n.intensity.TRAP_Background_median=median(I);
		analysis_n.intensity.TRAP_Background_max=max(I);
		analysis_n.intensity.TRAP_Background_min=min(I);
    else
		analysis.intensity.TRAP_Background_mean=NaN;
		analysis.intensity.TRAP_Background_median=NaN;
		analysis.intensity.TRAP_Background_max=NaN;
		analysis.intensity.TRAP_Background_min=NaN;
		analysis_n.intensity.TRAP_Background_mean=NaN;
		analysis_n.intensity.TRAP_Background_median=NaN;
		analysis_n.intensity.TRAP_Background_max=NaN;
		analysis_n.intensity.TRAP_Background_min=NaN;
    end
end
% analysis.intensity
% eval(['save ''',direct_for_old,'data\',c1,'_analysis_n_tr1.mat'' analysis_n']);
% eval(['save ''',direct_for_old,'data\',c1,'_analysis_tr1.mat'' analysis']);
% 
% return
%
%   AP and TRAP Intensity Analysis
%
%   March 24, 2016
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%










%if b_t==1
    c=phr1{b_t,:};
    %c1=[c,'\',c,'_',nnn,'FL\',c,'_',nnn,'FL.jpg_Files\',c,'_',nnn,'FL_c'];
    %c1=[c,nnn,'_c'];
    
    %if samples_per_bone==1
    %    c1=[c,nnn,'FL_c'];
    %else
    if strcmp(exp_type,'V')
        c1=[c,nnn,'_h',exp_type,'_s',num2str(s)];
    elseif strcmp(exp_type,'F')
        c1=[c,'L',num2str(s),'_s',nnn];
    end
%end
    
    
    %eval(['im=imrotate(imread(''',direc,c1,mb_com,''',''jpg''),rot_angle);']);
%    eval(['im=imrotate(imread(''',direc,c1,'_pseudo1.jpg'',''jpg''),rot_angle);']);

    eval(['im=imrotate(imread(''',direc,c1,'_shift3_NoDAPI.jpg'',''jpg''),rot_angle);']);
    im=imrotate(im,ang);
    
%elseif b_t==2
%    c=phr(b_t,:);
%    c1=[c,'\',c,'_',nnn,'FL\',c,'_',nnn,'FL.jpg_Files\',c,'_',nnn,'FL_c'];
%    eval(['im=imrotate(imread(''',direc,c1,mb_com,''',''jpg''),rot_angle);']);
%elseif b_t==3
%    c=phr(b_t,:);
%    c1=[c,'\',c,'_',nnn,'FL\',c,'_',nnn,'FL.jpg_Files\',c,'_',nnn,'FL_c'];
%    eval(['im=imrotate(imread(''',direc,c1,mb_com,''',''jpg''),rot_angle);']);
%elseif b_t==4
%    c=phr(b_t,:);
%    c1=[c,'\',c,'_',nnn,'FL\',c,'_',nnn,'FL.jpg_Files\',c,'_',nnn,'FL_c'];
%    eval(['im=imrotate(imread(''',direc,c1,mb_com,''',''jpg''),rot_angle);']);
%end






%if exp_type=='o'
%    eval(['im=imread(''',direc,b_t,b_t,b_n,exp_t,b_t1,b_n,common,'.jpg_Files',b_t1,b_n,common,'_c.jpg'',''jpg'');'])
%else
%    eval(['im=imread(''',direc,b_t,b_t,b_n,b_t1,b_n,common,'.jpg_Files',b_t1,b_n,common,'_c.jpg'',''jpg'');'])
%end

%imwrite(im,'C:\Users\shhong\Desktop\CJake\Analysis s2 3\JW3_WT\JW3_WT_5FR_v2\JW3_WT_5FR_fix_c.jpg_Files\JW3_WT_5FR_fix_c_pseudo.jpg','jpg','quality',100);
im=imresize(im,1/ratio);
%im=im(:,1:end-down_margin,:);

%im=im(1:end-down_margin,:,:);

%im=im(start_trabecula:end_point,:,:);



im3=im(start_point:end_point,:,:);
im4=im3;
im1=im;%imrotate(im3,-rot_angle);
%im1=imrotate(im,rot_angle);
im2=im1;

im1(start_point:end_point,:,1)=imadd(im1(start_point:end_point,:,1), uint8(imdilate(tr_frame,strel('disk',3))*128));
im1(start_point:end_point,:,2)=imadd(im1(start_point:end_point,:,2), uint8(imdilate(tr_frame,strel('disk',3))*128));
im1(start_point:end_point,:,3)=imadd(im1(start_point:end_point,:,3), uint8(imdilate(tr_frame,strel('disk',3))*128));
im3(:,:,1)=imadd(im3(:,:,1), uint8(imdilate(tr_frame,strel('disk',3))*128));
im3(:,:,2)=imadd(im3(:,:,2), uint8(imdilate(tr_frame,strel('disk',3))*128));
im3(:,:,3)=imadd(im3(:,:,3), uint8(imdilate(tr_frame,strel('disk',3))*128));

% if no_tr_region==2
%     im2(start_point:end_point,:,1)=imadd(im2(start_point:end_point,:,1), uint8(immultiply(imdilate(tr_frame2,strel('disk',3)),128)));
%     im2(start_point:end_point,:,2)=imadd(im2(start_point:end_point,:,2), uint8(immultiply(imdilate(tr_frame2,strel('disk',3)),128)));
%     im2(start_point:end_point,:,3)=imadd(im2(start_point:end_point,:,3), uint8(immultiply(imdilate(tr_frame2,strel('disk',3)),128)));
%     im4(:,:,1)=imadd(im4(:,:,1), uint8(immultiply(imdilate(tr_frame2,strel('disk',3)),128)));
%     im4(:,:,2)=imadd(im4(:,:,2), uint8(immultiply(imdilate(tr_frame2,strel('disk',3)),128)));
%     im4(:,:,3)=imadd(im4(:,:,3), uint8(immultiply(imdilate(tr_frame2,strel('disk',3)),128)));
% end

if no_tr_region==1
%     [It Jt Vt]=find(trabecula1);
    [It Jt Vt]=find(trabecula_for_cortex);      % March 15, 2016
    eval(['imwrite(im1,''',direc,c1,'_roi.jpg'',''jpg'',''quality'',100);']);
    eval(['imwrite(im3(max(1,',num2str(min(It)-50),'):min(size(im3,1),',num2str(max(It)+50),'),max(1,',num2str(min(Jt)-50),'):',num2str(max(Jt)+50),',:),''',direc,c1,'_roi_inside.jpg'',''jpg'',''quality'',100);']);
%     eval(['imwrite(im1,''',archive_d,delimeter,c,nnn,'_',num2str(s),'-ro.jpg'',''jpg'',''quality'',100);']);
    eval(['imwrite(im1,''',archive_d,c1,'_ro.jpg'',''jpg'',''quality'',100);']);
    eval(['imwrite(im3(',num2str(min(It)-50),':min(size(im3,1),',num2str(max(It)+50),'),',num2str(min(Jt)-50),':',num2str(max(Jt)+50),',:),''',archive_d,c1,'_ri.jpg'',''jpg'',''quality'',100);']);
end
% if no_tr_region==2
%     [It Jt Vt]=find(trabecula1);
%     eval(['imwrite(im1,''',direc,c1,'_roi.jpg'',''jpg'',''quality'',100);']);
%     eval(['imwrite(im3(',num2str(min(It)-50),':min(size(im3,1),',num2str(max(It)+50),'),',num2str(min(Jt)-50),':',num2str(max(Jt)+50),',:),''',direc,c1,'_roi_inside.jpg'',''jpg'',''quality'',100);']);
% %     eval(['imwrite(im1,''',archive_d,delimeter,c,nnn,exp_type,'-ro_',num2str(s),'.jpg'',''jpg'',''quality'',100);']);
% %     eval(['imwrite(im3(',num2str(min(It)-50),':min(size(im3,1),',num2str(max(It)+50),'),',num2str(min(Jt)-50),':',num2str(max(Jt)+50),',:),''',archive_d,delimeter,c,nnn,exp_type,'-ri_',num2str(s),'.jpg'',''jpg'',''quality'',100);']);
%     eval(['imwrite(im1,''',archive_d,delimeter,c,nnn,'-ro_',num2str(s),'.jpg'',''jpg'',''quality'',100);']);
%     eval(['imwrite(im3(',num2str(min(It)-50),':min(size(im3,1),',num2str(max(It)+50),'),',num2str(min(Jt)-50),':',num2str(max(Jt)+50),',:),''',archive_d,delimeter,c,nnn,'_',num2str(s),'-ri.jpg'',''jpg'',''quality'',100);']);
%     [It Jt Vt]=find(trabecula2);
%     eval(['imwrite(im2,''',direc,c1,'_roi2.jpg'',''jpg'',''quality'',100);']);
%     eval(['imwrite(im4(',num2str(min(It)-50),':min(size(im4,1),',num2str(max(It)+50),'),',num2str(min(Jt)-50),':',num2str(max(Jt)+50),',:),''',direc,c1,'_roi2_inside.jpg'',''jpg'',''quality'',100);']);
%     %eval(['imwrite(im2,''',archive_d,delimeter,c1,'_roi2.jpg'',''jpg'',''quality'',100);']);
%     %eval(['imwrite(im4(',num2str(min(It)-50),':min(size(im4,1),',num2str(max(It)+50),'),',num2str(min(Jt)-50),':',num2str(max(Jt)+50),',:),''',archive_d,delimeter,c1,'_roi2_inside.jpg'',''jpg'',''quality'',100);']);
% end


%if write==1
%    if exp_type=='o'
%        eval(['imwrite(im, ''',direc,b_t,b_t,b_n,exp_t,b_t1,b_n,common,'.jpg_Files',b_t1,b_n,common,'_c_roi.jpg'',''jpg'',''quality'',100);'])
%    else
%        eval(['imwrite(im, ''',direc,b_t,b_t,b_n,b_t1,b_n,common,'.jpg_Files',b_t1,b_n,common,'_c_roi.jpg'',''jpg'',''quality'',100);'])
%    end
%    %imwrite(im,'C:\Users\shhong\Desktop\CJake\Analysis s2 3\JW3_WT\JW3_WT_5FR_v2\JW3_WT_5FR_fix_c.jpg_Files\JW3_WT_5FR_fix_c_pseudo.jpg','jpg','quality',100);
%end

im=im3;
% clear im1 im2 im3


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   12/27/07
%   test for the comparison of real analysis
%trabecula=zeros(size(armt));
%trabecula(250:1050,300:1420)=ones(1050-250+1,1420-300+1);
%tissue_volume=(round(275*dis_5x/400))^2*16;             % round(275*dis_5x/400) : side length of a block in pixel(to search trabecula manually)
                                                        % (.)^2 : area of a block in pixel
                                                        % (.)^2 *16 : 16 blocks
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 
% % For corteces
% %
% % February 28, 2013
% [x y]=find(tr);
% m_cortical=cortical(1:max(x)-dis_10x,:);
% temp=imdilate(m_cortical,strel('disk',20));
% m_cortical_frame=temp-imerode(temp,strel('disk',3));
% 
% % % for full size im
% % im(start_point+dis_10x:start_point+size(m_cortical,1)-1+dis_10x,:,1)=imadd(im(start_point+dis_10x:start_point+size(m_cortical,1)-1+dis_10x,:,1),uint8(m_cortical_frame*128));
% % im(start_point+dis_10x:start_point+size(m_cortical,1)-1+dis_10x,:,2)=imadd(im(start_point+dis_10x:start_point+size(m_cortical,1)-1+dis_10x,:,2),uint8(m_cortical_frame*128));
% % im(start_point+dis_10x:start_point+size(m_cortical,1)-1+dis_10x,:,3)=imadd(im(start_point+dis_10x:start_point+size(m_cortical,1)-1+dis_10x,:,3),uint8(m_cortical_frame*128));
% 
% % for image from the start_point
% im(dis_10x:size(m_cortical,1)-1+dis_10x,:,1)=imadd(im(dis_10x:size(m_cortical,1)-1+dis_10x,:,1),uint8(m_cortical_frame*128));
% im(dis_10x:size(m_cortical,1)-1+dis_10x,:,2)=imadd(im(dis_10x:size(m_cortical,1)-1+dis_10x,:,2),uint8(m_cortical_frame*128));
% im(dis_10x:size(m_cortical,1)-1+dis_10x,:,3)=imadd(im(dis_10x:size(m_cortical,1)-1+dis_10x,:,3),uint8(m_cortical_frame*128));
% % 
% % For corteces
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


se=strel('disk',3);

'SignalToolbox License test start'
[l m]=license('checkout','signal_Toolbox');
while (~l)
   pause(1);
   [l m]=license('checkout','signal_Toolbox');
end
'SignalToolbox License test end'

r_break=break_signal(imfill(imclose(r,se),'holes'));
g_break=break_signal(imfill(imclose(g1,se),'holes'));

for tr_loop=1%:no_tr_region
    if no_tr_region==2
        if tr_loop==1
            trabecula=trabecula1;
        elseif tr_loop==2
            trabecula=trabecula2;
        end
    elseif no_tr_region==1
        trabecula=trabecula1;
    end
    
    
[It Jt Vt]=find(trabecula);


dilate_no=5;
dilate_no1=0;
test_DIC1=trabecula & result; test_DIC1=uint8(test_DIC1)*255;
test_DIC1=imclose(test_DIC1,strel('disk',dilate_no));
[L n_test_DIC1]=bwlabel(test_DIC1);
stats_test_DIC1=regionprops(L,'Area','Image','FilledImage','BoundingBox','Orientation'); clear L
test_DIC=zeros(size(test_DIC1));
for i=1:n_test_DIC1
    start_x=round(stats_test_DIC1(i).BoundingBox(2));
    start_y=round(stats_test_DIC1(i).BoundingBox(1));
    end_x=start_x+stats_test_DIC1(i).BoundingBox(4)-1;
    end_y=start_y+stats_test_DIC1(i).BoundingBox(3)-1;
%    test_DIC(start_x:end_x,start_y:end_y)=test_DIC(start_x:end_x,start_y:end_y)+stats_test_DIC1(i).FilledImage;
    test_DIC(start_x:end_x,start_y:end_y)=test_DIC(start_x:end_x,start_y:end_y)|stats_test_DIC1(i).Image;
end

for loop=1:1
    test_DIC=medfilt2(test_DIC,[5,5]);
end
%for loop=1:5
%    test_DIC=1-medfilt2(1-test_DIC,[5,5]);
%end


[L n_DIC]=bwlabel(test_DIC);
stats_DIC=regionprops(L,'Area','Image','FilledImage','BoundingBox','Orientation'); clear L


%
%   3/13
%
%J=find([stats_DIC.Area]/sum([stats_DIC.Area])*100>0.3);
image_in=false(size(test_DIC));
% gg=(g>graythresh(g));
gg=g;
%for i=1:n_DIC
for i=1:n_DIC
    start_x=round(stats_DIC(i).BoundingBox(2));
    start_y=round(stats_DIC(i).BoundingBox(1));
    end_x=start_x+stats_DIC(i).BoundingBox(4)-1;
    end_y=start_y+stats_DIC(i).BoundingBox(3)-1;
    %if (length(find(r(start_x:end_x,start_y:end_y)&stats_DIC(i).FilledImage))>0 ...
    %        | length(find(gg(start_x:end_x,start_y:end_y)&stats_DIC(i).FilledImage))>0) ...
    %        | stats_DIC(i).Area/sum([stats_DIC.Area])*100>0.1
    %    image_in(start_x:end_x,start_y:end_y)=image_in(start_x:end_x,start_y:end_y)+stats_DIC(i).FilledImage;
    if (length(find(r(start_x:end_x,start_y:end_y)&stats_DIC(i).Image))>0 ...
            | length(find(gg(start_x:end_x,start_y:end_y)&stats_DIC(i).Image))>0) ...
            | stats_DIC(i).Area/sum([stats_DIC.Area])*100>0.1
        image_in(start_x:end_x,start_y:end_y)=image_in(start_x:end_x,start_y:end_y)|stats_DIC(i).Image;
    end
end
%J=find([stats_DIC.Area]/sum([stats_DIC.Area])*100>0.3);
%image_in=zeros(size(test_DIC));
%%for i=1:n_DIC
%for i=1:length(J)
%    start_x=round(stats_DIC(J(i)).BoundingBox(2));
%    start_y=round(stats_DIC(J(i)).BoundingBox(1));
%    end_x=start_x+stats_DIC(J(i)).BoundingBox(4)-1;
%    end_y=start_y+stats_DIC(J(i)).BoundingBox(3)-1;
%    image_in(start_x:end_x,start_y:end_y)=image_in(start_x:end_x,start_y:end_y)+stats_DIC(J(i)).FilledImage;
%end
%
%   3/13
%

image_in(end+down_margin,end)=0;

test_DIC=image_in;

% %%%%%%%%%%%%%%%%%%%%%
% %
% %   For perimeter calculation
% %   August 30, 2013
% %
% se=[0 1 0;1 0 1;0 1 0];           % structure : diamond with empty center
% test_DIC_temp=double(test_DIC>0);
% test_DIC_temp=xcorr2(test_DIC_temp,se);
% test_DIC_temp=test_DIC_temp(2:end-1,2:end-1);
% %
% %%%%%%%%%%%%%%%%%%%%%
test_DIC_temp=test_DIC;

trabecula3=trabecula;
trabecula3(end+down_margin,end)=0;
if tr_loop==1
    r(end+down_margin,end)=0;
    g(end+down_margin,end)=0;
    r_break(end+down_margin,end)=0;
    g_break(end+down_margin,end)=0;
    G(end+down_margin,end)=0;
    if trap_flag==0
        TRAP=false(size(r));
    else
        TRAP(end+down_margin,end)=0;
    end
    if ap_flag==0
        AP=false(size(r));
    else
        AP(end+down_margin,end)=0;
    end    
end

%test_r=test_DIC & (r>graythresh(r)); test_r=immultiply(uint8(test_r),255);
% test_r=(test_DIC>0) & (r_break>0); test_r=immultiply(uint8(test_r),255);
% test_r=imclose(test_r,strel('disk',dilate_no1));
% test_g=(test_DIC>0) & (g_break>graythresh(g_break)); test_g=immultiply(uint8(test_g),255);
% test_g=imclose(test_g,strel('disk',dilate_no1));
% test_GFP=(test_DIC>0) & (G>0); test_GFP=immultiply(uint8(test_GFP),255);
% test_GFP=imclose(test_GFP,strel('disk',dilate_no1));
% test_TRAP=trabecula3 & (TRAP>0); test_TRAP=immultiply(uint8(test_TRAP),255);
% test_TRAP=imclose(test_TRAP,strel('disk',dilate_no1));
% test_AP=trabecula3 & (AP>0); test_AP=immultiply(uint8(test_AP),255);
% test_AP=imclose(test_AP,strel('disk',dilate_no1));
test_r=test_DIC & r_break; test_r=uint8(test_r)*255;
test_r=imclose(test_r,strel('disk',dilate_no1));
test_g=test_DIC & g_break; test_g=uint8(test_g)*255;
test_g=imclose(test_g,strel('disk',dilate_no1));
test_GFP=test_DIC & G; test_GFP=uint8(test_GFP)*255;
test_GFP=imclose(test_GFP,strel('disk',dilate_no1));
test_TRAP=trabecula3 & TRAP; test_TRAP=uint8(test_TRAP)*255;
test_TRAP=imclose(test_TRAP,strel('disk',dilate_no1));
test_AP=trabecula3 & AP; test_AP=uint8(test_AP)*255;
test_AP=imclose(test_AP,strel('disk',dilate_no1));

%im(:,:,3)=uint8(test_DIC*255);
%im(:,:,1)=test_r;
%im(:,:,2)=imadd(test_g,test_GFP);
%figure;imshow(uint8(im),[])
clear im
% im(:,:,3)=imadd(test_AP,uint8(test_DIC*128));
% 
% 
% %if b_t==1
%     im(:,:,1)=imadd(test_r,imadd(test_TRAP,uint8(test_DIC*128)));
%     im(:,:,2)=imadd(imadd(test_g,imadd(test_TRAP,immultiply(test_GFP,0.5))),uint8(test_DIC*128));
    
im(:,:,3)=imadd(test_GFP,uint8(test_DIC*64));


%if b_t==1
%     im(:,:,1)=imadd(test_r,imadd(test_TRAP,imadd(uint8(test_DIC*64),immultiply(test_AP,0.6))));
%     im(:,:,2)=imadd(imadd(test_g,imadd(test_TRAP,test_GFP)),imadd(uint8(test_DIC*64),immultiply(test_AP,0.3)));
%%%%%%%%%%%%%%%%%%%%%
%
%	June, 7th 2018
%
	if ~RG_label	% 1st label : red,  2nd label : green --> Reversed labels
		im(:,:,1)=imadd(uint8(g&test_DIC)*255,imadd(test_TRAP,imadd(uint8(test_DIC*64),test_AP)));
		im(:,:,2)=imadd(imadd(uint8(r&test_DIC)*255,imadd(test_TRAP,test_GFP)),imadd(uint8(test_DIC)*64,immultiply(test_AP,0.5)));
	else			% 1st label : green,  2nd label : red --> Normal labels
		im(:,:,1)=imadd(uint8(r&test_DIC)*255,imadd(test_TRAP,imadd(uint8(test_DIC*64),test_AP)));
		im(:,:,2)=imadd(imadd(uint8(g&test_DIC)*255,imadd(test_TRAP,test_GFP)),imadd(uint8(test_DIC)*64,immultiply(test_AP,0.5)));
	end
%
%	June, 7th 2018
%
%%%%%%%%%%%%%%%%%%%%%

%    im(:,:,1)=imadd(uint8(r&test_DIC)*255,imadd(test_TRAP,imadd(uint8(test_DIC*64),test_AP)));
%    im(:,:,2)=imadd(imadd(uint8(g&test_DIC)*255,imadd(test_TRAP,test_GFP)),imadd(uint8(test_DIC)*64,immultiply(test_AP,0.5)));
%elseif b_t==2
%    im(:,:,2)=imadd(imadd(test_r,test_g),test_GFP);
%    im(:,:,1)=test_r;
%elseif b_t==3
%    im(:,:,2)=imadd(test_g,test_GFP);
%    im(:,:,1)=test_r;
%elseif b_t==4
%    im(:,:,2)=imadd(test_r,test_GFP);
%    im(:,:,1)=test_g;
%end

%figure;imshow(im(1:end-down_margin,:,:),[])


if write==1
%    temp_image=imrotate(im(1:end-down_margin,:,:),-rot_angle);
    temp_image=im(1:end-down_margin,:,:);
    if tr_loop==1
        eval(['imwrite(temp_image, ''',direc,c1,'_pseudo1_tr1.jpg'',''jpg'',''quality'',100);']);
        eval(['imwrite(temp_image(',num2str(min(It)-50),':min(size(temp_image,1),',num2str(max(It)+50),'),',num2str(min(Jt)-50),':',num2str(max(Jt)+50),',:),''',direc,c1,'_pseudo1_tr1_inside.jpg'',''jpg'',''quality'',100);']);
        %eval(['imwrite(temp_image, ''',archive_d,delimeter,c1,'_pseudo1_tr1.jpg'',''jpg'',''quality'',100);']);
%         eval(['imwrite(temp_image(',num2str(min(It)-50),':min(size(temp_image,1),',num2str(max(It)+50),'),',num2str(min(Jt)-50),':',num2str(max(Jt)+50),',:),''',archive_d,delimeter,c,nnn,exp_type,'-ps_',num2str(s),'.jpg'',''jpg'',''quality'',100);']);
%         eval(['imwrite(temp_image(',num2str(min(It)-50),':min(size(temp_image,1),',num2str(max(It)+50),'),',num2str(min(Jt)-50),':',num2str(max(Jt)+50),',:),''',archive_d,delimeter,c,nnn,'_',num2str(s),'-ps.jpg'',''jpg'',''quality'',100);']);
        eval(['imwrite(temp_image(',num2str(min(It)-50),':min(size(temp_image,1),',num2str(max(It)+50),'),',num2str(min(Jt)-50),':',num2str(max(Jt)+50),',:),''',archive_d,c1,'_ps.jpg'',''jpg'',''quality'',100);']);
    elseif tr_loop==2
        eval(['imwrite(temp_image, ''',direc,c1,'_pseudo1_tr2.jpg'',''jpg'',''quality'',100);']);
        eval(['imwrite(temp_image(',num2str(min(It)-50),':min(size(temp_image,1),',num2str(max(It)+50),'),',num2str(min(Jt)-50),':',num2str(max(Jt)+50),',:),''',direc,c1,'_pseudo1_tr2_inside.jpg'',''jpg'',''quality'',100);']);
        %eval(['imwrite(temp_image, ''',archive_d,delimeter,c1,'_pseudo1_tr2.jpg'',''jpg'',''quality'',100);']);
        %eval(['imwrite(temp_image(',num2str(min(It)-50),':min(size(temp_image,1),',num2str(max(It)+50),'),',num2str(min(Jt)-50),':',num2str(max(Jt)+50),',:),''',archive_d,delimeter,c1,'_pseudo1_tr2_inside.jpg'',''jpg'',''quality'',100);']);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   imwrite middle point signals
%   March 02, 2013
im(:,:,3)=imadd(test_GFP,uint8(test_DIC)*64);
% im(:,:,1)=imadd(uint8(test_r&test_DIC)*255,imadd(test_TRAP,imadd(uint8(test_DIC)*64,test_AP)));
% im(:,:,2)=imadd(uint8(test_g&test_DIC)*255,imadd(test_TRAP,imadd(test_GFP,imadd(uint8(test_DIC)*64,immultiply(test_AP,0.3)))));
%%%%%%%%%%%%%%%%%%%%%
%
%	June, 7th 2018
%
if ~RG_label	% 1st label : red,  2nd label : green --> Reversed labels
	im(:,:,1)=imadd(uint8(test_g&test_DIC)*255,imadd(test_TRAP,imadd(uint8(test_DIC)*64,test_AP)));
	im(:,:,2)=imadd(uint8(test_r&test_DIC)*255,imadd(test_TRAP,imadd(test_GFP,imadd(uint8(test_DIC)*64,immultiply(test_AP,0.3)))));
else			% 1st label : green,  2nd label : red --> Normal labels
	im(:,:,1)=imadd(uint8(test_r&test_DIC)*255,imadd(test_TRAP,imadd(uint8(test_DIC)*64,test_AP)));
	im(:,:,2)=imadd(uint8(test_g&test_DIC)*255,imadd(test_TRAP,imadd(test_GFP,imadd(uint8(test_DIC)*64,immultiply(test_AP,0.3)))));
end
%
%	June, 7th 2018
%
%%%%%%%%%%%%%%%%%%%%%

if write==1
    temp_image=im(1:end-down_margin,:,:);
    if tr_loop==1
        eval(['imwrite(temp_image, ''',direc,c1,'_pseudo1_tr1_mid.jpg'',''jpg'',''quality'',100);']);
        eval(['imwrite(temp_image(',num2str(min(It)-50),':min(size(temp_image,1),',num2str(max(It)+50),'),',num2str(min(Jt)-50),':',num2str(max(Jt)+50),',:),''',direc,c1,'_pseudo1_tr1_mid_inside.jpg'',''jpg'',''quality'',100);']);
%         eval(['imwrite(temp_image(',num2str(min(It)-50),':min(size(temp_image,1),',num2str(max(It)+50),'),',num2str(min(Jt)-50),':',num2str(max(Jt)+50),',:),''',archive_d,delimeter,c,nnn,'_',num2str(s),'-ps.jpg'',''jpg'',''quality'',100);']);
    elseif tr_loop==2
        eval(['imwrite(temp_image, ''',direc,c1,'_pseudo1_tr2_mid.jpg'',''jpg'',''quality'',100);']);
        eval(['imwrite(temp_image(',num2str(min(It)-50),':min(size(temp_image,1),',num2str(max(It)+50),'),',num2str(min(Jt)-50),':',num2str(max(Jt)+50),',:),''',direc,c1,'_pseudo1_tr2_mid_inside.jpg'',''jpg'',''quality'',100);']);
    end
end
%
%   March 02, 2013
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

toc


%eval(['save ''C:\Users\shhong\Desktop\Jake New\11473-1\IP-ConvertImage-06\trabecula_home'' trabecula cortical im test_DIC']);
if write==1
    %if b_t==1
    %if exp_type=='o'
        if with_rg_label==1
            if tr_loop==1
                eval(['save ''',direc,c1,'_trabecula_bone1_tr1.mat'' trabecula cortical im test_DIC']);
            elseif tr_loop==2
                eval(['save ''',direc,c1,'_trabecula_bone1_tr2.mat'' trabecula cortical im test_DIC']);
            end
        elseif with_rg_label==0
            if tr_loop==1
                eval(['save ''',direc,c1,'_trabecula_bone1_tr1_f.mat'' trabecula cortical im test_DIC']);
            elseif tr_loop==2
                eval(['save ''',direc,c1,'_trabecula_bone1_tr2_f.mat'' trabecula cortical im test_DIC']);
            end
        end
    %else
    %    if with_rg_label==1
    %        eval(['save ''',direc,b_t,b_t,b_n,b_t1,b_n,common,'.jpg_Files\trabecula_bone1'' trabecula cortical im test_DIC -V6'])
    %    elseif with_rg_label==0
    %        eval(['save ''',direc,b_t,b_t,b_n,b_t1,b_n,common,'.jpg_Files\trabecula_bone1_f'' trabecula cortical im test_DIC -V6'])
    %    end
    %end
end

%eval(['save ''C:\Users\shhong\Desktop\CJake\Analysis s2 3\JW3_WT\JW3_WT_5FR_v2\JW3_WT_5FR_fix_c.jpg_Files\trabecula_bone'' trabecula cortical im test_DIC -V6']);

%clear im trabecula cortical a* endosteum periosteum osteum left* right*

%image_in=double(test_DIC);
%image_in=image_in/max(image_in(:));

% [average_height, average_width, m_cortical, osteocytes_number, osteocytes_area]=analyze_cortices(tr, cortical, DAPI, dis_10x);
[average_height, average_width, m_cortical, osteocytes_number, osteocytes_area]=analyze_cortices(trabecula_for_cortex, cortical, DAPI, dis_10x);    % March 15, 2016
analysis.roi.height=average_height*dis_10x_for_400micron_per_312pixels;
analysis.roi.width=average_width*dis_10x_for_400micron_per_312pixels;
analysis.roi.ratio=analysis.roi.height/analysis.roi.width;
analysis.cortex.area=length(find(m_cortical))*dis_10x_for_400micron_per_312pixels^2;
analysis.cortex.osteocytes.number=osteocytes_number;
analysis.cortex.osteocytes.number_per_cortex_area=osteocytes_number/analysis.cortex.area;
analysis.cortex.osteocytes.number_per_100micron_height=osteocytes_number/analysis.cortex.area*(100/dis_10x_for_400micron_per_312pixels)/analysis.roi.height;
analysis.cortex.osteocytes.area=length(find(osteocytes_area))*dis_10x_for_400micron_per_312pixels^2;
analysis.cortex.osteocytes.area_per_cortex_area=length(find(osteocytes_area))/length(find(m_cortical));
analysis.cortex.area=length(find(m_cortical))*dis_10x_for_400micron_per_312pixels^2;
analysis.cortex.thickness=analysis.cortex.area/analysis.roi.height/2; % average cortex width

temp=imdilate(m_cortical,strel('disk',20));
m_cortical_frame=temp-imerode(temp,strel('disk',3)); clear temp

% for full size im
im1(start_point+dis_10x:start_point+size(m_cortical,1)-1+dis_10x,:,1)=imadd(im1(start_point+dis_10x:start_point+size(m_cortical,1)-1+dis_10x,:,1),uint8(m_cortical_frame)*128);
im1(start_point+dis_10x:start_point+size(m_cortical,1)-1+dis_10x,:,2)=imadd(im1(start_point+dis_10x:start_point+size(m_cortical,1)-1+dis_10x,:,2),uint8(m_cortical_frame)*128);
im1(start_point+dis_10x:start_point+size(m_cortical,1)-1+dis_10x,:,3)=imadd(im1(start_point+dis_10x:start_point+size(m_cortical,1)-1+dis_10x,:,3),uint8(m_cortical_frame)*128);

% for image from the start_point
% im(dis_10x:size(m_cortical,1)-1+dis_10x,:,1)=imadd(im(dis_10x:size(m_cortical,1)-1+dis_10x,:,1),uint8(m_cortical_frame*128));
% im(dis_10x:size(m_cortical,1)-1+dis_10x,:,2)=imadd(im(dis_10x:size(m_cortical,1)-1+dis_10x,:,2),uint8(m_cortical_frame*128));
% im(dis_10x:size(m_cortical,1)-1+dis_10x,:,3)=imadd(im(dis_10x:size(m_cortical,1)-1+dis_10x,:,3),uint8(m_cortical_frame*128));
% im3(:,:,1)=imadd(im3(:,:,1), uint8(immultiply(imdilate(m_cortical_frame,strel('disk',3)),128)));
% im3(:,:,2)=imadd(im3(:,:,2), uint8(immultiply(imdilate(m_cortical_frame,strel('disk',3)),128)));
% im3(:,:,3)=imadd(im3(:,:,3), uint8(immultiply(imdilate(m_cortical_frame,strel('disk',3)),128)));

% [It Jt Vt]=find(trabecula1);
eval(['imwrite(im1,''',direc,c1,'_roi1.jpg'',''jpg'',''quality'',100);']);
% eval(['imwrite(im3(',num2str(min(It)-50),':min(size(im3,1),',num2str(max(It)+50),'),',num2str(min(Jt)-50),':',num2str(max(Jt)+50),',:),''',direc,c1,'_roi_inside1.jpg'',''jpg'',''quality'',100);']);
% eval(['imwrite(im1,''',archive_d,delimeter,c,nnn,'_',num2str(s),'-ro1.jpg'',''jpg'',''quality'',100);']);
% eval(['imwrite(im3(',num2str(min(It)-50),':min(size(im3,1),',num2str(max(It)+50),'),',num2str(min(Jt)-50),':',num2str(max(Jt)+50),',:),''',archive_d,delimeter,c,nnn,'_',num2str(s),'-ri1.jpg'',''jpg'',''quality'',100);']);




%
%
%   Jan 28. 2010
%
[analysis.perimeter.DIC, analysis.area.DIC] = perimeter_area1(test_DIC_temp, image_in-imerode(image_in,strel('disk',1)),dis_10x_for_400micron_per_312pixels, exact_boundary);

analysis.area.DIC=length(find(image_in))*dis_10x_for_400micron_per_312pixels^2;


%if red_flag==1
%    [out_label]=outside_trim_label(image_in, double(test_r)/255, 20, exact_boundary);
%    [in_label]=inside_trim_label(image_in, double(test_r)/255, loop_count, exact_boundary);
%    red_label=out_label | in_label; clear in_label out_label
%    
%    temp(:,:,1)=immultiply(uint8(red_label),127);
%%    temp(:,:,3)=;
%    % Mar/22 eval(['imwrite(temp, ''',bone_name,'red_label_trimmed_threshold_',num2str(loop_count),'.jpg'', ''jpg'', ''quality'', 100);']); clear temp
%    %eval(['imwrite(temp, ''',c_path,'\',bone_name,'\',bone_name,'_red_label_trimmed_threshold.jpg'', ''jpg'', ''quality'', 100);']); clear temp
%else
%    red_label=uint8(zeros(size(ref_bone)));
%end
%[analysis.perimeter.red, analysis.area.red] = perimeter_area1(test_DIC_temp, red_label);


roi=trabecula3-imerode(trabecula3,strel('disk',1));             % July 18, 2012
roi=(roi>0);                                                    % July 18, 2012

if GFP_flag==1
    [out_label]=outside_trim_label(image_in, double(test_GFP)/255, 20, exact_boundary);
    [in_label]=inside_trim_label(image_in, double(test_GFP)/255, loop_count, exact_boundary);
    GFP_label=out_label | in_label; clear in_label out_label
    
%     temp(:,:,2)=immultiply(uint8(GFP_label),127);
%     temp(:,:,3)=immultiply(uint8(GFP_label),127);
%    temp(:,:,3)=;
    % Mar/22 eval(['imwrite(temp, ''',bone_name,'red_label_trimmed_threshold_',num2str(loop_count),'.jpg'', ''jpg'', ''quality'', 100);']); clear temp
    %eval(['imwrite(temp, ''',c_path,'\',bone_name,'\',bone_name,'_red_label_trimmed_threshold.jpg'', ''jpg'', ''quality'', 100);']); clear temp
else
    GFP_label=(zeros(size(test_DIC)));
end
GFP_label=(GFP_label>0)&(1-roi);
GFP_label1=GFP_label;

if trap_flag==1
    [out_label]=outside_trim_label(image_in, double(test_TRAP)/255, 20, exact_boundary);
    [in_label]=inside_trim_label(image_in, double(test_TRAP)/255, loop_count, exact_boundary);
    TRAP_label=out_label | in_label; clear in_label out_label
    
%     temp(:,:,2)=immultiply(uint8(TRAP_label),127);
%     temp(:,:,3)=immultiply(uint8(TRAP_label),127);
%    temp(:,:,3)=;
    % Mar/22 eval(['imwrite(temp, ''',bone_name,'red_label_trimmed_threshold_',num2str(loop_count),'.jpg'', ''jpg'', ''quality'', 100);']); clear temp
    %eval(['imwrite(temp, ''',c_path,'\',bone_name,'\',bone_name,'_red_label_trimmed_threshold.jpg'', ''jpg'', ''quality'', 100);']); clear temp
else
    TRAP_label=(zeros(size(test_DIC)));
end
TRAP_label=(TRAP_label>0)&(1-roi);
%TRAP_label1=TRAP_label;
[analysis.perimeter.T, analysis.area.T] = perimeter_area1(test_DIC_temp, TRAP_label, dis_10x_for_400micron_per_312pixels, exact_boundary);

if ap_flag==1
    [out_label]=outside_trim_label(image_in, double(test_AP)/255, 20, exact_boundary);
    [in_label]=inside_trim_label(image_in, double(test_AP)/255, loop_count, exact_boundary);
    AP_label=out_label | in_label; clear in_label out_label
    
%     temp(:,:,2)=immultiply(uint8(AP_label),127);
%     temp(:,:,3)=immultiply(uint8(AP_label),127);
%    temp(:,:,3)=;
    % Mar/22 eval(['imwrite(temp, ''',bone_name,'red_label_trimmed_threshold_',num2str(loop_count),'.jpg'', ''jpg'', ''quality'', 100);']); clear temp
    %eval(['imwrite(temp, ''',c_path,'\',bone_name,'\',bone_name,'_red_label_trimmed_threshold.jpg'', ''jpg'', ''quality'', 100);']); clear temp
else
    AP_label=(zeros(size(test_DIC)));
end
AP_label=(AP_label>0)&(1-roi);
%AP_label1=AP_label;
[analysis.perimeter.A, analysis.area.A] = perimeter_area1(test_DIC_temp, AP_label, dis_10x_for_400micron_per_312pixels, exact_boundary);
AT=AP_label & TRAP_label;
[analysis.perimeter.AT, analysis.area.AT]=perimeter_area1(test_DIC_temp, AT, dis_10x_for_400micron_per_312pixels, exact_boundary);
AGFP=AP_label & GFP_label;
[analysis.perimeter.AGFP, analysis.area.AGFP]=perimeter_area1(test_DIC_temp, AGFP, dis_10x_for_400micron_per_312pixels, exact_boundary);
TGFP=TRAP_label & GFP_label;
[analysis.perimeter.TGFP, analysis.area.TGFP]=perimeter_area1(test_DIC_temp, TGFP, dis_10x_for_400micron_per_312pixels, exact_boundary);
toc

%[analysis.perimeter.GFP, analysis.area.GFP] = perimeter_area1(test_DIC_temp, GFP_label);

g_trim=zeros(size(test_DIC));
g_trim1=g_trim;
r_trim=zeros(size(test_DIC));
r_trim1=g_trim;
g_trim_wo_red=g_trim;
test1=test_DIC & test_g;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   label separation in the end of the bone and multiple labels in one area
%   December 15, 2009
%
%test1=label_separation(test1);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[L_t1 n_t1]=bwlabel(test1);
stats_t1=regionprops(L_t1,'Image','BoundingBox','Orientation','Area','Centroid'); clear L_t1
test1=zeros(size(test1));
area=[stats_t1.Area];
I=find(area>5);

for i=1:length(I)
    start_x=round(stats_t1(I(i)).BoundingBox(2));
    start_y=round(stats_t1(I(i)).BoundingBox(1));
    end_x=start_x+stats_t1(I(i)).BoundingBox(4)-1;
    end_y=start_y+stats_t1(I(i)).BoundingBox(3)-1;
    test1(start_x:end_x,start_y:end_y)=test1(start_x:end_x,start_y:end_y)+stats_t1(I(i)).Image;
end
clear stats_t1

test2=test_DIC & test_r;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   label separation in the end of the bone and multiple labels in one area
%   December 15, 2009
%
%test2=label_separation(test2);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[L_t2 n_t2]=bwlabel(test2);
stats_t2=regionprops(L_t2,'Image','BoundingBox','Orientation','Area','Centroid'); clear L_t2
test2=zeros(size(test2));
area=[stats_t2.Area];
I=find(area>5);

for i=1:length(I)
    start_x=round(stats_t2(I(i)).BoundingBox(2));
    start_y=round(stats_t2(I(i)).BoundingBox(1));
    end_x=start_x+stats_t2(I(i)).BoundingBox(4)-1;
    end_y=start_y+stats_t2(I(i)).BoundingBox(3)-1;
    test2(start_x:end_x,start_y:end_y)=test2(start_x:end_x,start_y:end_y)+stats_t2(I(i)).Image;
end
clear stats_t2

%if b_t==1 & red_first_green_last==1          % red first, green last
%    temp_test1=test1;
%    test1=test2;
%    test2=temp_test1;
%    clear temp_test1;
%end

[L n_DIC]=bwlabel(test_DIC);
stats_DIC=regionprops(L,'Area','Image','FilledImage','BoundingBox','Orientation'); clear L
J=find([stats_DIC.Area]);
thickness_g=0;
thickness_r=0;
TB=struct;
for i=1:length(J)
    %i
    if i==12
        11;
    end

    sdx=round(stats_DIC(J(i)).BoundingBox(4)/2);
    sdy=round(stats_DIC(J(i)).BoundingBox(3)/2);
    edx=round(stats_DIC(J(i)).BoundingBox(4)/2);
    edy=round(stats_DIC(J(i)).BoundingBox(3)/2);
    sd=[sdx, edx, sdy, edy];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   December, 22 2009
%   
%
%     start_x=round(stats_DIC(J(i)).BoundingBox(2))-sdx;
%     if start_x<1
%         sdx=sdx-start_x;
%         start_x=1;
%         %sxf=1;
%     end
%     start_y=round(stats_DIC(J(i)).BoundingBox(1))-sdy;
%     if start_y<1
%         sdy=sdy-start_y;
%         start_y=1;
%         %syf=1;
%     end
%     end_x=start_x+stats_DIC(J(i)).BoundingBox(4)-1+sdx+edx;
%     if end_x>size(test1,1)
%         edx=edx-(end_x-size(test1,1));
%         end_x=size(test1,1);
%         %exf=1;
%     end
%     end_y=start_y+stats_DIC(J(i)).BoundingBox(3)-1+sdy+edy;
%     if end_y>size(test1,2)
%         edy=edy-(end_y-size(test1,2));
%         end_y=size(test1,2);
%         %efy=1;
%     end
%     im=zeros(end_x-start_x+1,end_y-start_y+1);
%     im(sdx+1:end-edx,sdy+1:end-edy)=stats_DIC(J(i)).FilledImage;
%     g_test=test1(start_x:end_x,start_y:end_y)&im;
%     r_test=test2(start_x:end_x,start_y:end_y)&im;
%
    sx=round(stats_DIC(J(i)).BoundingBox(2));
    ex=sx+stats_DIC(J(i)).BoundingBox(4)-1;
    sy=round(stats_DIC(J(i)).BoundingBox(1));
    ey=sy+stats_DIC(J(i)).BoundingBox(3)-1;
    ss=[sx, ex, sy, ey];
    
    im=zeros(stats_DIC(J(i)).BoundingBox(4)+2*sdx,stats_DIC(J(i)).BoundingBox(3)+2*sdy);
    g_test=im;
    r_test=im;
%%%%%%%%%%%%%%%%%%%%%
%
%   June 3, 2010
%     im(sdx+1:sdx+stats_DIC(J(i)).BoundingBox(4),sdy+1:sdy+stats_DIC(J(i)).BoundingBox(3))=stats_DIC(J(i)).FilledImage;
%     g_test(sdx+1:sdx+stats_DIC(J(i)).BoundingBox(4),sdy+1:sdy+stats_DIC(J(i)).BoundingBox(3))=test1(sx:ex,sy:ey) & stats_DIC(J(i)).FilledImage;
%     r_test(sdx+1:sdx+stats_DIC(J(i)).BoundingBox(4),sdy+1:sdy+stats_DIC(J(i)).BoundingBox(3))=test2(sx:ex,sy:ey) & stats_DIC(J(i)).FilledImage;
%
%%%%%%%%%%%%%%%%%%%%%
    im(sdx+1:sdx+stats_DIC(J(i)).BoundingBox(4),sdy+1:sdy+stats_DIC(J(i)).BoundingBox(3))=stats_DIC(J(i)).Image;
    g_test(sdx+1:sdx+stats_DIC(J(i)).BoundingBox(4),sdy+1:sdy+stats_DIC(J(i)).BoundingBox(3))=test1(sx:ex,sy:ey) & stats_DIC(J(i)).Image;
    r_test(sdx+1:sdx+stats_DIC(J(i)).BoundingBox(4),sdy+1:sdy+stats_DIC(J(i)).BoundingBox(3))=test2(sx:ex,sy:ey) & stats_DIC(J(i)).Image;
    
%   December, 22 2009
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    trim_template=zeros(size(test_DIC));
    if red_first_green_last==1 & b_t==1         % red_first, green last
        r_test_temp=r_test;
        r_test=g_test;
        g_test=r_test_temp;
    end
    [r_trim_t, r_trim1_t, thick_r, g_trim_t, g_trim1_t, thick_g, TB]=move_control(i, ss, sd, trim_template, im, r_test, g_test, GFP_label, GFP_label1, stats_DIC, J, TB);
%     if isstruct(T)==1
%         TB(i)=T;
%     else
%         i;
%     end
    %if i==9
    %    i
    %end

    r_trim=r_trim + r_trim_t;
    r_trim1=r_trim1 + r_trim1_t;
    g_trim=g_trim + g_trim1_t;
    g_trim1=g_trim1 + g_trim1_t;
    thickness_r(i,1:size(thick_r,2))=thick_r;
    thickness_g(i,1:size(thick_g,2))=thick_g;       

end
if red_first_green_last==1 & b_t==1         % red_first, green last
    r_t=r_trim; clear t_trim
    r_trim=g_trim; clear g_trim
    g_trim=r_t;
    r_t=r_trim1; clear t_trim1
    r_trim1=g_trim1; clear g_trim1
    g_trim1=r_t;
    
    thickness_r_temp=thickness_r; clear thickness_r 
    thickness_r=thickness_g; clear thickness_g
    thickness_g=thickness_r_temp;
end

toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   March 17, 2016
%
analysis.thickness.Tb_Th=0;
if isempty(find(thickness_g))==0
    thickness_g=thickness_g(find(1-isnan(thickness_g)));
%     [Id Jd Vd]=find(thickness_g);[x y z]=find(Vd>0);
    [Idg Jdg Vdg]=find(thickness_g);[xg yg zg]=find(Vdg>0 & Vdg<120);       % March 17, 2016
%     analysis.thickness.Tb_Th=mean(Vd(x))*dis_10x_for_400micron_per_312pixels;
%     clear Id Jd Vd
else
    Vdg=[];
end
if isempty(find(thickness_r))==0
    thickness_r=thickness_r(find(1-isnan(thickness_r)));
%     [Id Jd Vd]=find(thickness_r);[x y z]=find(Vd>0);
    [Idr Jdr Vdr]=find(thickness_r);[xr yr zr]=find(Vdr>0 & Vdr<120);       % March 17, 2016
%     analysis.thickness.Tb_Th=analysis.thickness.Tb_Th+mean(Vd(x))*dis_10x_for_400micron_per_312pixels;
%     analysis.thickness.Tb_Th=analysis.thickness.Tb_Th/2;                            % January 5, 2009
%     clear thickness_g thickness_r Id Jd Vd
else
    Vdr=[];
end
if isempty(Vdr) & isempty(Vdg)
    analysis.thickness.Tb_Th=NaN;
elseif isempty(Vdr) & ~isempty(Vdg)
    analysis.thickness.Tb_Th=mean([Vdg(xg)])*dis_10x_for_400micron_per_312pixels;
elseif ~isempty(Vdr) & isempty(Vdg)
    analysis.thickness.Tb_Th=mean([Vdr(xr)])*dis_10x_for_400micron_per_312pixels;
elseif ~isempty(Vdr) & ~isempty(Vdg)
    tempg=Vdg(xg); tempr=Vdr(xr);
    analysis.thickness.Tb_Th=mean([tempg(:);tempr(:)])*dis_10x_for_400micron_per_312pixels;
end
clear thickness_g thickness_r Idr Jdr Vdr Idg Jdg Vdg
%
%   March 17, 2016
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[analysis.perimeter.GFP, analysis.area.GFP] = perimeter_area1(test_DIC_temp, GFP_label, dis_10x_for_400micron_per_312pixels, exact_boundary);

%clear g g2 g2t r  result

%[out_label]=outside_trim_label(image_in, g_trim_wo_red, 20, exact_boundary);
%[in_label]=inside_trim_label(image_in, g_trim_wo_red, loop_count, exact_boundary);
%g_trim_wo_red_label=out_label;% | in_label; clear in_label out_label
%g_label=double(in_label)+g_trim;

% roi=trabecula3-imerode(trabecula3,strel('disk',1));             % July 18, 2012
r_trim=(r_trim>0)&(1-roi);                                          % July 18, 2012
g_trim=(g_trim>0)&(1-roi);                                          % July 18, 2012
% r_trim1=(r_trim1>0)-roi;                                         % July 18, 2012
% g_trim1=(g_trim1>0)-roi;                                         % July 18, 2012

r_label=r_trim;
r_label=imdilate(r_label,strel('disk',1));
test_DIC=logical(double(test_DIC)/255);
r_label=r_label & (test_DIC-imerode(test_DIC,strel('disk',1)));
[analysis.perimeter.r, analysis.area.r] = perimeter_area1(test_DIC_temp, r_label, dis_10x_for_400micron_per_312pixels, exact_boundary);

g_label=g_trim;
g_label=uint8(logical(g_label));
g_label=imdilate(g_label,strel('disk',1));
g_label=g_label & (test_DIC-imerode(test_DIC,strel('disk',1)));
[analysis.perimeter.g, analysis.area.g] = perimeter_area1(test_DIC_temp, logical(g_label), dis_10x_for_400micron_per_312pixels, exact_boundary);

rg=r_label & logical(g_label);
[analysis.perimeter.rg, analysis.area.rg] = perimeter_area1(test_DIC_temp, rg, dis_10x_for_400micron_per_312pixels, exact_boundary);

%r_only=r_label & (1-rg);
%[analysis.perimeter.r_only, analysis.area.r_only] = perimeter_area1(test_DIC_temp, r_only);
analysis.perimeter.r_only=analysis.perimeter.r-analysis.perimeter.rg;
analysis.area.r_only=analysis.area.r-analysis.area.rg;
%g_only=g_label & (1-rg);
%[analysis.perimeter.g_only, analysis.area.g_only] = perimeter_area1(test_DIC_temp, g_only);
analysis.perimeter.g_only=analysis.perimeter.g-analysis.perimeter.rg;
analysis.area.g_only=analysis.area.g-analysis.area.rg;
%sLS=r_only | g_only ;
%[analysis.perimeter.sLS, analysis.area.sLS] = perimeter_area1(test_DIC_temp, sLS);
analysis.perimeter.sLS=analysis.perimeter.r_only+analysis.perimeter.g_only;
analysis.area.sLS=analysis.area.r_only+analysis.area.g_only;
%analysis.perimeter.LS=analysis.perimeter.sLS+analysis.perimeter.rg

rGFP=r_label & GFP_label;
[analysis.perimeter.rGFP, analysis.area.rGFP] = perimeter_area1(test_DIC_temp, rGFP, dis_10x_for_400micron_per_312pixels, exact_boundary);
gGFP=logical(g_label) & GFP_label;
[analysis.perimeter.gGFP, analysis.area.gGFP] = perimeter_area1(test_DIC_temp, gGFP, dis_10x_for_400micron_per_312pixels, exact_boundary);
rgGFP=rGFP & gGFP;
[analysis.perimeter.rgGFP, analysis.area.rgGFP] = perimeter_area1(test_DIC_temp, rgGFP, dis_10x_for_400micron_per_312pixels, exact_boundary);
analysis.surface.BS=analysis.perimeter.DIC;                                           % Bone Surface
analysis.surface.LS=analysis.perimeter.r+analysis.perimeter.g-analysis.perimeter.rg;                    % Labeled surface
analysis.surface.r_BS=analysis.perimeter.r/analysis.surface.BS*100;
analysis.surface.g_BS=analysis.perimeter.g/analysis.surface.BS*100;
analysis.surface.rg_BS=analysis.perimeter.rg/analysis.surface.BS*100;                                  % dLS/BS (double labeled surface / bone surface)
analysis.surface.rGFP_BS=analysis.perimeter.rGFP/analysis.surface.BS*100;
analysis.surface.gGFP_BS=analysis.perimeter.gGFP/analysis.surface.BS*100;
analysis.surface.rgGFP_BS=analysis.perimeter.rgGFP/analysis.surface.BS*100;
analysis.surface.GFP_BS=analysis.perimeter.GFP/analysis.surface.BS*100;                                 % GFP/BS
analysis.surface.AGFP_BS=analysis.perimeter.AGFP/analysis.surface.BS*100;                                 % GFP/BS
analysis.surface.TGFP_BS=analysis.perimeter.TGFP/analysis.surface.BS*100;                                 % GFP/BS
%analysis.surface.single_BS=(analysis.perimeter.r+analysis.perimeter.g-analysis.perimeter.rg*2)/analysis.surface.BS*100;   % sLS/BS (single labeled surface / bone surface)
%analysis.surface.single_LS=(analysis.perimeter.r+analysis.perimeter.g-analysis.perimeter.rg*2)/analysis.surface.LS*100;   % sLS/LS (single labeled surface / labeled surface)
analysis.surface.r_only_BS=analysis.perimeter.r_only/analysis.surface.BS*100;   % sLS(red only)/BS (red only labeled surface / bone surface)
analysis.surface.g_only_BS=analysis.perimeter.g_only/analysis.surface.BS*100;   % sLS(green only)/BS (green only labeled surface / bone surface)
analysis.surface.sLS_BS=analysis.perimeter.sLS/analysis.surface.BS*100;   % sLS/BS (single labeled surface / bone surface)
if analysis.surface.LS~=0
    analysis.surface.sLS_LS=analysis.perimeter.sLS/analysis.surface.LS*100;   % sLS/LS (single labeled surface / labeled surface)
    analysis.surface.rg_LS=analysis.perimeter.rg/analysis.surface.LS*100;                                   % dLS/LS (double labeled surface / labeled surface)
elseif analysis.surface.LS==0
    analysis.surface.sLS_LS=0;   % sLS/LS (single labeled surface / labeled surface)
    analysis.surface.rg_LS=0;                                   % dLS/LS (double labeled surface / labeled surface)
elseif isnan(analysis.surface.LS)
    analysis.surface.sLS_LS=NaN;   % sLS/LS (single labeled surface / labeled surface)
    analysis.surface.rg_LS=NaN;
end
analysis.surface.LS_BS=analysis.surface.LS/analysis.surface.BS*100;
analysis.surface.dLS_sLS=analysis.perimeter.rg/analysis.perimeter.sLS*100;   % sLS/dLS (single labeled surface / double labeled surface)

analysis.surface.MS_BS=(analysis.surface.r_only_BS+analysis.surface.g_only_BS)/2+analysis.surface.rg_BS;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   modify after middle point is introduced
%   March 2, 2013
%   
% t=g_trim1 & r_trim1;
% if red_first_green_last==1 & b_t==1
%     tt=t.*r_trim1-t.*g_trim1;               % Jan. 5, 2010
% else
%     tt=t.*g_trim1-t.*r_trim1;
% end
% [dI dJ dV]=find(tt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   New MAR using only double label
%   July 14th, 2015
%
% if red_first_green_last==1 & b_t==1
%     tt=r_trim1-g_trim1;               % Jan. 5, 2010
% else
%     tt=g_trim1-r_trim1;
% end
trim_template=(g_trim1>0)&(r_trim1>0);      % double label
if red_first_green_last==1
    tt=(r_trim1-g_trim1).*double(trim_template);
else
    tt=(g_trim1-r_trim1).*double(trim_template);
end
%
%   New MAR using only double_label
%   July 14th, 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[dI dJ dV]=find(tt);
sdV=sort(dV);%*400/dis_10x;
%
%   modify after middle point is introduced
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%isdV=find(sdV<0);
%test_dv=sdV(length(isdV)+1:end);
%middle=round(length(test_dv)/2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   June 3, 2010 : ignore the part after the big jump of test_dv1
%test_dv=sdV(find(sdV>0));
test_dv1=sdV(find(sdV>0));
%end_dv=find(diff(test_dv1)>mean(diff(test_dv1))*5,1);
end_dv=find(diff(test_dv1)>5,1);
if isempty(end_dv)==1
    test_dv=test_dv1*dis_10x_for_400micron_per_312pixels;
else
    test_dv=test_dv1(1:end_dv)*dis_10x_for_400micron_per_312pixels;   % added
end%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   modify after middle point is introduced
%   March 2, 2013
%   
% middle=round(length(test_dv)/2);
[x y]=find(test_dv1==mode(test_dv1));
if ~isempty(x)
    middle=min(x(round(length(x)/2)), round(length(test_dv)/2));
%
%   modify after middle point is introduced
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %analysis.thickness.Ir_L_Th=mean(sdV(length(isdV)+1:end));          % Inter Label Thickness
    %analysis.thickness.Ir_L_Th=median(sdV(length(isdV)+1:end));          % Inter Label Thickness
    if middle<4
    %     analysis.thickness.Ir_L_Th=mean(sdV(length(isdV)+1:end));
    %     Ir_L_Th_std=std(sdV(length(isdV)+1:end));
        analysis.thickness.Ir_L_Th=mean(test_dv);
        Ir_L_Th_std=std(test_dv);
    else
    %    analysis.thickness.Ir_L_Th=mean(sdV(middle-floor(middle/2):middle+ceil(middle/2)));
    %    Ir_L_Th_std=std(sdV(middle-floor(middle/2):middle+ceil(middle/2)));
    %    analysis.thickness.Ir_L_Th=mean(test_dv(middle-floor(middle/2):middle+ceil(middle/2)));
    %    Ir_L_Th_std=std(test_dv(middle-floor(middle/2):middle+ceil(middle/2)));
    %     analysis.thickness.Ir_L_Th=mean(test_dv(middle-floor(middle/2):middle+ceil(middle/2)));
        analysis.thickness.Ir_L_Th=mean(test_dv(middle-floor(middle/2):middle+ceil(middle/4)));       %     2/3/2016
        Ir_L_Th_std=std(test_dv(middle-floor(middle/2):middle+ceil(middle/4)));                       %     2/3/2016
        analysis.thickness.Ir_L_Th=mean(test_dv(1:middle+ceil(middle/2)));                              %     2/3/2016
        Ir_L_Th_std=std(test_dv(1:middle+ceil(middle/2)));                                              %     2/3/2016
    end
    analysis.thickness.Mar=analysis.thickness.Ir_L_Th/days2_1;                                         % Mineral Apposition Rate (Inter Label Thickness / days)
else
    analysis.thickness.Ir_L_Th=NaN;
    analysis.thickness.Mar=NaN;
    Ir_L_Th_std=0;
end

eval(['analysis.volume.TV=tissue_volume',num2str(tr_loop),';'])                 % Tissue Volume (in pixel area) Dec. 08, 2009
analysis.volume.BV=analysis.area.DIC;                                           % Bone Volume (in pixel area)
analysis.volume.BV_TV=analysis.volume.BV/analysis.volume.TV*100;                % Bone volume / Tissue Volume (area ratio)

if analysis.thickness.Tb_Th~=0
    analysis.thickness.Tb_N=analysis.volume.BV_TV/100/analysis.thickness.Tb_Th;
    analysis.thickness.Tb_Sp=(1/analysis.thickness.Tb_N)-analysis.thickness.Tb_Th;
else
    analysis.thickness.Tb_N=0;
    analysis.thickness.Tb_Sp=NaN;
end
analysis.thickness.BFR=analysis.thickness.Mar*analysis.surface.MS_BS/100;


rgT=r_label & g_label & logical(TRAP_label);
[analysis.perimeter.rgT, analysis.area.rgT] = perimeter_area1(test_DIC_temp, rgT, dis_10x_for_400micron_per_312pixels, exact_boundary);
rT=r_label & logical(TRAP_label); rT=rT-rgT;
[analysis.perimeter.rT, analysis.area.rT] = perimeter_area1(test_DIC_temp, rT, dis_10x_for_400micron_per_312pixels, exact_boundary);
gT=g_label & logical(TRAP_label); gT=gT-rgT;
[analysis.perimeter.gT, analysis.area.gT] = perimeter_area1(test_DIC_temp, gT, dis_10x_for_400micron_per_312pixels, exact_boundary);
T_only=logical(TRAP_label)-(rgT+rT+gT);
T_only1=logical(TRAP_label)-(rT);
[analysis.perimeter.T_only, analysis.area.T_only] = perimeter_area1(test_DIC_temp, T_only, dis_10x_for_400micron_per_312pixels, exact_boundary);
[analysis.perimeter.T_only1, analysis.area.T_only1] = perimeter_area1(test_DIC_temp, T_only1, dis_10x_for_400micron_per_312pixels, exact_boundary);
%analysis.perimeter.T=analysis.perimeter.rgT+analysis.perimeter.rT+analysis.perimeter.gT+analysis.perimeter.T_only;

analysis.surface.T_BS=analysis.perimeter.T/analysis.surface.BS*100;
analysis.surface.T_only_BS=analysis.perimeter.T_only/analysis.surface.BS*100;
analysis.surface.T_only1_BS=analysis.perimeter.T_only1/analysis.surface.BS*100;
analysis.surface.rT_BS=analysis.perimeter.rT/analysis.surface.BS*100;
analysis.surface.gT_BS=analysis.perimeter.gT/analysis.surface.BS*100;
analysis.surface.rgT_BS=analysis.perimeter.rgT/analysis.surface.BS*100;

rgA=r_label & g_label & logical(AP_label);
[analysis.perimeter.rgA, analysis.area.rgA] = perimeter_area1(test_DIC_temp, rgA, dis_10x_for_400micron_per_312pixels, exact_boundary);
rA=r_label & logical(AP_label); rA=rA-rgA;
[analysis.perimeter.rA, analysis.area.rA] = perimeter_area1(test_DIC_temp, rA, dis_10x_for_400micron_per_312pixels, exact_boundary);
gA=g_label & logical(AP_label); gA=gA-rgA;
[analysis.perimeter.gA, analysis.area.gA] = perimeter_area1(test_DIC_temp, gA, dis_10x_for_400micron_per_312pixels, exact_boundary);
A_only=logical(AP_label)-(rgA+rA+gA);
A_only1=logical(AP_label)-(rA);
[analysis.perimeter.A_only, analysis.area.A_only] = perimeter_area1(test_DIC_temp, A_only, dis_10x_for_400micron_per_312pixels, exact_boundary);
[analysis.perimeter.A_only1, analysis.area.A_only1] = perimeter_area1(test_DIC_temp, A_only1, dis_10x_for_400micron_per_312pixels, exact_boundary);
%analysis.perimeter.A=analysis.perimeter.rgA+analysis.perimeter.rA+analysis.perimeter.gA+analysis.perimeter.A_only;

analysis.surface.A_BS=analysis.perimeter.A/analysis.surface.BS*100;
analysis.surface.A_only_BS=analysis.perimeter.A_only/analysis.surface.BS*100;
analysis.surface.A_only1_BS=analysis.perimeter.A_only1/analysis.surface.BS*100;
analysis.surface.rA_BS=analysis.perimeter.rA/analysis.surface.BS*100;
analysis.surface.gA_BS=analysis.perimeter.gA/analysis.surface.BS*100;
analysis.surface.rgA_BS=analysis.perimeter.rgA/analysis.surface.BS*100;

rgAT=r_label & g_label & logical(AP_label) & logical(TRAP_label);
[analysis.perimeter.rgAT, analysis.area.rgAT] = perimeter_area1(test_DIC_temp, rgAT, dis_10x_for_400micron_per_312pixels, exact_boundary);
rAT=r_label & logical(AP_label) & logical(TRAP_label);rAT=rAT-rgAT;
[analysis.perimeter.rAT, analysis.area.rAT] = perimeter_area1(test_DIC_temp, rAT, dis_10x_for_400micron_per_312pixels, exact_boundary);
gAT=g_label & logical(AP_label) & logical(TRAP_label);gAT=gAT-rgAT;
[analysis.perimeter.gAT, analysis.area.gAT] = perimeter_area1(test_DIC_temp, gAT, dis_10x_for_400micron_per_312pixels, exact_boundary);

AT_only=logical(AT)-(rAT+gAT+rgAT);
[analysis.perimeter.AT_only, analysis.area.AT_only] = perimeter_area1(test_DIC_temp, AT_only, dis_10x_for_400micron_per_312pixels, exact_boundary);
analysis.surface.rAT_BS=analysis.perimeter.rAT/analysis.surface.BS*100;
analysis.surface.gAT_BS=analysis.perimeter.gAT/analysis.surface.BS*100;
analysis.surface.rgAT_BS=analysis.perimeter.rgAT/analysis.surface.BS*100;
analysis.surface.AT_only_BS=analysis.perimeter.AT_only/analysis.surface.BS*100;

analysis.surface.AT_BS=analysis.perimeter.AT/analysis.surface.BS*100;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   GFP analysis similar to AP
%   Added February 07, 2013
%
rgGFP=r_label & g_label & logical(GFP_label);
[analysis.perimeter.rgGFP, analysis.area.rgGFP] = perimeter_area1(test_DIC_temp, rgGFP, dis_10x_for_400micron_per_312pixels, exact_boundary);
rGFP=r_label & logical(GFP_label); rGFP=rGFP-rgGFP;
[analysis.perimeter.rGFP, analysis.area.rGFP] = perimeter_area1(test_DIC_temp, rGFP, dis_10x_for_400micron_per_312pixels, exact_boundary);
gGFP=g_label & logical(GFP_label); gGFP=gGFP-rgGFP;
[analysis.perimeter.gGFP, analysis.area.gGFP] = perimeter_area1(test_DIC_temp, gGFP, dis_10x_for_400micron_per_312pixels, exact_boundary);
GFP_only=logical(GFP_label)-(rgGFP+rGFP+gGFP);
GFP_only1=logical(GFP_label)-(rGFP);
[analysis.perimeter.GFP_only, analysis.area.GFP_only] = perimeter_area1(test_DIC_temp, GFP_only, dis_10x_for_400micron_per_312pixels, exact_boundary);
[analysis.perimeter.GFP_only1, analysis.area.GFP_only1] = perimeter_area1(test_DIC_temp, GFP_only1, dis_10x_for_400micron_per_312pixels, exact_boundary);
%analysis.perimeter.GFP=analysis.perimeter.rgGFP+analysis.perimeter.rGFP+analysis.perimeter.gGFP+analysis.perimeter.GFP_only;

analysis.surface.GFP_BS=analysis.perimeter.GFP/analysis.surface.BS*100;
analysis.surface.GFP_only_BS=analysis.perimeter.GFP_only/analysis.surface.BS*100;
analysis.surface.GFP_only1_BS=analysis.perimeter.GFP_only1/analysis.surface.BS*100;
analysis.surface.rGFP_BS=analysis.perimeter.rGFP/analysis.surface.BS*100;
analysis.surface.gGFP_BS=analysis.perimeter.gGFP/analysis.surface.BS*100;
analysis.surface.rgGFP_BS=analysis.perimeter.rgGFP/analysis.surface.BS*100;

rgGFPT=r_label & g_label & logical(GFP_label) & logical(TRAP_label);
[analysis.perimeter.rgGFPT, analysis.area.rgGFPT] = perimeter_area1(test_DIC_temp, rgGFPT, dis_10x_for_400micron_per_312pixels, exact_boundary);
rGFPT=r_label & logical(GFP_label) & logical(TRAP_label);rGFPT=rGFPT-rgGFPT;
[analysis.perimeter.rGFPT, analysis.area.rGFPT] = perimeter_area1(test_DIC_temp, rGFPT, dis_10x_for_400micron_per_312pixels, exact_boundary);
gGFPT=g_label & logical(GFP_label) & logical(TRAP_label);gGFPT=gGFPT-rgGFPT;
[analysis.perimeter.gGFPT, analysis.area.gGFPT] = perimeter_area1(test_DIC_temp, gGFPT, dis_10x_for_400micron_per_312pixels, exact_boundary);

GFPT=logical(GFP_label) & logical(TRAP_label);
[analysis.perimeter.GFPT, analysis.area.GFPT] = perimeter_area1(test_DIC_temp, GFPT, dis_10x_for_400micron_per_312pixels, exact_boundary);
GFPT_only=logical(GFPT)-(rGFPT+gGFPT+rgGFPT);
[analysis.perimeter.GFPT_only, analysis.area.GFPT_only] = perimeter_area1(test_DIC_temp, GFPT_only, dis_10x_for_400micron_per_312pixels, exact_boundary);
analysis.surface.rGFPT_BS=analysis.perimeter.rGFPT/analysis.surface.BS*100;
analysis.surface.gGFPT_BS=analysis.perimeter.gGFPT/analysis.surface.BS*100;
analysis.surface.rgGFPT_BS=analysis.perimeter.rgGFPT/analysis.surface.BS*100;
analysis.surface.GFPT_only_BS=analysis.perimeter.GFPT_only/analysis.surface.BS*100;

analysis.surface.GFPT_BS=analysis.perimeter.GFPT/analysis.surface.BS*100;
%
%   Added February 07, 2013
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%


t_DIC=imdilate(test_DIC, strel('disk',20));
test=zeros(size(t_DIC));
TRAP_on_surface=test;
TRAP_off_surface=test;
[L n]=bwlabel(test_TRAP);
stats=regionprops(L, 'BoundingBox','Image'); clear L
for search_loop=1:n
    sx=round(stats(search_loop).BoundingBox(2));
    sy=round(stats(search_loop).BoundingBox(1));
    ex=sx+stats(search_loop).BoundingBox(4)-1;
    ey=sy+stats(search_loop).BoundingBox(3)-1;
    test(sx:ex,sy:ey)=stats(search_loop).Image;
    if size(find(t_DIC & test),1)>0
        TRAP_on_surface=TRAP_on_surface+test;
    else
        TRAP_off_surface=TRAP_off_surface+test;
    end
    test=zeros(size(t_DIC));
end
clear t_DIC test
analysis.volume.T_on=size(find(TRAP_on_surface),1)*dis_10x_for_400micron_per_312pixels^2;
analysis.volume.T_off=size(find(TRAP_off_surface),1)*dis_10x_for_400micron_per_312pixels^2;
analysis.volume.T_on_TV=analysis.volume.T_on/analysis.volume.TV*100;
analysis.volume.T_off_TV=analysis.volume.T_off/analysis.volume.TV*100;
analysis.volume.T_on_T=analysis.volume.T_on/(analysis.volume.T_on+analysis.volume.T_off)*100;
analysis.volume.T_off_T=100-analysis.volume.T_on_T;


if write==1
    if with_rg_label==1
        if tr_loop==1
            eval(['save ''',direc_ana,c1,'_analysis_tr1.mat'' analysis']);
        elseif tr_loop==2
            eval(['save ''',direc_ana,c1,'_analysis_tr2.mat'' analysis']);
        end
    elseif with_rg_label==0
        if tr_loop==1
            eval(['save ''',direc_ana,c1,'_analysis_tr1_f.mat'' analysis']);
        elseif tr_loop==2
            eval(['save ''',direc_ana,c1,'_analysis_tr2_f.mat'' analysis']);
        end
    end
end

%if write==1
%    if exp_type=='o'
%        if with_rg_label==1
%            eval(['save ''',direc,b_t,b_t,b_n,exp_t,b_t1,b_n,common,'.jpg_Files\analysis'' analysis -V6'])
%        elseif with_rg_label==0
%            eval(['save
%            ''',direc,b_t,b_t,b_n,exp_t,b_t1,b_n,common,'.jpg_Files\analysis_f'' analysis -V6'])
%        end
%    else
%        if with_rg_label==1
%            eval(['save ''',direc,b_t,b_t,b_n,b_t1,b_n,common,'.jpg_Files\analysis'' analysis -V6'])
%        elseif woth_rg_label==0
%            eval(['save ''',direc,b_t,b_t,b_n,b_t1,b_n,common,'.jpg_Files\analysis_f'' analysis -V6'])
%        end
%    end
%    %eval(['save ''C:\Users\shhong\Desktop\CJake\Analysis s2
%end

%3\JW3_WT\JW3_WT_5FR_v2\JW3_WT_5FR_fix_c.jpg_Files\analysis'' analysis -V6']);
    
clear temp temp1
% temp1=imerode(image_in,strel('disk',2));
% temp=temp1 + imdilate(r_label,strel('disk',1)) + imdilate(TRAP_label,strel('disk',1)) + imdilate(AP_label*0.6,strel('disk',1));
% temp(:,:,2)=temp1 + imdilate(g_label,strel('disk',1)) + imdilate(TRAP_label, strel('disk',1)) + imdilate(GFP_label, strel('disk',1)) + imdilate(AP_label*0.3,strel('disk',1));
% temp(:,:,3)=temp1 + imdilate(AP_label*0.3,strel('disk',1)) + imdilate(GFP_label, strel('disk',1));
% temp=temp/max(temp(:));
se1=strel('disk',1);
se=strel('disk',2);
temp1=imerode(image_in,se1);
%%%%%%%%%%%%%%%%%%%%%
%
%	June, 7th 2018
%
if ~RG_label	% 1st label : red,  2nd label : green --> Reversed labels
	temp=temp1*0.5 + (imdilate(g_label,se) + imdilate(TRAP_label,se) + imdilate(AP_label,se)*0.6).*(1-temp1);
	temp(:,:,2)=temp1*0.5 + (imdilate(r_label,se) + imdilate(TRAP_label,se) + imdilate(GFP_label,se) + imdilate(AP_label,se)*0.3).*(1-temp1);
	temp(:,:,3)=temp1*0.5 + (imdilate(AP_label,se)*0.3 + imdilate(GFP_label,se)).*(1-temp1);
else			% 1st label : green,  2nd label : red --> Normal labels
	temp=temp1*0.5 + (imdilate(r_label,se) + imdilate(TRAP_label,se) + imdilate(AP_label,se)*0.6).*(1-temp1);
	temp(:,:,2)=temp1*0.5 + (imdilate(g_label,se) + imdilate(TRAP_label,se) + imdilate(GFP_label,se) + imdilate(AP_label,se)*0.3).*(1-temp1);
	temp(:,:,3)=temp1*0.5 + (imdilate(AP_label,se)*0.3 + imdilate(GFP_label,se)).*(1-temp1);
end
%
%	June, 7th 2018
%
%%%%%%%%%%%%%%%%%%%%%
temp=temp/max(temp(:));

% ttt=test_DIC-(test_DIC&r_label&g_label);
% temp(:,:,3)=immultiply(uint8(ttt),128);
% 
% %
% %       January 5, 2010
% %
% % if red_first_green_last==1
% %     if b_t==1
% %         temp(:,:,1)=immultiply(uint8(g_label),256);
% %         %temp(:,:,2)=immultiply(uint8(GFP_label),255);
% %         %temp(:,:,3)=imadd(temp(:,:,3),immultiply(uint8(GFP_label),255));
% %         temp(:,:,2)=imadd(temp(:,:,2),immultiply(uint8(r_label),255));
% %         %temp(:,:,1)=imadd(temp(:,:,1),immultiply(uint8(r_label),255));
% %     elseif b_t==2
% %         temp(:,:,1)=immultiply(uint8(r_label),256);
% %         %temp(:,:,2)=immultiply(uint8(GFP_label),255);
% %         %temp(:,:,3)=imadd(temp(:,:,3),immultiply(uint8(GFP_label),255));
% %         temp(:,:,2)=imadd(temp(:,:,2),immultiply(uint8(g_label),255));
% %         %temp(:,:,1)=imadd(temp(:,:,1),immultiply(uint8(r_label),255));
% %     end
% % else
%      temp(:,:,1)=immultiply(uint8(r_label),256);
%      %temp(:,:,2)=immultiply(uint8(GFP_label),255);
%      %temp(:,:,3)=imadd(temp(:,:,3),immultiply(uint8(GFP_label),255));
%      temp(:,:,2)=imadd(temp(:,:,2),immultiply(uint8(g_label),255));
%      %temp(:,:,1)=imadd(temp(:,:,1),immultiply(uint8(r_label),255));    
% % end
% %       January 5, 2010
% %


%if b_t==1
%    temp(:,:,1)=immultiply(uint8(g_label),256);
%    %temp(:,:,2)=immultiply(uint8(GFP_label),255);
%    %temp(:,:,3)=imadd(temp(:,:,3),immultiply(uint8(GFP_label),255));
%    temp(:,:,2)=imadd(temp(:,:,2),immultiply(uint8(r_label),255));
%    %temp(:,:,1)=imadd(temp(:,:,1),immultiply(uint8(r_label),255));
%elseif b_t==2
%    temp(:,:,1)=immultiply(uint8(r_label),256);
%    %temp(:,:,2)=immultiply(uint8(GFP_label),255);
%    %temp(:,:,3)=imadd(temp(:,:,3),immultiply(uint8(GFP_label),255));
%    temp(:,:,2)=imadd(temp(:,:,2),immultiply(uint8(g_label),255));
%    %temp(:,:,1)=imadd(temp(:,:,1),immultiply(uint8(r_label),255));
%end
% figure;imshow(temp(1:end-down_margin,:,:))

%if write==1
%    if tr_loop==1
%        eval(['imwrite(imrotate(temp(1:end-down_margin,:,:),-rot_angle), ''',direc,c1,'_trim11_tr1.jpg'',''jpg'',''quality'',100);']);
%    elseif tr_loop==2
%        eval(['imwrite(imrotate(temp(1:end-down_margin,:,:),-rot_angle), ''',direc,c1,'_trim11_tr2.jpg'',''jpg'',''quality'',100);']);
%    end
%end
if write==1
%    temp_image=imrotate(temp(1:end-down_margin,:,:),-rot_angle);
    temp_image=temp(1:end-down_margin,:,:);
    if tr_loop==1
        eval(['imwrite(temp_image, ''',direc,c1,'_trim11_tr1.jpg'',''jpg'',''quality'',100);']);
        eval(['imwrite(temp_image(',num2str(min(It)-50),':min(size(temp_image,1),',num2str(max(It)+50),'),',num2str(min(Jt)-50),':',num2str(max(Jt)+50),',:),''',direc,c1,'_trim11_tr1_inside.jpg'',''jpg'',''quality'',100);']);
        %eval(['imwrite(temp_image, ''',archive_d,delimeter,c1,'_trim11_tr1.jpg'',''jpg'',''quality'',100);']);
        %eval(['imwrite(temp_image(',num2str(min(It)-50),':min(size(temp_image,1),',num2str(max(It)+50),'),',num2str(min(Jt)-50),':',num2str(max(Jt)+50),',:),''',archive_d,delimeter,c1,'_trim11_tr1_inside.jpg'',''jpg'',''quality'',100);']);
    elseif tr_loop==2
        eval(['imwrite(temp_image, ''',direc,c1,'_trim11_tr2.jpg'',''jpg'',''quality'',100);']);
        eval(['imwrite(temp_image(',num2str(min(It)-50),':min(size(temp_image,1),',num2str(max(It)+50),'),',num2str(min(Jt)-50),':',num2str(max(Jt)+50),',:),''',direc,c1,'_trim11_tr2_inside.jpg'',''jpg'',''quality'',100);']);
        %eval(['imwrite(temp_image, ''',archive_d,delimeter,c1,'_trim11_tr2.jpg'',''jpg'',''quality'',100);']);
        %eval(['imwrite(temp_image(',num2str(min(It)-50),':min(size(temp_image,1),',num2str(max(It)+50),'),',num2str(min(Jt)-50),':',num2str(max(Jt)+50),',:),''',archive_d,delimeter,c1,'_trim11_tr2_inside.jpg'',''jpg'',''quality'',100);']);
    end
end
%if write==1
%    if exp_type=='o'
%        eval(['imwrite(imrotate(temp(1:end-down_margin,:,:),-rot_angle), ''',direc,b_t,b_t,b_n,exp_t,b_t1,b_n,common,'.jpg_Files',b_t1,b_n,common,'_trim11.jpg'',''jpg'',''quality'',100);'])
%    else
%        eval(['imwrite(imrotate(temp(1:end-down_margin,:,:),-rot_angle), ''',direc,b_t,b_t,b_n,b_t1,b_n,common,'.jpg_Files',b_t1,b_n,common,'_trim11.jpg'',''jpg'',''quality'',100);'])
%    end
%end


% temp(:,:,3)=immultiply(uint8(ttt),128);
% temp(:,:,1)=immultiply(uint8(r_label),128);
% %temp(:,:,2)=immultiply(uint8(GFP_label),128);
% %temp(:,:,3)=imadd(temp(:,:,3),immultiply(uint8(GFP_label),128));
% temp(:,:,2)=imadd(temp(:,:,2),immultiply(uint8(g_label),128));
% %temp(:,:,1)=imadd(temp(:,:,1),immultiply(uint8(r_label),128));
% %figure;imshow(temp(1:end-down_margin,:,:))


g_trim_new=zeros(size(g_trim));
g_trim2=zeros(size(g_trim));
if no_green==1 | no_red==1
    g_trim_new=g_trim;
else

n=whos;

% if isempty(regexp([n.name],'TB'))==0        %do this routine when TB exists
if isfield(TB,'G_label')        %do this routine when TB exists
for i=1:length(TB)
    %i

    sdx=round(stats_DIC(J(i)).BoundingBox(4)/2);
    sdy=round(stats_DIC(J(i)).BoundingBox(3)/2);
    edx=round(stats_DIC(J(i)).BoundingBox(4)/2);
    edy=round(stats_DIC(J(i)).BoundingBox(3)/2);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   December, 22 2009
%   
%
%    start_x=round(stats_DIC(J(i)).BoundingBox(2))-sdx;
%     if start_x<1
%         sdx=sdx-start_x;
%         start_x=1;
%     end
%     start_y=round(stats_DIC(J(i)).BoundingBox(1))-sdy;
%     if start_y<1
%         sdy=sdy-start_y;
%         start_y=1;
%     end
%     end_x=start_x+stats_DIC(J(i)).BoundingBox(4)-1+sdx+edx;
%     if end_x>size(test1,1)
%         edx=edx-(end_x-size(test1,1));
%         end_x=size(test1,1);
%     end
%     end_y=start_y+stats_DIC(J(i)).BoundingBox(3)-1+sdy+edy;
%     if end_y>size(test1,2)
%         edy=edy-(end_y-size(test1,2));
%         end_y=size(test1,2);
%     end
%     im=zeros(end_x-start_x+1,end_y-start_y+1);
%     im(sdx+1:end-edx,sdy+1:end-edy)=stats_DIC(J(i)).FilledImage;
%     g_test=test1(start_x:end_x,start_y:end_y)&im;
%     r_test=test2(start_x:end_x,start_y:end_y)&im;

    sx=round(stats_DIC(J(i)).BoundingBox(2));
    ex=sx+stats_DIC(J(i)).BoundingBox(4)-1;
    sy=round(stats_DIC(J(i)).BoundingBox(1));
    ey=sy+stats_DIC(J(i)).BoundingBox(3)-1;
    
    im=zeros(stats_DIC(J(i)).BoundingBox(4)+2*sdx,stats_DIC(J(i)).BoundingBox(3)+2*sdy);
    g_test=im;
    r_test=im;
    im(sdx+1:sdx+stats_DIC(J(i)).BoundingBox(4),sdy+1:sdy+stats_DIC(J(i)).BoundingBox(3))=stats_DIC(J(i)).FilledImage;
    g_test(sdx+1:sdx+stats_DIC(J(i)).BoundingBox(4),sdy+1:sdy+stats_DIC(J(i)).BoundingBox(3))=test1(sx:ex,sy:ey) & stats_DIC(J(i)).FilledImage;
    r_test(sdx+1:sdx+stats_DIC(J(i)).BoundingBox(4),sdy+1:sdy+stats_DIC(J(i)).BoundingBox(3))=test2(sx:ex,sy:ey) & stats_DIC(J(i)).FilledImage;
    
%   December, 22 2009
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    
    if red_first_green_last==1 & b_t==1
        r_t=r_test;                                     % Jan. 05 2010
        r_test=g_test;                                  % Jan. 05 2010
        g_test=r_t;                                     % Jan. 05 2010
    end
        
    red_im=r_test;
    
    if isempty(find(g_test)==1)==0 && isempty(find(r_test)==1)==0        % green and red in DIC
        g_trim_temp=zeros(size(g_test));
        for j=1:length(TB(i).G_label)
            
            if isempty(TB(i).G_label(j).trim)==0
                if length(TB(i).G_label(j).trim.col)>0
% 12/22/09                    [g_trim_temp]=move1_2(TB(i).G_label(j), r_trim1(start_x:end_x,start_y:end_y), analysis.thickness.Ir_L_Th+Ir_L_Th_std*std_const);       % 1/23/2008
                    if red_first_green_last==1 & b_t==1
% 06/-3/2010                        r_test(sdx+1:sdx+stats_DIC(J(i)).BoundingBox(4),sdy+1:sdy+stats_DIC(J(i)).BoundingBox(3))=g_trim1(sx:ex,sy:ey) & stats_DIC(J(i)).FilledImage;  % Jan. 05 2010
                        r_test(sdx+1:sdx+stats_DIC(J(i)).BoundingBox(4),sdy+1:sdy+stats_DIC(J(i)).BoundingBox(3))=g_trim1(sx:ex,sy:ey) & stats_DIC(J(i)).Image;  % Jan. 05 2010
                    else
% 06/-3/2010                        r_test(sdx+1:sdx+stats_DIC(J(i)).BoundingBox(4),sdy+1:sdy+stats_DIC(J(i)).BoundingBox(3))=r_trim1(sx:ex,sy:ey) & stats_DIC(J(i)).FilledImage;
                        r_test(sdx+1:sdx+stats_DIC(J(i)).BoundingBox(4),sdy+1:sdy+stats_DIC(J(i)).BoundingBox(3))=r_trim1(sx:ex,sy:ey) & stats_DIC(J(i)).Image;
                    end
                    [g_trim_temp]=move1_2(TB(i).G_label(j), r_test, analysis.thickness.Ir_L_Th+Ir_L_Th_std*std_const);       % 1/23/2008
                    
                    tmp=double(g_trim_temp>0);
                    g_trim_new(sx:ex,sy:ey)=g_trim_new(sx:ex,sy:ey)+tmp(sdx+1:sdx+stats_DIC(J(i)).BoundingBox(4),sdy+1:sdy+stats_DIC(J(i)).BoundingBox(3));
% 12/22/09                    g_trim_new(start_x:end_x,start_y:end_y)=g_trim_new(start_x:end_x,start_y:end_y)+tmp;
                    % 1/7/08                g_trim1(start_x:end_x,start_y:end_y)=g_trim1(start_x:end_x,start_y:end_y)+g_trim_temp;
                    [g_trim2(sx:ex,sy:ey)]=find_min(g_trim2(sx:ex,sy:ey), g_trim_temp(sdx+1:sdx+stats_DIC(J(i)).BoundingBox(4),sdy+1:sdy+stats_DIC(J(i)).BoundingBox(3)));
% 12/22/09                    [g_trim2(start_x:end_x,start_y:end_y)]=find_min(g_trim2(start_x:end_x,start_y:end_y), g_trim_temp);
                end
            end
        end
% 1/3/08        g_trim_wo_red(start_x:end_x,start_y:end_y)=g_trim_wo_red(start_x:end_x,start_y:end_y)+g_trim_wo_red_temp;
        
        if red_first_green_last==1 & b_t==1
            r_trim_new=g_trim_new;                                  % Jan. 05 2010
            r_trim2=g_trim2;                                        % Jan. 05 2010
        end
            
    elseif isempty(find(g_test)==1)==0 && isempty(find(r_test)==1)==1        % green, but no red in DIC
        for j=1:length(TB(i).G_label)
            %start_r_x=round(stats_g(j).BoundingBox(2));
            %end_r_x=start_r_x+stats_g(j).BoundingBox(4)-1;
            %start_r_y=round(stats_g(j).BoundingBox(1));
            %end_r_y=start_r_y+stats_g(j).BoundingBox(3)-1;
            if isempty(TB(i).G_label(j).trim)==0
                if length(TB(i).G_label(j).trim.col)>0
                    %               [g_trim_wo_red_temp]=move_4(im,stats_g(j));
                    [g_trim_wo_red_temp]=move1_2(TB(i).G_label(j));       % 1/23/2008
                    
                    tmp=double(g_trim_wo_red_temp>0);
                    g_trim_new(sx:ex,sy:ey)=g_trim_new(sx:ex,sy:ey)+tmp(sdx+1:sdx+stats_DIC(J(i)).BoundingBox(4),sdy+1:sdy+stats_DIC(J(i)).BoundingBox(3));
% 12/22/09                    g_trim_new(start_x:end_x,start_y:end_y)=g_trim_new(start_x:end_x,start_y:end_y)+tmp;
                    % 1/7/08                g_trim1(start_x:end_x,start_y:end_y)=g_trim1(start_x:end_x,start_y:end_y)+g_trim_wo_red_temp;
                    [g_trim2(sx:ex,sy:ey)]=find_min(g_trim2(sx:ex,sy:ey), g_trim_wo_red_temp(sdx+1:sdx+stats_DIC(J(i)).BoundingBox(4),sdy+1:sdy+stats_DIC(J(i)).BoundingBox(3)));
% 12/22/09                    [g_trim2(start_x:end_x,start_y:end_y)]=find_min(g_trim2(start_x:end_x,start_y:end_y), g_trim_wo_red_temp);
                    %g_trim_wo_red(start_x:end_x,start_y:end_y)=g_trim_wo_red(start_x:end_x,start_y:end_y)+g_test;%(start_x:end_x,start_y:end_y);
                end
            end
        end
        if red_first_green_last==1 & b_t==1
            r_trim_new=g_trim_new;                                  % Jan. 05 2010
            r_trim2=g_trim2;                                        % Jan. 05 2010
        end
    end
end     % for i
end % for isfield
% end     % for isempty
end % if no_green
toc
analysis_n.intensity=analysis.intensity;
analysis_n.roi=analysis.roi;
analysis_n.cortex=analysis.cortex;

analysis_n.perimeter.DIC=analysis.perimeter.DIC;
%analysis_n.area.DIC=analysis.area.DIC;
analysis_n.area.DIC=analysis.volume.BV;
analysis_n.perimeter.GFP=analysis.perimeter.GFP;
analysis_n.area.GFP=analysis.area.GFP;
analysis_n.perimeter.AGFP=analysis.perimeter.AGFP;
analysis_n.perimeter.TGFP=analysis.perimeter.TGFP;
analysis_n.area.AGFP=analysis.area.AGFP;
analysis_n.area.TGFP=analysis.area.TGFP;
%clear g g2 g2t r  result

%[out_label]=outside_trim_label(image_in, g_trim_wo_red, 20, exact_boundary);
%[in_label]=inside_trim_label(image_in, g_trim_wo_red, loop_count, exact_boundary);
%g_trim_wo_red_label=out_label;% | in_label; clear in_label out_label
%g_label=double(in_label)+g_trim;

%r_label=r_trim;
%r_label=imdilate(r_label,strel('disk',1));
%test_DIC=logical(double(test_DIC)/255);
%r_label=r_label & (test_DIC-imerode(test_DIC,strel('disk',1)));
%[analysis_n.perimeter.r, analysis_n.area.r] = perimeter_area1(test_DIC_temp, r_label);
analysis_n.perimeter.r=analysis.perimeter.r;
analysis_n.area.r=analysis.area.r;

if red_first_green_last==1 & b_t==1
    r_label=r_trim_new;                                                 % Jan. 05 2010
    r_label=uint8(logical(r_label));                                    % Jan. 05 2010
    r_label=imdilate(r_label,strel('disk',1));                          % Jan. 05 2010
    r_label=r_label & (test_DIC-imerode(test_DIC,strel('disk',1)));     % Jan. 05 2010
else
    g_label=g_trim_new;
    g_label=uint8(logical(g_label));
    g_label=imdilate(g_label,strel('disk',1));
    g_label=g_label & (test_DIC-imerode(test_DIC,strel('disk',1)));
end

[analysis_n.perimeter.g, analysis_n.area.g] = perimeter_area1(test_DIC_temp, logical(g_label), dis_10x_for_400micron_per_312pixels, exact_boundary);

rg=r_label & logical(g_label);
[analysis_n.perimeter.rg, analysis_n.area.rg] = perimeter_area1(test_DIC_temp, rg, dis_10x_for_400micron_per_312pixels, exact_boundary);

%r_only=r_label & (1-rg);
%[analysis.perimeter.r_only, analysis.area.r_only] = perimeter_area1(test_DIC_temp, r_only);
analysis_n.perimeter.r_only=analysis_n.perimeter.r-analysis_n.perimeter.rg;
analysis_n.area.r_only=analysis_n.area.r-analysis_n.area.rg;
%g_only=g_label & (1-rg);
%[analysis.perimeter.g_only, analysis.area.g_only] = perimeter_area1(test_DIC_temp, g_only);
analysis_n.perimeter.g_only=analysis_n.perimeter.g-analysis_n.perimeter.rg;
analysis_n.area.g_only=analysis_n.area.g-analysis_n.area.rg;
%sLS=r_only | g_only ;
%[analysis.perimeter.sLS, analysis.area.sLS] = perimeter_area1(test_DIC_temp, sLS);
analysis_n.perimeter.sLS=analysis_n.perimeter.r_only+analysis_n.perimeter.g_only;
analysis_n.area.sLS=analysis_n.area.r_only+analysis_n.area.g_only;
%analysis.perimeter.LS=analysis.perimeter.sLS+analysis.perimeter.rg

rGFP=r_label & GFP_label;
[analysis_n.perimeter.rGFP, analysis_n.area.rGFP] = perimeter_area1(test_DIC_temp, rGFP, dis_10x_for_400micron_per_312pixels, exact_boundary);
gGFP=logical(g_label) & GFP_label;
[analysis_n.perimeter.gGFP, analysis_n.area.gGFP] = perimeter_area1(test_DIC_temp, gGFP, dis_10x_for_400micron_per_312pixels, exact_boundary);
rgGFP=rGFP & gGFP;
[analysis_n.perimeter.rgGFP, analysis_n.area.rgGFP] = perimeter_area1(test_DIC_temp, rgGFP, dis_10x_for_400micron_per_312pixels, exact_boundary);
analysis_n.surface.BS=analysis_n.perimeter.DIC;                                           % Bone Surface
analysis_n.surface.LS=analysis_n.perimeter.r+analysis_n.perimeter.g-analysis_n.perimeter.rg;                    % Labeled surface
analysis_n.surface.LS_BS=analysis_n.surface.LS/analysis_n.surface.BS*100;
analysis_n.surface.r_BS=analysis_n.perimeter.r/analysis_n.surface.BS*100;
analysis_n.surface.g_BS=analysis_n.perimeter.g/analysis_n.surface.BS*100;
analysis_n.surface.rg_BS=analysis_n.perimeter.rg/analysis_n.surface.BS*100;                                  % dLS/BS (double labeled surface / bone surface)
analysis_n.surface.rGFP_BS=analysis_n.perimeter.rGFP/analysis_n.surface.BS*100;
analysis_n.surface.gGFP_BS=analysis_n.perimeter.gGFP/analysis_n.surface.BS*100;
analysis_n.surface.rgGFP_BS=analysis_n.perimeter.rgGFP/analysis_n.surface.BS*100;
analysis_n.surface.GFP_BS=analysis_n.perimeter.GFP/analysis_n.surface.BS*100;                                 % GFP/BS
analysis_n.surface.AGFP_BS=analysis_n.perimeter.AGFP/analysis_n.surface.BS*100;                                 % GFP/BS
analysis_n.surface.TGFP_BS=analysis_n.perimeter.TGFP/analysis_n.surface.BS*100;                                 % GFP/BS
analysis_n.surface.r_only_BS=analysis_n.perimeter.r_only/analysis_n.surface.BS*100;   % sLS(red only)/BS (red only labeled surface / bone surface)
analysis_n.surface.g_only_BS=analysis_n.perimeter.g_only/analysis_n.surface.BS*100;   % sLS(green only)/BS (green only labeled surface / bone surface)
analysis_n.surface.sLS_BS=analysis_n.perimeter.sLS/analysis_n.surface.BS*100;   % sLS/BS (single labeled surface / bone surface)
analysis_n.surface.sLS_LS=analysis_n.perimeter.sLS/analysis_n.surface.LS*100;   % sLS/LS (single labeled surface / labeled surface)
if analysis.surface.LS~=0
    analysis_n.surface.sLS_LS=analysis_n.perimeter.sLS/analysis_n.surface.LS*100;   % sLS/LS (single labeled surface / labeled surface)
    analysis_n.surface.rg_LS=analysis_n.perimeter.rg/analysis_n.surface.LS*100;                                   % dLS/LS (double labeled surface / labeled surface)
elseif analysis.surface.LS==0
    analysis_n.surface.sLS_LS=0;   % sLS/LS (single labeled surface / labeled surface)
    analysis_n.surface.rg_LS=0;                                   % dLS/LS (double labeled surface / labeled surface)
elseif isnan(analysis.surface.LS)
    analysis_n.surface.sLS_LS=NaN;   % sLS/LS (single labeled surface / labeled surface)
    analysis_n.surface.rg_LS=NaN;                                   % dLS/LS (double labeled surface / labeled surface)
end
% analysis_n.surface.rg_LS=analysis_n.perimeter.rg/analysis_n.surface.LS*100;                                   % dLS/LS (double labeled surface / labeled surface)
% analysis_n.surface.LS_BS=analysis_n.surface.LS/analysis_n.surface.BS*100;
analysis_n.surface.dLS_sLS=analysis_n.perimeter.rg/analysis_n.perimeter.sLS*100;   % sLS/dLS (single labeled surface / double labeled surface)
analysis_n.surface.MS_BS=(analysis_n.surface.r_only_BS+analysis_n.surface.g_only_BS)/2+analysis_n.surface.rg_BS;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   modify after middle point is introduced
%   March 2, 2013
%   
% if red_first_green_last==1 & b_t==1
%     t=g_trim1 & r_trim2;
%     tt=t.*r_trim2-t.*g_trim1;                   % Jan. 5, 2010
% else
%     t=g_trim2 & r_trim1;
%     tt=t.*g_trim2-t.*r_trim1;
% end
% [dI dJ dV]=find(tt);
% sdV=sort(dV);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   New MAR using only double label
%   July 14th, 2015
%
% if red_first_green_last==1 & b_t==1
%     tt=r_trim1-g_trim1;               % Jan. 5, 2010
% else
%     tt=g_trim1-r_trim1;
% end
trim_template=(g_trim1>0)&(r_trim1>0);      % double label
if red_first_green_last==1
    tt=(r_trim1-g_trim1).*double(trim_template);
else
    tt=(g_trim1-r_trim1).*double(trim_template);
end
%
%   New MAR using only double_label
%   July 14th, 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[dI dJ dV]=find(tt);
sdV=sort(dV);%*400/dis_10x;
%
%   modify after middle point is introduced
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%isdV=find(sdV<0);
%analysis_n.thickness.Ir_L_Th=mean(sdV(length(isdV)+1:end));          % Inter Label Thickness
%analysis_n.thickness.Ir_L_Th=median(sdV(length(isdV)+1:end));          % Inter Label Thickness
%test_dv=sdV(length(isdV)+1:end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   June 3, 2010 : ignore the part after the big jump of test_dv1
%test_dv=sdV(find(sdV>0));
test_dv1=sdV(find(sdV>0));
%end_dv=find(diff(test_dv1)>mean(diff(test_dv1))*5,1);
end_dv=find(diff(test_dv1)>5,1);
if isempty(end_dv)==1
    test_dv=test_dv1*dis_10x_for_400micron_per_312pixels;
else
    test_dv=test_dv1(1:end_dv)*dis_10x_for_400micron_per_312pixels;   % added
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   modify after middle point is introduced
%   March 2, 2013
%   
% middle=round(length(test_dv)/2);
[x y]=find(test_dv1==mode(test_dv1));
if ~isempty(x)
    middle=min(x(round(length(x)/2)), round(length(test_dv)/2));
%
%   modify after middle point is introduced
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if middle<4
    %    analysis_n.thickness.Ir_L_Th=mean(sdV(length(isdV)+1:end));
        analysis_n.thickness.Ir_L_Th=mean(test_dv);
        Ir_L_Th_std=std(test_dv);
    else
    %    analysis_n.thickness.Ir_L_Th=mean(sdV(middle-floor(middle/2):middle+ceil(middle/2)));
    %     analysis_n.thickness.Ir_L_Th=mean(test_dv(middle-floor(middle/2):middle+ceil(middle/2)));
        analysis_n.thickness.Ir_L_Th=mean(test_dv(middle-floor(middle/2):middle+ceil(middle/4)));
        Ir_L_Th_std=std(test_dv(middle-floor(middle/2):middle+ceil(middle/4)));
    end
    analysis_n.thickness.Mar=analysis_n.thickness.Ir_L_Th/days2_1;                                         % Mineral Apposition Rate (Inter Label Thickness / days)
    %MAR_std=std(sdV(length(isdV)+1:end));
else
    analysis_n.thickness.Ir_L_Th=NaN;
    analysis_n.thickness.Mar=NaN;
    Ir_L_Th_std=0;
end
analysis_n.thickness.Tb_Th=analysis.thickness.Tb_Th;


%analysis_n.volume.TV=tissue_volume;                                           % Tissue Volume (in pixel area)
eval(['analysis_n.volume.TV=tissue_volume',num2str(tr_loop),';'])                 % Tissue Volume (in pixel area) Dec. 08, 2009
analysis_n.volume.BV=analysis_n.area.DIC;                                                % Bone Volume (in pixel area)
analysis_n.volume.BV_TV=analysis_n.volume.BV/analysis_n.volume.TV*100;                                             % Bone volume / Tissue Volume (area ratio)

if analysis_n.thickness.Tb_Th~=0
    analysis_n.thickness.Tb_N=analysis_n.volume.BV_TV/100/analysis_n.thickness.Tb_Th;
    analysis_n.thickness.Tb_Sp=(1/analysis_n.thickness.Tb_N)-analysis_n.thickness.Tb_Th;
else
    analysis_n.thickness.Tb_N=0;
    analysis_n.thickness.Tb_Sp=NaN;
end
% analysis_n.thickness.Tb_N=analysis_n.volume.BV_TV/100/analysis_n.thickness.Tb_Th;
% analysis_n.thickness.Tb_Sp=(1/analysis_n.thickness.Tb_N)-analysis_n.thickness.Tb_Th;
analysis_n.thickness.BFR=analysis_n.thickness.Mar*analysis_n.surface.MS_BS/100;
analysis.thickness

analysis_n.perimeter.T=analysis.perimeter.T;
rgT=r_label & g_label & logical(TRAP_label);
[analysis_n.perimeter.rgT, analysis_n.area.rgT] = perimeter_area1(test_DIC_temp, rgT, dis_10x_for_400micron_per_312pixels, exact_boundary);
rT=r_label & logical(TRAP_label); rT=rT-rgT;
[analysis_n.perimeter.rT, analysis_n.area.rT] = perimeter_area1(test_DIC_temp, rT, dis_10x_for_400micron_per_312pixels, exact_boundary);
gT=g_label & logical(TRAP_label); gT=gT-rgT;
[analysis_n.perimeter.gT, analysis_n.area.gT] = perimeter_area1(test_DIC_temp, gT, dis_10x_for_400micron_per_312pixels, exact_boundary);
T_only=logical(TRAP_label)-(rgT+rT+gT);
[analysis_n.perimeter.T_only, analysis_n.area.T_only] = perimeter_area1(test_DIC_temp, T_only, dis_10x_for_400micron_per_312pixels, exact_boundary);
%analysis.perimeter.T=analysis.perimeter.rgT+analysis.perimeter.rT+analysis.perimeter.gT+analysis.perimeter.T_only;

analysis_n.surface.T_BS=analysis_n.perimeter.T/analysis_n.surface.BS*100;
analysis_n.surface.T_only_BS=analysis_n.perimeter.T_only/analysis_n.surface.BS*100;
analysis_n.surface.rT_BS=analysis_n.perimeter.rT/analysis_n.surface.BS*100;
analysis_n.surface.gT_BS=analysis_n.perimeter.gT/analysis_n.surface.BS*100;
analysis_n.surface.rgT_BS=analysis_n.perimeter.rgT/analysis_n.surface.BS*100;

rgA=r_label & g_label & logical(AP_label);
[analysis_n.perimeter.rgA, analysis_n.area.rgA] = perimeter_area1(test_DIC_temp, rgA, dis_10x_for_400micron_per_312pixels, exact_boundary);
rA=r_label & logical(AP_label); rA=rA-rgA;
[analysis_n.perimeter.rA, analysis_n.area.rA] = perimeter_area1(test_DIC_temp, rA, dis_10x_for_400micron_per_312pixels, exact_boundary);
gA=g_label & logical(AP_label); gA=gA-rgA;
[analysis_n.perimeter.gA, analysis_n.area.gA] = perimeter_area1(test_DIC_temp, gA, dis_10x_for_400micron_per_312pixels, exact_boundary);
A_only=logical(AP_label)-(rgA+rA+gA);
[analysis_n.perimeter.A_only, analysis_n.area.A_only] = perimeter_area1(test_DIC_temp, A_only, dis_10x_for_400micron_per_312pixels, exact_boundary);
%analysis.perimeter.A=analysis.perimeter.rgA+analysis.perimeter.rA+analysis.perimeter.gA+analysis.perimeter.A_only;

analysis_n.perimeter.A=analysis.perimeter.A;
analysis_n.surface.A_BS=analysis_n.perimeter.A/analysis_n.surface.BS*100;
analysis_n.surface.A_only_BS=analysis_n.perimeter.A_only/analysis_n.surface.BS*100;
analysis_n.surface.rA_BS=analysis_n.perimeter.rA/analysis_n.surface.BS*100;
analysis_n.surface.gA_BS=analysis_n.perimeter.gA/analysis_n.surface.BS*100;
analysis_n.surface.rgA_BS=analysis_n.perimeter.rgA/analysis_n.surface.BS*100;

rgAT=r_label & g_label & logical(AP_label) & logical(TRAP_label);
[analysis_n.perimeter.rgAT, analysis_n.area.rgAT] = perimeter_area1(test_DIC_temp, rgAT, dis_10x_for_400micron_per_312pixels, exact_boundary);
rAT=r_label & logical(AP_label) & logical(TRAP_label); rAT=rAT-rgAT;
[analysis_n.perimeter.rAT, analysis_n.area.rAT] = perimeter_area1(test_DIC_temp, rAT, dis_10x_for_400micron_per_312pixels, exact_boundary);
gAT=g_label & logical(AP_label) & logical(TRAP_label); gAT=gAT-rgAT;
[analysis_n.perimeter.gAT, analysis_n.area.gAT] = perimeter_area1(test_DIC_temp, gAT, dis_10x_for_400micron_per_312pixels, exact_boundary);
AT_only=logical(AT)-(rgAT+rAT+gAT);
[analysis_n.perimeter.AT_only, analysis_n.area.AT_only] = perimeter_area1(test_DIC_temp, AT_only, dis_10x_for_400micron_per_312pixels, exact_boundary);

analysis_n.surface.rAT_BS=analysis_n.perimeter.rAT/analysis_n.surface.BS*100;
analysis_n.surface.gAT_BS=analysis_n.perimeter.gAT/analysis_n.surface.BS*100;
analysis_n.surface.rgAT_BS=analysis_n.perimeter.rgAT/analysis_n.surface.BS*100;

analysis_n.surface.AT_BS=analysis.surface.AT_BS;
analysis_n.perimeter.AT=analysis.perimeter.AT;
analysis_n.surface.AT_only_BS=analysis_n.perimeter.AT_only/analysis_n.surface.BS*100;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   GFP analysis similar to AP
%   Added February 07, 2013
%
rgGFP=r_label & g_label & logical(GFP_label);
[analysis_n.perimeter.rgGFP, analysis_n.area.rgGFP] = perimeter_area1(test_DIC_temp, rgGFP, dis_10x_for_400micron_per_312pixels, exact_boundary);
rGFP=r_label & logical(GFP_label); rGFP=rGFP-rgGFP;
[analysis_n.perimeter.rGFP, analysis_n.area.rGFP] = perimeter_area1(test_DIC_temp, rGFP, dis_10x_for_400micron_per_312pixels, exact_boundary);
gGFP=g_label & logical(GFP_label); gGFP=gGFP-rgGFP;
[analysis_n.perimeter.gGFP, analysis_n.area.gGFP] = perimeter_area1(test_DIC_temp, gGFP, dis_10x_for_400micron_per_312pixels, exact_boundary);
GFP_only=logical(GFP_label)-(rgGFP+rGFP+gGFP);
[analysis_n.perimeter.GFP_only, analysis_n.area.GFP_only] = perimeter_area1(test_DIC_temp, GFP_only, dis_10x_for_400micron_per_312pixels, exact_boundary);
%analysis.perimeter.A=analysis.perimeter.rgA+analysis.perimeter.rA+analysis.perimeter.gA+analysis.perimeter.A_only;

analysis_n.perimeter.GFP=analysis.perimeter.GFP;
analysis_n.surface.GFP_BS=analysis_n.perimeter.GFP/analysis_n.surface.BS*100;
analysis_n.surface.GFP_only_BS=analysis_n.perimeter.GFP_only/analysis_n.surface.BS*100;
analysis_n.surface.rGFP_BS=analysis_n.perimeter.rGFP/analysis_n.surface.BS*100;
analysis_n.surface.gGFP_BS=analysis_n.perimeter.gGFP/analysis_n.surface.BS*100;
analysis_n.surface.rgGFP_BS=analysis_n.perimeter.rgGFP/analysis_n.surface.BS*100;

rgGFPT=r_label & g_label & logical(GFP_label) & logical(TRAP_label);
[analysis_n.perimeter.rgGFPT, analysis_n.area.rgGFPT] = perimeter_area1(test_DIC_temp, rgGFPT, dis_10x_for_400micron_per_312pixels, exact_boundary);
rGFPT=r_label & logical(GFP_label) & logical(TRAP_label); rGFPT=rGFPT-rgGFPT;
[analysis_n.perimeter.rGFPT, analysis_n.area.rGFPT] = perimeter_area1(test_DIC_temp, rGFPT, dis_10x_for_400micron_per_312pixels, exact_boundary);
gGFPT=g_label & logical(GFP_label) & logical(TRAP_label); gGFPT=gGFPT-rgGFPT;
[analysis_n.perimeter.gGFPT, analysis_n.area.gGFPT] = perimeter_area1(test_DIC_temp, gGFPT, dis_10x_for_400micron_per_312pixels, exact_boundary);

GFPT=logical(GFP_label) & logical(TRAP_label);
[analysis_n.perimeter.GFPT, analysis_n.area.GFPT] = perimeter_area1(test_DIC_temp, GFPT, dis_10x_for_400micron_per_312pixels, exact_boundary);
GFPT_only=logical(GFPT)-(rgGFPT+rGFPT+gGFPT);
[analysis_n.perimeter.GFPT_only, analysis_n.area.GFPT_only] = perimeter_area1(test_DIC_temp, GFPT_only, dis_10x_for_400micron_per_312pixels, exact_boundary);

analysis_n.surface.rGFPT_BS=analysis_n.perimeter.rGFPT/analysis_n.surface.BS*100;
analysis_n.surface.gGFPT_BS=analysis_n.perimeter.gGFPT/analysis_n.surface.BS*100;
analysis_n.surface.rgGFPT_BS=analysis_n.perimeter.rgGFPT/analysis_n.surface.BS*100;

analysis_n.surface.GFPT_BS=analysis.surface.GFPT_BS;
analysis_n.perimeter.GFPT=analysis.perimeter.GFPT;
analysis_n.surface.GFPT_only_BS=analysis_n.perimeter.GFPT_only/analysis_n.surface.BS*100;
%
%   Added February 07, 2013
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%


analysis_n.surface.AGFP_BS=analysis.surface.AGFP_BS;
analysis_n.surface.TGFP_BS=analysis.surface.TGFP_BS;

analysis_n.volume.T_on=analysis.volume.T_on;
analysis_n.volume.T_off=analysis.volume.T_off;
analysis_n.volume.T_on_TV=analysis.volume.T_on_TV/analysis_n.volume.TV*100;
analysis_n.volume.T_off_TV=analysis.volume.T_off_TV/analysis_n.volume.TV*100;
analysis_n.volume.T_on_T=analysis.volume.T_on_T;
analysis_n.volume.T_off_T=analysis.volume.T_off_T;

if write==1
    if with_rg_label==1
        if tr_loop==1
            eval(['save ''',direc_ana,c1,'_analysis_n_tr1.mat'' analysis_n']);
        elseif tr_loop==2
            eval(['save ''',direc_ana,c1,'_analysis_n_tr2.mat'' analysis_n']);
        end
    elseif with_rg_label==0
        if tr_loop==1
            eval(['save ''',direc_ana,c1,'_analysis_n_f_tr1.mat'' analysis_n']);
        elseif tr_loop==2
            eval(['save ''',direc_ana,c1,'_analysis_n_f_tr2.mat'' analysis_n']);
        end
    end
end
%if write==1
%    if exp_type=='o'
%        if with_rg_label==1
%            eval(['save ''',direc,b_t,b_t,b_n,exp_t,b_t1,b_n,common,'.jpg_Files\analysis_n'' analysis_n -V6'])
%        elseif with_rg_label==0
%            eval(['save ''',direc,b_t,b_t,b_n,exp_t,b_t1,b_n,common,'.jpg_Files\analysis_n_f'' analysis -V6'])
%        end
%    else
%        if with_rg_label==1
%            eval(['save ''',direc,b_t,b_t,b_n,b_t1,b_n,common,'.jpg_Files\analysis_n'' analysis_n -V6'])
%        elseif woth_rg_label==0
%            eval(['save ''',direc,b_t,b_t,b_n,b_t1,b_n,common,'.jpg_Files\analysis_n_f'' analysis -V6'])
%        end
%    end
%    %eval(['save ''C:\Users\shhong\Desktop\CJake\Analysis s2 3\JW3_WT\JW3_WT_5FR_v2\JW3_WT_5FR_fix_c.jpg_Files\analysis_n'' analysis_n -V6']);
%end


clear temp temp1
% temp1=imerode(image_in,strel('disk',2));
% temp=temp1 + imdilate(r_label,strel('disk',1)) + imdilate(TRAP_label,strel('disk',1)) + imdilate(AP_label*0.6,strel('disk',1));
% temp(:,:,2)=temp1 + imdilate(g_label,strel('disk',1)) + imdilate(TRAP_label,strel('disk',1)) + imdilate(GFP_label, strel('disk',1)) + imdilate(AP_label*0.3,strel('disk',1));
% temp(:,:,3)=temp1 + imdilate(AP_label*0.3,strel('disk',1)) + imdilate(GFP_label, strel('disk',1));
% temp=temp/max(temp(:));
se1=strel('disk',1);
se=strel('disk',2);
temp1=imerode(image_in,se1);
%%%%%%%%%%%%%%%%%%%%%
%
%	June, 7th 2018
%
if ~RG_label	% 1st label : red,  2nd label : green --> Reversed labels
	temp=temp1*0.5 + (imdilate(g_label,se) + imdilate(TRAP_label,se) + imdilate(AP_label,se)*0.6).*(1-temp1);
	temp(:,:,2)=temp1*0.5 + (imdilate(r_label,se) + imdilate(TRAP_label,se) + imdilate(GFP_label,se) + imdilate(AP_label,se)*0.3).*(1-temp1);
	temp(:,:,3)=temp1*0.5 + (imdilate(AP_label,se)*0.3 + imdilate(GFP_label,se)).*(1-temp1);
else			% 1st label : green,  2nd label : red --> Normal labels
	temp=temp1*0.5 + (imdilate(r_label,se) + imdilate(TRAP_label,se) + imdilate(AP_label,se)*0.6).*(1-temp1);
	temp(:,:,2)=temp1*0.5 + (imdilate(g_label,se) + imdilate(TRAP_label,se) + imdilate(GFP_label,se) + imdilate(AP_label,se)*0.3).*(1-temp1);
	temp(:,:,3)=temp1*0.5 + (imdilate(AP_label,se)*0.3 + imdilate(GFP_label,se)).*(1-temp1);
end
%
%	June, 7th 2018
%
%%%%%%%%%%%%%%%%%%%%%

temp=temp/max(temp(:));


% % if red_first_green_last==1
% %     if b_t==1
% %         temp(:,:,1)=immultiply(uint8(g_label),256);
% %         %temp(:,:,2)=immultiply(uint8(GFP_label),255);
% %         %temp(:,:,3)=imadd(temp(:,:,3),immultiply(uint8(GFP_label),255));
% %         temp(:,:,2)=imadd(temp(:,:,2),immultiply(uint8(r_label),255));
% %         %temp(:,:,1)=imadd(temp(:,:,1),immultiply(uint8(r_label),255));
% %     elseif b_t==2
% %         temp(:,:,1)=immultiply(uint8(r_label),256);
% %         %temp(:,:,2)=immultiply(uint8(GFP_label),255);
% %         %temp(:,:,3)=imadd(temp(:,:,3),immultiply(uint8(GFP_label),255));
% %         temp(:,:,2)=imadd(temp(:,:,2),immultiply(uint8(g_label),255));
% %         %temp(:,:,1)=imadd(temp(:,:,1),immultiply(uint8(r_label),255));
% %     end
% % else
%     temp(:,:,1)=immultiply(uint8(r_label),256);
%     %temp(:,:,2)=immultiply(uint8(GFP_label),255);
%     %temp(:,:,3)=imadd(temp(:,:,3),immultiply(uint8(GFP_label),255));
%     temp(:,:,2)=imadd(temp(:,:,2),immultiply(uint8(g_label),255));
%     %temp(:,:,1)=imadd(temp(:,:,1),immultiply(uint8(r_label),255));    
% % end
% figure;imshow(temp(1:end-down_margin,:,:))

%if write==1
%    if tr_loop==1
%        eval(['imwrite(imrotate(temp(1:end-down_margin,:,:),-rot_angle), ''',direc,c1,'_trim12_tr1.jpg'',''jpg'',''quality'',100);']);
%    elseif tr_loop==2
%        eval(['imwrite(imrotate(temp(1:end-down_margin,:,:),-rot_angle), ''',direc,c1,'_trim12_tr2.jpg'',''jpg'',''quality'',100);']);
%    end
%end
if write==1
%    temp_image=imrotate(temp(1:end-down_margin,:,:),-rot_angle);
    temp_image=temp(1:end-down_margin,:,:);
    if tr_loop==1
        eval(['imwrite(temp_image, ''',direc,c1,'_trim12_tr1.jpg'',''jpg'',''quality'',100);']);
        eval(['imwrite(temp_image(',num2str(min(It)-50),':min(size(temp_image,1),',num2str(max(It)+50),'),',num2str(min(Jt)-50),':',num2str(max(Jt)+50),',:),''',direc,c1,'_trim12_tr1_inside.jpg'',''jpg'',''quality'',100);']);
        %eval(['imwrite(temp_image, ''',archive_d,delimeter,c1,'_trim12_tr1.jpg'',''jpg'',''quality'',100);']);
%         eval(['imwrite(temp_image(',num2str(min(It)-50),':min(size(temp_image,1),',num2str(max(It)+50),'),',num2str(min(Jt)-50),':',num2str(max(Jt)+50),',:),''',archive_d,delimeter,c,nnn,exp_type,'-tm_',num2str(s),'.jpg'',''jpg'',''quality'',100);']);
%         eval(['imwrite(temp_image(',num2str(min(It)-50),':min(size(temp_image,1),',num2str(max(It)+50),'),',num2str(min(Jt)-50),':',num2str(max(Jt)+50),',:),''',archive_d,delimeter,c,nnn,'_',num2str(s),'-tm.jpg'',''jpg'',''quality'',100);']);
        eval(['imwrite(temp_image(',num2str(min(It)-50),':min(size(temp_image,1),',num2str(max(It)+50),'),',num2str(min(Jt)-50),':',num2str(max(Jt)+50),',:),''',archive_d,c1,'_tm.jpg'',''jpg'',''quality'',100);']);
    elseif tr_loop==2
        eval(['imwrite(temp_image, ''',direc,c1,'_trim12_tr2.jpg'',''jpg'',''quality'',100);']);
        eval(['imwrite(temp_image(',num2str(min(It)-50),':min(size(temp_image,1),',num2str(max(It)+50),'),',num2str(min(Jt)-50),':',num2str(max(Jt)+50),',:),''',direc,c1,'_trim12_tr2_inside.jpg'',''jpg'',''quality'',100);']);
        %eval(['imwrite(temp_image, ''',archive_d,delimeter,c1,'_trim12_tr2.jpg'',''jpg'',''quality'',100);']);
        %eval(['imwrite(temp_image(',num2str(min(It)-50),':min(size(temp_image,1),',num2str(max(It)+50),'),',num2str(min(Jt)-50),':',num2str(max(Jt)+50),',:),''',archive_d,delimeter,c1,'_trim12_tr2_inside.jpg'',''jpg'',''quality'',100);']);
    end
end

%if write==1
%    if exp_type=='o'
%        eval(['imwrite(imrotate(temp(1:end-down_margin,:,:),-rot_angle), ''',direc,b_t,b_t,b_n,exp_t,b_t1,b_n,common,'.jpg_Files',b_t1,b_n,common,'_trim12.jpg'',''jpg'',''quality'',100);'])
%    else
%        eval(['imwrite(imrotate(temp(1:end-down_margin,:,:),-rot_angle), ''',direc,b_t,b_t,b_n,b_t1,b_n,common,'.jpg_Files',b_t1,b_n,common,'_trim12.jpg'',''jpg'',''quality'',100);'])
%    end
%end

% temp(:,:,3)=immultiply(uint8(ttt),128);
% temp(:,:,1)=immultiply(uint8(r_label),128);
% %temp(:,:,2)=immultiply(uint8(GFP_label),128);
% %temp(:,:,3)=imadd(temp(:,:,3),immultiply(uint8(GFP_label),128));
% temp(:,:,2)=imadd(temp(:,:,2),immultiply(uint8(g_label),128));
% %temp(:,:,1)=imadd(temp(:,:,1),immultiply(uint8(r_label),128));
% %figure;imshow(temp(1:end-down_margin,:,:))


toc
%v_version
%analysis.surface
analysis_n.surface
analysis_n.volume


close all
end     % for no_tr_region

clear G* T* X Y g* I J col* et height image* mm ou* rG* re* r_* rg* st stats* s_* t t_* te* tm* tt x y no_green no_red

    
        end %s
    end % for number_bone
    
end % for type_bone

%clear % TB test* temp* stats* out* r_* g* G GFP_l* GFP_t* rg*
%matlab_ver=7.1;
%warning('off')
%tic

%ratio=2;
%std_const=2;
%days2_1=7;                         % days between injection 2 and injection 1 (for labels)
%dis_5x=312;                    % 312 pixels in 5x --> 400 um
%dis_10x=dis_5x*2/ratio;              % 312*2 pixels in 10x --> 400 um (based on the image 'C:\Users\shhong\Desktop\CJake\scale.jpg')
%margin_5x=200;                  % 230 pixel

%start_trabecula=250;            % 250 um
%sp_tb=dis_10x; %*start_trabecula/400;      % start point of trabecula from the growth plate (400 um)
%lr_sp_tb=dis_10x*start_trabecula/400;      % start point of trabecula from the endosteum (250 um)
%margi=margin_5x;
%margi_lr=156;
%start_point=550/ratio;
%end_point=4800/ratio;

%red_flag=1;
%GFP_flag=0;
%exact_boundary=1;
%loop_count=20;
%rot_angle=-90;
%thresh=0.0001;

%b_no='3';
%exp_type='o';       % 'n' : normal exposure, 'o' : over_exposure
%with_rg_label=1;    % 1 : using red label and green label to segment DIC
                    % 0 : not using red label and green label to segment DIC


%if str2num(b_no)<7
%    b_type='WT';       % 'WT' : wild type, 'KO' : Knock out
%    b_t='\JW3_WT';
%elseif str2num(b_no)>6
%    b_type='KO';       % 'WT' : wild type, 'KO' : Knock out
%    b_t='\JW3_KO';
%    b_t1='\JW3_KO';     % only for 'analysis1\JW3_KO'
%end

%b_n=['_',b_no,'FR'];
%common='_c';
%direc='C:\Users\shhong\Desktop\CJake\Analysis1';

%end

%select_move_images(direct, d, phr1, samples_per_bone, bone_number1)


%write_data1(direct, phr1, phr2, exp_type, bone_number1, bone_number2)
%write_data3_cell(direct, phr1, phr2, exp_type, bone_number1, bone_number2)
return

















function read_write_images(varargin)
%function read_write_images(direc, dis, color_order, bubble_threshold, phr, samples_per_bone, bone_n1, bone_n2, bone_n3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       Read individual images and Save combined image
%
%direc='C:\Users\shhong\Desktop\C-A labe\Analysis\';


if length(varargin)<7
    error('Not enough argments');
end


directories=varargin{1};
direct=cell2mat(directories(1));
dis=varargin{2};
color_order=varargin{3};
bubble_threshold=varargin{4};
for k = 7:length(varargin)
    eval(['bone_n',num2str(k-6),' = varargin{',num2str(k),'};']);
end

max_dis=max(dis);
no_bone_type=length(varargin)-6;

eval(['a=dir(''', direct,''');'])
if isempty(regexpi([a.name],'write_done.mat'))==0
    return
end

no_flag=0;

for i=1:no_bone_type  
    eval(['bone_number=bone_n',num2str(i),';']);
    %switch i
    %    case 1      % i=1 : Red label is inside(first injected), Green label is outside(last injected)
    %        bone_number=bone_n1;
    %    case 2      % i=2 : Green label is inside(first injected), Red label is outside(last injected) --> normal setting
    %        bone_number=bone_n2;
    %end
    for j=1:length(bone_number)
        %if i==1
        %    break
        %end
        s=1;
        while s<varargin{6}+1 %for s=1:varargin{6}
            for dummy=1
%             if varargin{6}==1
%                 phr=[varargin{5}(i,:), num2str(bone_number(j)), 'V.jpg_Files\'];
%             else
                phr=[varargin{5}(i,:), num2str(bone_number(j)), 'V_',num2str(s),'.jpg_Files\'];
%            end
            direc=[direct, phr];
            nn=num2str(bone_number(j));
            
            c=varargin{5}(i,:);
%                 if varargin{6}==1
%                     c1=[c,nn,'V_c'];
%                 else
            c1=[c,nn,'V_',num2str(s),'_c'];
%                end
                %c1=[c,nn,'_8wk_',num2str(j),'FL_',num2str(s),'_c'];
            eval(['test=dir(''',direc,c1,'0.jpg'');']);
            if isempty(test)==1
                break
            end


            eval(['g=imread(''',direc,c1,'0.jpg'',''jpg'');']);
            %eval(['y=imread(''',direc,c1,'3.jpg'',''jpg'');']);
            eval(['r=imread(''',direc,c1,'1.jpg'',''jpg'');']);
            eval(['a=imread(''',direc,c1,'2.jpg'',''jpg'');']);
            
            no_flag=0;
            a=a(max_dis(2)+1:end,max_dis(1)+1:end);
            %a=a(dis_y+1:end,dis_x+1:end);
            g=g(max_dis(2)-dis(1,2)+1:end-(dis(1,2)),max_dis(1)-dis(1,1)+1:end-(dis(1,1)));
            %y=y(max_dis(2)-dis(2,2)+1:end-(dis(2,2)),max_dis(1)-dis(2,1)+1:end-(dis(2,1)));
            r=r(max_dis(2)-dis(2,2)+1:end-(dis(2,2)),max_dis(1)-dis(2,1)+1:end-(dis(2,1)));
            
            if i==1 | i==2 | i==3
                c=varargin{5}(i,:);
                %c1=[c,'\',c,'_',nn,'FL\',c,'_',nn,'FL.jpg_Files\',c,'_',nn,'FL_c'];
                %c1=[c,nn,'\',c,'Femur_',nn,'_c'];
%                 if varargin{6}==1
%                     c1=[c,nn,'FL_c'];
%                 else
                    c1=[c,nn,'V_',num2str(s),'_c'];
%                end
               
                eval(['imwrite(g,''',direc,c1,'0_shift.jpg'',''jpg'',''quality'',100);']);
                %eval(['imwrite(y,''',direc,c1,'3_shift.jpg'',''jpg'',''quality'',100);']);
                eval(['imwrite(r,''',direc,c1,'1_shift.jpg'',''jpg'',''quality'',100);']);
                eval(['imwrite(a,''',direc,c1,'2_shift.jpg'',''jpg'',''quality'',100);']);
                
                no_flag=0;
            elseif i==4
                c=varargin{5}(i,:);
                %c1=[c,'\',c,'_',nn,'FL\',c,'_',nn,'FL.jpg_Files\',c,'_',nn,'FL_c'];
                c1=[c,nn,'\',c,'Femur_',nn,'_c'];
                eval(['imwrite(a,''',direc,c1,'0_shift.jpg'',''jpg'',''quality'',100);']);
                eval(['imwrite(g,''',direc,c1,'1_shift.jpg'',''jpg'',''quality'',100);']);
                %eval(['y=imread(''',direc,c1,'2.jpg'',''jpg'');']);
                eval(['imwrite(r,''',direc,c1,'2_shift.jpg'',''jpg'',''quality'',100);']);
                
                no_flag=0;
            end

            %aa=imadd(imadd(immultiply(a,0.5),r),y);
            %aa(:,:,2)=imadd(imadd(immultiply(a,0.5),g),y);
            %aa(:,:,3)=immultiply(a,0.5);
            %aa=imadd(imadd(immultiply(a,0.5),r),y);
            %aa(:,:,2)=imadd(immultiply(a,0.5),g);
            %aa(:,:,3)=immultiply(a,0.5);
            aa=imadd(immultiply(a,0.5),r);
            aa(:,:,2)=imadd(immultiply(a,0.5),g);
            aa(:,:,3)=immultiply(a,0.5);
            
            
            %if i==1
            eval(['imwrite(aa,''',direc,c1,'_shift.jpg'',''jpg'',''quality'',100);']);
            %elseif i==2
            %    eval(['imwrite(aa,''',direc,c1,'_shift.jpg'',''jpg'',''quality'',100);']);
            %elseif i==3
            %    eval(['imwrite(aa,''',direc,c1,'_shift.jpg'',''jpg'',''quality'',100);']);
            %end
            clear aa
            
            end % dummy loop
            s=s+1;
        end
    end
end
image_write=1;
eval(['save ''',direct, 'write_done.mat'' ','image_write -v6']); 
return

function [perimeter, area] = perimeter_area (image_in,ratio)
%function [perimeter, area, pp] = perimeter_area (image_interest)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Calculate the perimeter and the area of the input image
%
%   img_interest        : input bw image 
%
%   perimeter           : calculated perimeter
%   area                : calculated area
%   pp                  : boundary image
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    se=[0 1 0;1 0 1;0 1 0];           % structure : diamond with empty center
    %%e=imerode(img,se);                  % erosion
    %%d=img-e                            % boundary of img
    %%if size(img,1)<4 && size(img,2)<4
    %%    d=img
    %%end
    %%figure;imagesc(d)
    %img=imfill(image_interest, 'holes');
    image_interest=double(image_in);
    temp=xcorr2(image_interest,se);
    e=temp(2:end-1,2:end-1);
    d=(e<4).*image_interest;
    %a=diff(img);
    %b=(a==-1);
    %d=zeros(size(img));
    %d(2:end,:)=(a==1);
    %d(1:end-1,:)=d(1:end-1,:)+b;
    %a=diff(img,1,2);
    %b=(a==-1);
    %d(:,2:end)=d(:,2:end)+(a==1);
    %d(:,1:end-1)=d(:,1:end-1)+b;
    %d=double(d>0);
    
    t=[10 2 10;2 1 2;10 2 10];          % matrix to generate patterns for the perimeter calculation
    para=conv2(d,t);                    % 2D convolution
    p=para(2:end-1,2:end-1);            % find the original area of the image
    pp=p.*d;                            % consider only the boundary
    %figure;imagesc(pp)
    
    t1=length(find(pp==5|pp==7|pp==15|pp==17|pp==25|pp==27));   % perimeter contributed by pattern 5, 7, 15, 17, 25, 27
                                                                % it's perimeter is '1' for each pixel
    t2=length(find(pp==21|pp==33|pp==31));                      % perimeter contributed by pattern 21, 33
                                                                % it's perimeter is 'sqrt(2)' for each pixel
    t3=length(find(pp==13|pp==23));                             % perimeter contributed by pattern 13, 23
                                                                % it's perimeter is '(1+sqrt(2))/2' for each pixel
    t4=length(find(pp==3));                                     % perimeter contributed by pattern 3 (only 2 pixels)
                                                                % it's perimeter is '1/2' for each pixel
    t5=length(find(pp==11));                                    % perimeter contributed by pattern 11 (only 2 pixels)
                                                                % it's perimeter is 'sqrt(2)/2' for each pixel
    size(find(pp~=0));                                          % total number of boundary pixels
    
    perimeter=(t1+t2*sqrt(2)+t3*(1+sqrt(2))/2+t4/2+t5*sqrt(2)/2)*ratio;          % perimeter of the boundary of the rotated square
    %area=size(find(image_interest==1),1);                          % area of the rotated square
    area=(length(find(image_interest(:)==1)))*ratio^2;

    %perimeter=area;             % Jan 28, 2010 
%     perimeter=perimeter*rat_mouse;
%     area=area*(rat_mouse^2);
return


function [perimeter, area] = perimeter_area1 (dic, image_in, ratio, exact_boundary)
%function [perimeter, area, pp] = perimeter_area (image_interest)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Modified on August 30, 2013 (surface length (that does not touch DIC) of each pixel 

%   Calculate the perimeter and the area of the input image
%
%   img_interest        : input bw image 
%
%   perimeter           : calculated perimeter
%   area                : calculated area
%   pp                  : boundary image
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    se=[0 1 0;1 0 1;0 1 0];           % structure : diamond with empty center
    dic=double(dic>0);
    d1=xcorr2(dic,se);
    e=d1(2:end-1,2:end-1);

    image_in=double(image_in>0);
    if exact_boundary==1
        e1=dic.*e.*image_in;        % signals are on the edge of DIC
    else
        e1=e.*image_in;             % signals are outside of DIC (1 pixel away)
    end
    t1=length(find(e1==1));         % total # of surfaces that touch DIC on 1 side of pixel
    t3=length(find(e1==3));         % total # of surfaces that touch DIC on 2 sides of pixel
    t2=length(find(e1==2));         % total # of surfaces that touch DIC on 3 sides of pixel
    perimeter=(t1*3+t2*2+t3)*ratio;
    area=length(find(image_in))*ratio^2;
    
return

function [perimeter, area] = perimeter_area2 (dic, image_in, ratio)
%function [perimeter, area, pp] = perimeter_area (image_interest)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Modified on August 30, 2013 (surface length (that does not touch DIC) of each pixel 

%   Calculate the perimeter and the area of the input image
%
%   img_interest        : input bw image 
%
%   perimeter           : calculated perimeter
%   area                : calculated area
%   pp                  : boundary image
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    se=[0 1 0;1 0 1;0 1 0];           % structure : diamond with empty center
    dic=double(dic>0);
    d1=xcorr2(dic,se);
    e=d1(2:end-1,2:end-1);

    image_in=double(image_in>0);
    e1=e.*image_in;
    t1=length(find(e1==1));         % total # of surfaces that touch DIC on 1 side of pixel
    t3=length(find(e1==3));         % total # of surfaces that touch DIC on 2 sides of pixel
    t2=length(find(e1==2));         % total # of surfaces that touch DIC on 3 sides of pixel
    perimeter=(t1*3+t2*2+t3)*ratio;
    area=length(find(image_in))*ratio^2;
    
return

function [r_trim, r_trim1, thickness_r, g_trim, g_trim1, thickness_g, TB]=move_control(i, s, sd, trim_template, im, r_test, g_test, GFP_label, GFP_label1, stats_DIC, J, TB)
sx=s(1); ex=s(2); sy=s(3); ey=s(4);
sdx=sd(1); edx=sd(2); sdy=sd(3); edy=sd(4);
r_trim=trim_template;
r_trim1=trim_template;
g_trim=trim_template;
g_trim1=trim_template;
thickness_r=0;
thickness_g=0;
% TB(i)=struct;

[L_g n_g]=bwlabel(g_test);
    stats_g=regionprops(L_g,'Image','BoundingBox','Orientation','Area','Centroid'); clear L_g
    ori_g=[stats_g.Orientation];
    [L_r n_r]=bwlabel(r_test);
    stats_r=regionprops(L_r,'Image','BoundingBox','Orientation','Area','Centroid'); clear L_r
    ori_r=[stats_r.Orientation];
    
    if stats_DIC(J(i)).Orientation<0
        ref_ang=180+stats_DIC(J(i)).Orientation;
    else
        ref_ang=stats_DIC(J(i)).Orientation;
    end
    tt_r=zeros(1,length(ori_r));
    tt_r(ori_r<0)=180+ori_r(ori_r<0);
    tt_r(ori_r>=0)=ori_r(ori_r>=0);
    test_r_ori=tt_r-ref_ang;
    r_orient_order=sort(abs(test_r_ori));
    
    tt_g=zeros(1,length(ori_g));
    tt_g(ori_g<0)=180+ori_g(ori_g<0);
    tt_g(ori_g>=0)=ori_g(ori_g>=0);
    test_g_ori=tt_g-ref_ang;
    g_orient_order=sort(abs(test_g_ori));
    
    red_im=r_test;
    
    if isempty(find(g_test)==1)==0 && isempty(find(r_test)==1)==0        % green and red in DIC
        r_trim_temp=zeros(size(r_test));
        for j=1:n_r
            %start_r_x=round(stats_r(j).BoundingBox(2));
            %end_r_x=start_r_x+stats_r(j).BoundingBox(4)-1;
            %start_r_y=round(stats_r(j).BoundingBox(1));
            %end_r_y=start_r_y+stats_r(j).BoundingBox(3)-1;
            %if j==7
            ['n_r = ',num2str(j)];
            %    j;
            %end
            if stats_r(j).Area~=1
                
%       February 28, 2013
%                 % June 9, 2008 [r_trim_temp]=move5(im,stats_r(j));
%                 [r_trim_temp, r_trim_alt_temp]=move_red_7(im,stats_r(j));
%                 [Id Jd Vd]=find(r_trim_temp);
%                 thickness_r(j)=mean(Vd);
%                 [Id Jd Vd]=find(r_trim_alt_temp);
%                 thickness_r(j)=thickness_r(j)+mean(Vd);
                ind=find(abs(test_r_ori)==r_orient_order(j));
                for ind_loop=1:length(ind)
                    [r_trim_temp, r_trim_alt_temp]=move_red_7(im,stats_r(ind(ind_loop)));
%                     [r_trim_temp, r_trim_alt_temp]=move_7(im,stats_r(ind(ind_loop)));           % 2/3/2016
                    [Id Jd Vd]=find(r_trim_temp);
                    thickness_r(ind(ind_loop))=mean(Vd);
                    [Id Jd Vd]=find(r_trim_alt_temp);
                    thickness_r(ind(ind_loop))=thickness_r(ind(ind_loop))+mean(Vd);
                end

%       February 28, 2013


                tmp=double(r_trim_temp>0);
% 12/22/09              r_trim(start_x:end_x,start_y:end_y)=r_trim(start_x:end_x,start_y:end_y)+tmp;
                r_trim(sx:ex,sy:ey)=r_trim(sx:ex,sy:ey)+tmp(sdx+1:sdx+stats_DIC(J(i)).BoundingBox(4),sdy+1:sdy+stats_DIC(J(i)).BoundingBox(3));
% 1/7/08                r_trim1(start_x:end_x,start_y:end_y)=r_trim1(start_x:end_x,start_y:end_y)+r_trim_temp;
% 12/22/09              [r_trim1(start_x:end_x,start_y:end_y)]=find_min(r_trim1(start_x:end_x,start_y:end_y), r_trim_temp);
                [r_trim1(sx:ex,sy:ey)]=find_min(r_trim1(sx:ex,sy:ey), r_trim_temp(sdx+1:sdx+stats_DIC(J(i)).BoundingBox(4),sdy+1:sdy+stats_DIC(J(i)).BoundingBox(3)));
            end
        end
        
        g_trim_wo_red_temp=zeros(size(g_test));
        for j=1:n_g
            %start_r_x=round(stats_g(j).BoundingBox(2));
            %end_r_x=start_r_x+stats_g(j).BoundingBox(4)-1;
            %start_r_y=round(stats_g(j).BoundingBox(1));
            %end_r_y=start_r_y+stats_g(j).BoundingBox(3)-1;
            ['n_g = ',num2str(j)];
            ang=abs(ori_r-ori_g(j));
            ind=find(ang<30);
            
%            red_im=zeros(size(r_test));
%            if isempty(ind)==0
%                for p=1:length(ind)
%                    st_x=round(stats_r(ind(p)).BoundingBox(2));
%                    st_y=round(stats_r(ind(p)).BoundingBox(1));
%                    e_x=st_x+stats_r(ind(p)).BoundingBox(4)-1;
%                    e_y=st_y+stats_r(ind(p)).BoundingBox(3)-1;
%                    red_im(st_x:e_x,st_y:e_y)=red_im(st_x:e_x,st_y:e_y)+stats_r(ind(p)).Image;
%                end
%            end
            
            %if length(ind)>1
            %    ind=find(ang==min(ang));
            %end
            
            
            if isempty(ind)==0 && stats_g(j).Area~=1
                %[direction]=find_direction(stats_g,j,stats_r,ind);
%                [g_trim_temp]=move_4_temp(stats_g(j),stats_r(ind),stats_DIC(i).Image);
%                [g_trim_temp]=move_4(im,stats_g(j),stats_r(ind));
%                [g_trim_temp, GFP_trim_temp]=move_4(im,stats_g(j),r_trim1(start_x:end_x,start_y:end_y));
%                [g_trim_temp, g_trim_alt_temp, GFP_trim_temp]=move_5(im,stats_g(j),r_trim1(start_x:end_x,start_y:end_y));       % 1/23/2008
                [g_trim_temp, g_trim_alt_temp, GFP_trim_temp]=move_7(im,stats_g(j),r_test);       % 6/5/2008
%                 [g_trim_temp, g_trim_alt_temp, GFP_trim_temp]=move_7(im,stats_g(j));       % 2/3/2016
                [Id Jd Vd]=find(g_trim_temp);
                thickness_g(j)=mean(Vd);
                [Id Jd Vd]=find(g_trim_alt_temp);
                thickness_g(j)=thickness_g(j)+mean(Vd);

                [col_trim, row_trim, val_trim]=find(g_trim_temp);
                TB(i).G_label(j).trim.col=col_trim;
                TB(i).G_label(j).trim.row=row_trim;
                TB(i).G_label(j).trim.val=val_trim;
                TB(i).G_label(j).trim.size_trim_ver=size(g_trim_temp,1);
                TB(i).G_label(j).trim.size_trim_hor=size(g_trim_temp,2);
%                 TB(i).G_label(j).trim.start_x=start_x;
%                 TB(i).G_label(j).trim.end_x=end_x;
%                 TB(i).G_label(j).trim.start_y=start_y;
%                 TB(i).G_label(j).trim.end_y=end_y;
                TB(i).G_label(j).trim.start_x=1;
                TB(i).G_label(j).trim.end_x=size(g_trim_temp,1);
                TB(i).G_label(j).trim.start_y=1;
                TB(i).G_label(j).trim.end_y=size(g_trim_temp,2);
                
                
                [col_trim, row_trim, val_trim]=find(g_trim_alt_temp);
                TB(i).G_label(j).trim_alt.col=col_trim;
                TB(i).G_label(j).trim_alt.row=row_trim;
                TB(i).G_label(j).trim_alt.val=val_trim;
                TB(i).G_label(j).trim_alt.size_trim_ver=size(g_trim_alt_temp,1);
                TB(i).G_label(j).trim_alt.size_trim_hor=size(g_trim_alt_temp,2);
%                 TB(i).G_label(j).trim_alt.start_x=start_x;
%                 TB(i).G_label(j).trim_alt.end_x=end_x;
%                 TB(i).G_label(j).trim_alt.start_y=start_y;
%                 TB(i).G_label(j).trim_alt.end_y=end_y;
                TB(i).G_label(j).trim_alt.start_x=1;
                TB(i).G_label(j).trim_alt.end_x=size(g_trim_alt_temp,1);
                TB(i).G_label(j).trim_alt.start_y=1;
                TB(i).G_label(j).trim_alt.end_y=size(g_trim_alt_temp,2);

                tmp=double(g_trim_temp>0);
% 12/22/09                g_trim(start_x:end_x,start_y:end_y)=g_trim(start_x:end_x,start_y:end_y)+tmp;
                g_trim(sx:ex,sy:ey)=g_trim(sx:ex,sy:ey)+tmp(sdx+1:sdx+stats_DIC(J(i)).BoundingBox(4),sdy+1:sdy+stats_DIC(J(i)).BoundingBox(3));
% 1/7/08                g_trim1(start_x:end_x,start_y:end_y)=g_trim1(start_x:end_x,start_y:end_y)+g_trim_temp;
% 12/22/09                [g_trim1(start_x:end_x,start_y:end_y)]=find_min(g_trim1(start_x:end_x,start_y:end_y), g_trim_temp);
                [g_trim1(sx:ex,sy:ey)]=find_min(g_trim1(sx:ex,sy:ey), g_trim_temp(sdx+1:sdx+stats_DIC(J(i)).BoundingBox(4),sdy+1:sdy+stats_DIC(J(i)).BoundingBox(3)));
                tmp=double(GFP_trim_temp>0);
% 12/22/09                GFP_label(start_x:end_x,start_y:end_y)=GFP_label(start_x:end_x,start_y:end_y)+tmp;
                GFP_label(sx:ex,sy:ey)=GFP_label(sx:ex,sy:ey)+tmp(sdx+1:sdx+stats_DIC(J(i)).BoundingBox(4),sdy+1:sdy+stats_DIC(J(i)).BoundingBox(3));
% 1/7/08                GFP_label1(start_x:end_x,start_y:end_y)=GFP_label1(start_x:end_x,start_y:end_y)+GFP_trim_temp;
% 12/22/09                [GFP_label1(start_x:end_x,start_y:end_y)]=find_min(GFP_label1(start_x:end_x,start_y:end_y), GFP_trim_temp);
                [GFP_label1(sx:ex,sy:ey)]=find_min(GFP_label1(sx:ex,sy:ey), GFP_trim_temp(sdx+1:sdx+stats_DIC(J(i)).BoundingBox(4),sdy+1:sdy+stats_DIC(J(i)).BoundingBox(3)));
            elseif isempty(ind)==1 && stats_g(j).Area~=1
%                [g_trim_wo_red_temp]=move_4(im,stats_g(j));
%                [g_trim_wo_red_temp,GFP_trim_temp]=move_4(im,stats_g(j),r_trim1(start_x:end_x,start_y:end_y));
%                [g_trim_wo_red_temp, g_trim_wo_red_alt_temp, GFP_trim_temp]=move_5(im,stats_g(j),r_trim1(start_x:end_x,start_y:end_y));       % 1/23/2008
                [g_trim_wo_red_temp, g_trim_wo_red_alt_temp, GFP_trim_temp]=move_7(im,stats_g(j),r_test);       % 6/5/2008
%                 [g_trim_wo_red_temp, g_trim_wo_red_alt_temp, GFP_trim_temp]=move_7(im,stats_g(j));       % 2/3/2016
                [Id Jd Vd]=find(g_trim_wo_red_temp);
                thickness_g(j)=mean(Vd);
                [Id Jd Vd]=find(g_trim_wo_red_alt_temp);
                thickness_g(j)=thickness_g(j)+mean(Vd);

                [col_trim, row_trim, val_trim]=find(g_trim_wo_red_temp);
                TB(i).G_label(j).trim.col=col_trim;
                TB(i).G_label(j).trim.row=row_trim;
                TB(i).G_label(j).trim.val=val_trim;
                TB(i).G_label(j).trim.size_trim_ver=size(g_trim_wo_red_temp,1);
                TB(i).G_label(j).trim.size_trim_hor=size(g_trim_wo_red_temp,2);
%                 TB(i).G_label(j).trim.start_x=start_x;
%                 TB(i).G_label(j).trim.end_x=end_x;
%                 TB(i).G_label(j).trim.start_y=start_y;
%                 TB(i).G_label(j).trim.end_y=end_y;
                TB(i).G_label(j).trim.start_x=1;
                TB(i).G_label(j).trim.end_x=size(g_trim_wo_red_temp,1);
                TB(i).G_label(j).trim.start_y=1;
                TB(i).G_label(j).trim.end_y=size(g_trim_wo_red_temp,2);
                
                
                [col_trim, row_trim, val_trim]=find(g_trim_wo_red_alt_temp);
                TB(i).G_label(j).trim_alt.col=col_trim;
                TB(i).G_label(j).trim_alt.row=row_trim;
                TB(i).G_label(j).trim_alt.val=val_trim;
                TB(i).G_label(j).trim_alt.size_trim_ver=size(g_trim_wo_red_alt_temp,1);
                TB(i).G_label(j).trim_alt.size_trim_hor=size(g_trim_wo_red_alt_temp,2);
%                 TB(i).G_label(j).trim_alt.start_x=start_x;
%                 TB(i).G_label(j).trim_alt.end_x=end_x;
%                 TB(i).G_label(j).trim_alt.start_y=start_y;
%                 TB(i).G_label(j).trim_alt.end_y=end_y;
                TB(i).G_label(j).trim_alt.start_x=1;
                TB(i).G_label(j).trim_alt.end_x=size(g_trim_wo_red_alt_temp,1);
                TB(i).G_label(j).trim_alt.start_y=1;
                TB(i).G_label(j).trim_alt.end_y=size(g_trim_wo_red_alt_temp,2);


                tmp=double(g_trim_wo_red_temp>0);
                g_trim(sx:ex,sy:ey)=g_trim(sx:ex,sy:ey)+tmp(sdx+1:sdx+stats_DIC(J(i)).BoundingBox(4),sdy+1:sdy+stats_DIC(J(i)).BoundingBox(3));
% 12/22/09                g_trim(start_x:end_x,start_y:end_y)=g_trim(start_x:end_x,start_y:end_y)+tmp;
% 1/7/08                g_trim1(start_x:end_x,start_y:end_y)=g_trim1(start_x:end_x,start_y:end_y)+g_trim_wo_red_temp;
                [g_trim1(sx:ex,sy:ey)]=find_min(g_trim1(sx:ex,sy:ey), g_trim_wo_red_temp(sdx+1:sdx+stats_DIC(J(i)).BoundingBox(4),sdy+1:sdy+stats_DIC(J(i)).BoundingBox(3)));
% 12/22/09                [g_trim1(start_x:end_x,start_y:end_y)]=find_min(g_trim1(start_x:end_x,start_y:end_y), g_trim_wo_red_temp);
                tmp=double(GFP_trim_temp>0);
                GFP_label(sx:ex,sy:ey)=GFP_label(sx:ex,sy:ey)+tmp(sdx+1:sdx+stats_DIC(J(i)).BoundingBox(4),sdy+1:sdy+stats_DIC(J(i)).BoundingBox(3));
% 12/22/09                GFP_label(start_x:end_x,start_y:end_y)=GFP_label(start_x:end_x,start_y:end_y)+tmp;
% 1/7/08                GFP_label1(start_x:end_x,start_y:end_y)=GFP_label1(start_x:end_x,start_y:end_y)+GFP_trim_temp;
                [GFP_label1(sx:ex,sy:ey)]=find_min(GFP_label1(sx:ex,sy:ey), GFP_trim_temp(sdx+1:sdx+stats_DIC(J(i)).BoundingBox(4),sdy+1:sdy+stats_DIC(J(i)).BoundingBox(3)));
% 12/22/09                [GFP_label1(start_x:end_x,start_y:end_y)]=find_min(GFP_label1(start_x:end_x,start_y:end_y), GFP_trim_temp);
%                g_trim_wo_red_temp(start_r_x:end_r_x,start_r_y:end_r_y)=g_trim_wo_red_temp(start_r_x:end_r_x,start_r_y:end_r_y)+stats_g(j).Image;
            end
        end
% 1/3/08        g_trim_wo_red(start_x:end_x,start_y:end_y)=g_trim_wo_red(start_x:end_x,start_y:end_y)+g_trim_wo_red_temp;
        
    elseif isempty(find(g_test)==1)==1 && isempty(find(r_test)==1)==0        % red, but no green in DIC
        r_trim_temp=zeros(size(r_test));
        for j=1:n_r
            %start_r_x=round(stats_r(j).BoundingBox(2));
            %end_r_x=start_r_x+stats_r(j).BoundingBox(4)-1;
            %start_r_y=round(stats_r(j).BoundingBox(1));
            %end_r_y=start_r_y+stats_r(j).BoundingBox(3)-1;
            ['n_r with no green = ',num2str(j)];
            if stats_r(j).Area~=1
                
%       February 28, 2013
% %                [r_trim_temp]=move_4(im,stats_r(j));
%                 [r_trim_temp, r_trim_alt_temp]=move_5(im,stats_r(j));
%                 [Id Jd Vd]=find(r_trim_temp);
%                 thickness_r(j)=mean(Vd);
%                 [Id Jd Vd]=find(r_trim_alt_temp);
%                 thickness_r(j)=thickness_r(j)+mean(Vd);
                ind=find(abs(test_r_ori)==r_orient_order(j));
                for ind_loop=1:length(ind)
%                     [r_trim_temp, r_trim_alt_temp]=move_5(im,stats_r(ind(ind_loop)));
                    [r_trim_temp, r_trim_alt_temp]=move_7(im,stats_r(ind(ind_loop)));           % 2/3/2016
                    [Id Jd Vd]=find(r_trim_temp);
                    thickness_r(ind(ind_loop))=mean(Vd);
                    [Id Jd Vd]=find(r_trim_alt_temp);
                    thickness_r(ind(ind_loop))=thickness_r(ind(ind_loop))+mean(Vd);
                end
%       February 28, 2013

                tmp=double(r_trim_temp>0);
                r_trim(sx:ex,sy:ey)=r_trim(sx:ex,sy:ey)+tmp(sdx+1:sdx+stats_DIC(J(i)).BoundingBox(4),sdy+1:sdy+stats_DIC(J(i)).BoundingBox(3));
% 12/22/09                r_trim(start_x:end_x,start_y:end_y)=r_trim(start_x:end_x,start_y:end_y)+tmp;
% 1/7/08                r_trim1(start_x:end_x,start_y:end_y)=r_trim1(start_x:end_x,start_y:end_y)+r_trim_temp;
                [r_trim1(sx:ex,sy:ey)]=find_min(r_trim1(sx:ex,sy:ey), r_trim_temp(sdx+1:sdx+stats_DIC(J(i)).BoundingBox(4),sdy+1:sdy+stats_DIC(J(i)).BoundingBox(3)));
% 12/22/09                [r_trim1(start_x:end_x,start_y:end_y)]=find_min(r_trim1(start_x:end_x,start_y:end_y), r_trim_temp);
            end
        end
        
    elseif isempty(find(g_test)==1)==0 && isempty(find(r_test)==1)==1        % green, but no red in DIC
        for j=1:n_g
            %start_r_x=round(stats_g(j).BoundingBox(2));
            %end_r_x=start_r_x+stats_g(j).BoundingBox(4)-1;
            %start_r_y=round(stats_g(j).BoundingBox(1));
            %end_r_y=start_r_y+stats_g(j).BoundingBox(3)-1;
            ['n_g with no green = ',num2str(j)];
            if stats_g(j).Area~=1
%       February 28, 2013
%  %               [g_trim_wo_red_temp]=move_4(im,stats_g(j));
%                 [g_trim_wo_red_temp, g_trim_wo_red_alt_temp]=move_5(im,stats_g(j));       % 1/23/2008
%                 [Id Jd Vd]=find(g_trim_wo_red_temp);
%                 thickness_g(j)=mean(Vd);
%                 [Id Jd Vd]=find(g_trim_wo_red_alt_temp);
%                 thickness_g(j)=thickness_g(j)+mean(Vd);
                ind=find(abs(test_g_ori)==g_orient_order(j));
                for ind_loop=1:length(ind)
%                     [g_trim_wo_red_temp, g_trim_wo_red_alt_temp]=move_5(im,stats_g(ind(ind_loop)));       % 1/23/2008
                    [g_trim_wo_red_temp, g_trim_wo_red_alt_temp]=move_7(im,stats_g(ind(ind_loop)));       % 2/3/2016
                    [Id Jd Vd]=find(g_trim_wo_red_temp);
                    thickness_g(ind(ind_loop))=mean(Vd);
                    [Id Jd Vd]=find(g_trim_wo_red_alt_temp);
                    thickness_g(ind(ind_loop))=thickness_g(ind(ind_loop))+mean(Vd);
                end
%       February 28, 2013

                [col_trim, row_trim, val_trim]=find(g_trim_wo_red_temp);
                TB(i).G_label(j).trim.col=col_trim;
                TB(i).G_label(j).trim.row=row_trim;
                TB(i).G_label(j).trim.val=val_trim;
                TB(i).G_label(j).trim.size_trim_ver=size(g_trim_wo_red_temp,1);
                TB(i).G_label(j).trim.size_trim_hor=size(g_trim_wo_red_temp,2);
%                 TB(i).G_label(j).trim.start_x=start_x;
%                 TB(i).G_label(j).trim.end_x=end_x;
%                 TB(i).G_label(j).trim.start_y=start_y;
%                 TB(i).G_label(j).trim.end_y=end_y;
                TB(i).G_label(j).trim.start_x=1;
                TB(i).G_label(j).trim.end_x=size(g_trim_wo_red_temp,1);
                TB(i).G_label(j).trim.start_y=1;
                TB(i).G_label(j).trim.end_y=size(g_trim_wo_red_temp,2);

                
                
                [col_trim, row_trim, val_trim]=find(g_trim_wo_red_alt_temp);
                TB(i).G_label(j).trim_alt.col=col_trim;
                TB(i).G_label(j).trim_alt.row=row_trim;
                TB(i).G_label(j).trim_alt.val=val_trim;
                TB(i).G_label(j).trim_alt.size_trim_ver=size(g_trim_wo_red_alt_temp,1);
                TB(i).G_label(j).trim_alt.size_trim_hor=size(g_trim_wo_red_alt_temp,2);
%                 TB(i).G_label(j).trim_alt.start_x=start_x;
%                 TB(i).G_label(j).trim_alt.end_x=end_x;
%                 TB(i).G_label(j).trim_alt.start_y=start_y;
%                 TB(i).G_label(j).trim_alt.end_y=end_y;
                TB(i).G_label(j).trim_alt.start_x=1;
                TB(i).G_label(j).trim_alt.end_x=size(g_trim_wo_red_alt_temp,1);
                TB(i).G_label(j).trim_alt.start_y=1;
                TB(i).G_label(j).trim_alt.end_y=size(g_trim_wo_red_alt_temp,2);


                tmp=double(g_trim_wo_red_temp>0);
                g_trim(sx:ex,sy:ey)=g_trim(sx:ex,sy:ey)+tmp(sdx+1:sdx+stats_DIC(J(i)).BoundingBox(4),sdy+1:sdy+stats_DIC(J(i)).BoundingBox(3));
% 12/22/09                g_trim(start_x:end_x,start_y:end_y)=g_trim(start_x:end_x,start_y:end_y)+tmp;
% 1/7/08                g_trim1(start_x:end_x,start_y:end_y)=g_trim1(start_x:end_x,start_y:end_y)+g_trim_wo_red_temp;
                [g_trim1(sx:ex,sy:ey)]=find_min(g_trim1(sx:ex,sy:ey), g_trim_wo_red_temp(sdx+1:sdx+stats_DIC(J(i)).BoundingBox(4),sdy+1:sdy+stats_DIC(J(i)).BoundingBox(3)));
% 12/22/09                [g_trim1(start_x:end_x,start_y:end_y)]=find_min(g_trim1(start_x:end_x,start_y:end_y), g_trim_wo_red_temp);
                %g_trim_wo_red(start_x:end_x,start_y:end_y)=g_trim_wo_red(start_x:end_x,start_y:end_y)+g_test;%(start_x:end_x,start_y:end_y);
            end
        end
    end
    return



function [trim, trim_alt]=move_red_7(DIC_image,stats_g,red1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Modified on June 5, 2008
%
%   ignore red1
%   regard stats_g as red
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

area=sum(sum(stats_g.Image));
dx=size(DIC_image,1);
dy=size(DIC_image,2);

if nargin==3
    red=(red1>0);
    boundary=DIC_image-imerode(DIC_image,strel('disk',1));
    [L n]=bwlabel(red);
    stats_red=regionprops(L, 'Orientation', 'Centroid', 'Area'); clear L
    ang=abs([stats_red.Orientation]-stats_g.Orientation);
    ind=find(ang==min(ang));
    
    %d=sum(((cat(1,stats_red.Centroid)-(stats_g.Centroid'*ones(1,n))').^2),2);
    %d1=find(d==min(d(:)));
    ori_r=[stats_g.Orientation];
    %ori_r=stats_red(ind).Orientation;
elseif nargin==2
    red=zeros(size(DIC_image));
    boundary=DIC_image-imerode(DIC_image,strel('disk',1));
    ori_g=[stats_g.Orientation];
    ori_r=[stats_g.Orientation];
end


%ori_temp=abs(ori_g-ori_r);
%if ori_temp>90
%    ori_temp=ori_temp-90;
%end
%if ori_temp>30
%    ori=ori_g;
%else
%    ori=ori_r;
%end

out1=zeros(size(boundary)); out_test1=zeros(size(boundary));
out2=zeros(size(boundary)); out_test2=zeros(size(boundary)); trim=zeros(size(boundary));
GFP_trim=zeros(size(boundary));

%if nargin==3
%    red=out1;
%end

ori=ori_r;

if ori<=0
    direction=(ori+90)*pi/180;
    if ori+90>45
        x_first=1;
    else
        x_first=0;
    end
elseif ori>0
    direction=(ori-90)*pi/180;
    if ori-90<-45
        x_first=1;
    else
        x_first=0;
    end
end
ratio=tan(direction);


%if nargin==3
%    start_x=round(stats_r.BoundingBox(2));
%    end_x=start_x+stats_r.BoundingBox(4)-1;
%    start_y=round(stats_r.BoundingBox(1));
%    end_y=start_y+stats_r.BoundingBox(3)-1;
%    red(start_x:end_x,start_y:end_y)=stats_r.Image;
%end

template=stats_g.Image;
start_x=round(stats_g.BoundingBox(2));
end_x=start_x+stats_g.BoundingBox(4)-1;
start_y=round(stats_g.BoundingBox(1));
end_y=start_y+stats_g.BoundingBox(3)-1;

new_g=zeros(size(boundary));
new_g(start_x:end_x,start_y:end_y)=template;
out=boundary & new_g;

test=1;
count=0;
inside_bone_1=0;
while test==1
    new_g=zeros(size(out1));
    count=count+1;
    if x_first==0
        st_y=start_y+(count);
        ed_y=end_y+(count);
        st_x=start_x-round(ratio*(count));
        ed_x=end_x-round(ratio*(count));
    else
        st_x=start_x+(count);
        ed_x=end_x+(count);
        st_y=start_y-round((count)/ratio);
        ed_y=end_y-round((count)/ratio);
    end
    if st_x>size(boundary,1) | ed_x>size(boundary,1) | ed_y>size(boundary,2) | st_y>size(boundary,2) | st_x<1 | st_y<1 | ed_x<1 | ed_y<1
        if count==1
            dis1_1(count)=0;
            dis1(1:2)=[0 0];
            if nargin==3
                dis1(1:2)=[0 0];
            end
            break
        end
        break
    end
    new_g(st_x:ed_x,st_y:ed_y)=template;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Modified on June 5, 2008
%
    %if nargin==2
        new_g=new_g & DIC_image;
        if sum(sum(new_g))/area>0.70
            inside_bone_1=1;
        end
        template=new_g(st_x:ed_x,st_y:ed_y);
    %elseif nargin==3
    %    new_g=new_g & DIC_image & red;
    %end
    if (x_first==1 & (count)>round(dx/3)) |(x_first==0 & (count)>round(dy/3)) | sum(sum(new_g))/area==0%<0.30
        if count==1
            dis1_1(count)=0;
            if nargin==2
                dis1(1:2)=[0 0];
            end
        end
        break
    end
    %if sum(sum(new_g))/area<0.20 %isempty(find(new_g))==1
    %    if count==1
    %        dis1_1(count)=0;
    %        if nargin==3
    %            dis1=0;
    %        end
    %    end
    %    break
    %end
%    out1=min(out1,count*(new_g & boundary));
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   find the mid-point of the labels instead of leading edge
%
%   December 01, 2010
%
   out1=out1+(1-double(out1>0)).*((count)*double(new_g & boundary));    % 2/3/2016
%     out1=out1+double(new_g & boundary);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if nargin==3
        out_test1=out_test1+(new_g & red);
        dis1(count)=size(find((new_g & red)==1),1);
        dis1_1(count)=size(find((new_g & boundary)==1),1);
    elseif nargin==2
        %out_test1=out_test1+(new_g);
        dis1(count)=size(find((new_g & red)==1),1);
        dis1_1(count)=size(find((new_g & boundary)==1),1);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   find the mid-point of the labels instead of leading edge
%
%   December 01, 2010
%
% out1=round(out1/2);    % 2/3/2016
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


template=stats_g.Image;
count=0;
inside_bone_2=0;
while test==1
    new_g=zeros(size(out1));
    count=count-1;
    if x_first==0
        st_y=start_y+(count);
        ed_y=end_y+(count);
        st_x=start_x-round(ratio*(count));
        ed_x=end_x-round(ratio*(count));
    else
        st_x=start_x+(count);
        ed_x=end_x+(count);
        st_y=start_y-round((count)/ratio);
        ed_y=end_y-round((count)/ratio);
    end
    if st_x>size(boundary,1) | ed_x>size(boundary,1) | ed_y>size(boundary,2) | st_y>size(boundary,2) | st_x<1 | st_y<1 | ed_x<1 | ed_y<1
        if -count==1
%            dis1_1(count)=0;
            dis2_1(-count)=0;
            dis2(1:2)=[0 0];
            if nargin==3
                dis2(1:2)=[0 0];
            end
            break
        end
        break
    end
    new_g(st_x:ed_x,st_y:ed_y)=template;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Modified on June 5, 2008
%
    %if nargin==2
        new_g=new_g & DIC_image;
        if sum(sum(new_g))/area>0.70
            inside_bone_2=1;
        end
        template=new_g(st_x:ed_x,st_y:ed_y);
    %elseif nargin==3
    %    new_g=new_g & DIC_image & red;
    %end
    if (x_first==1 & -(count)>round(dx/3)) |(x_first==0 & -(count)>round(dy/3)) | sum(sum(new_g))/area==0%<0.30
        if count==-1
            dis2_1(-count)=0;
            if nargin==2
                dis2(1:2)=[0 0];
            end
        end
        break
    end
    %if sum(sum(new_g))/area<0.05 %isempty(find(new_g))==1
    %    if count==-1
    %        dis2_1(-count)=0;
    %        if nargin==3
    %            dis2=0;
    %        end
    %    end
    %    break
    %end
    
%    out2=min(out2,-count*(new_g & boundary));
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   find the mid-point of the labels instead of leading edge
%
%   December 01, 2010
%
  out2=out2+(1-double(out2>0)).*(-(count)*double(new_g & boundary));    % 2/3/2016
%     out2=out2+double(new_g & boundary);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if nargin==3
        out_test2=out_test2+(new_g & red);
        dis2(-count)=size(find((new_g & red)==1),1);
        dis2_1(-count)=size(find((new_g & boundary)==1),1);
    elseif nargin==2
        %out_test2=out_test2+(new_g);
        dis2(-count)=size(find((new_g & red)==1),1);
        dis2_1(-count)=size(find((new_g & boundary)==1),1);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   find the mid-point of the labels instead of leading edge
%
%   December 01, 2010
%
% out2=round(out2/2);       % 2/3/2016
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Modified on June 5, 2008
%
if nargin==4                % originally it was 'if nargin==3'
%
%                           % bypass this routine
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    common_1=out1 & red1;
    test_out_1=double(common_1).*out1;
    test_red_1=double(common_1).*red1;
%    if isempty(find(test_out_1==1))==1              % Matlab 6.5
    if isempty(find(test_out_1,1))==1              % Matlab 7.1
        ave_dis_out_1=sum(test_out_1(:));
    else
        ave_dis_out_1=sum(test_out_1(:))/length(find(test_out_1));
    end
%    if isempty(find(test_red_1==1))==1               % Matlab 6.5
    if isempty(find(test_red_1,1))==1               % Matlab 7.1
        ave_dis_red_1=sum(test_red_1(:));
    else
        ave_dis_red_1=sum(test_red_1(:))/length(find(test_red_1));
    end
    common_2=out2 & red1;
    test_out_2=double(common_2).*out2;
    test_red_2=double(common_2).*red1;
%    if isempty(find(test_out_2==1))==1               % Matlab 6.5
    if isempty(find(test_out_2,1))==1               % Matlab 7.1
        ave_dis_out_2=sum(test_out_2(:));
    else
        ave_dis_out_2=sum(test_out_2(:))/length(find(test_out_2));
    end
%    if isempty(find(test_red_2==1))==1               % Matlab 6.5
    if isempty(find(test_red_2,1))==1               % Matlab 7.1
        ave_dis_red_2=sum(test_red_2(:));
    else
        ave_dis_red_2=sum(test_red_2(:))/length(find(test_red_2));
    end
    
    if ave_dis_out_1~=0 & sum(test_out_1(:))>sum(test_out_2(:))
        if ave_dis_out_1>ave_dis_red_1
            trim=out1;
            trim_alt=out2;                          % alternate trim candidate
        else
            trim=zeros(size(out1));
            trim_alt=zeros(size(out1));
            GFP_trim=out1;
        end
    elseif ave_dis_out_2~=0 & sum(test_out_2(:))>sum(test_out_1(:))
        if ave_dis_out_2>ave_dis_red_2
            trim=out2;
            trim_alt=out1;                          % alternate trim candidate
        else
            trim=zeros(size(out1));
            trim_alt=zeros(size(out1));
            GFP_trim=out2;
        end
    else
        trim=zeros(size(out1));
        trim_alt=zeros(size(out1));
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Modified on June 5, 2008
%
elseif nargin==2 | nargin==3    % originally it was 'elseif nargin==2'
%
%                               % all process goes this routine
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if max(dis1_1)==0
        J1=0;
    else
        [I1 J1]=find(dis1_1==max(dis1_1));
    end
    if max(dis2_1)==0
        J2=0;
    else
        [I2 J2]=find(dis2_1==max(dis2_1));
    end
    if max(dis1(2:end))==0
        Jr1=0;
    else
        [Ir1 Jr1]=find(dis1==max(dis1));    %   added on June 5, 2008
    end
    if max(dis2)==0
        Jr2=0;
    else
        [Ir2 Jr2]=find(dis2==max(dis2));    %   added on June 5, 2008
    end
    
    if max(dis1)>0 & max(dis2)>0
        if max(dis2)/max(dis1)>1.5            
            Jr1=0;
        elseif max(dis1)/max(dis2)>1.5        
            Jr2=0;
        end
    end
    
    if min(Jr1)>0 & min(Jr2)>0
        if min(Jr1)>min(Jr2)
            trim=out2;
            trim_alt=out1;
        else
            trim=out1;
            trim_alt=out2;
        end
    elseif min(Jr1)>0 & min(Jr2)==0
        if min(Jr1)-min(J2)>3
            trim=out2;
            trim_alt=out1;
        else
            trim=out1;
            trim_alt=out2;
        end
    elseif min(Jr2)>0 & min(Jr1)==0
        if min(Jr2)-min(J1)>3
            trim=out1;
            trim_alt=out2;
        else
            trim=out2;
            trim_alt=out1;
        end
    else
        if min(J2)>0 & min(J1)>0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       2/3/2016
%   
%             o1=abs(length(dis1_1)-max((1-x_first)*stats_g.BoundingBox(3),x_first*stats_g.BoundingBox(4)));
%             o2=abs(length(dis2_1)-max((1-x_first)*stats_g.BoundingBox(3),x_first*stats_g.BoundingBox(4)));
%             if abs(o1-o2)<16 | max(dis1_1)/max(dis2_1)<2 | max(dis2_1)/max(dis1_1)<2
%                 if min(J2)>min(J1)
%                     trim=out1;
%                     trim_alt=out2;
%                 else
%                     trim=out2;
%                     trim_alt=out1;
%                 end
%             elseif o1>o2
%                 trim=out2;
%                 trim_alt=out1;
%             elseif o2>01
%                 trim=out1;
%                 trim_alt=out2;
%             end
            [tI1, tJ1, tV1]=find(out1);
            [tI2, tJ2, tV2]=find(out2);
            mean1=mean(tV1); std1=std(tV1);% /mean1);       % normalized standard deviation
            mean2=mean(tV2); std2=std(tV2);% /mean2);       % normalized standard deviation
            if mean1*std1/sqrt(length(tV1))-mean2*std2/sqrt(length(tV2))>0                  % mean(distance) * std(distance) / total size --> the smaller, the better
                trim=out2;
                trim_alt=out1;
            else
                trim=out1;
                trim_alt=out2;
            end
%
%       2/3/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        elseif min(J2)>0 & min(J1)==0
            if inside_bone_1==1
                trim=out2;
                trim_alt=out1;
            else
                trim=out;
                trim_alt=out2;
            end
        elseif min(J1)>0 & min(J2)==0
            if inside_bone_2==1
                trim=out1;
                trim_alt=out2;
            else
                trim=out;
                trim_alt=out1;
           end
        else                
            if inside_bone_1==1 & inside_bone_2==1  % error --> just put out1 as trim
                trim=out;
                trim_alt=out2;
            elseif inside_bone_1==0 & inside_bone_2==1
                trim=out;
                trim_alt=out2;
            elseif inside_bone_1==1 & inside_bone_2==0
                trim=out;
                trim_alt=out1;
            else
                trim=out;
                trim_alt=out1;
            end
        end
    end
        
    
    
    %if abs(min(J1)-min(J2))<round((length(dis2_1)+length(dis1_1))/2)
    %    if max(dis1_1)>max(dis2_1)
    %        trim=out1;
    %    else
    %        trim=out2;
    %    end
    %else
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Modified on June 5, 2008
%
%                               compare distance between green and red if red is there
%
%%%    measure=[min(J1), min(J2), (1-isempty(find(dis1)))*min(Jr1), (1-isempty(find(dis2)))*min(Jr2)];
%%%    [I J]=min(measure(find(measure)));
%%%    switch J
%%%        case 1
%%%            trim=out1;
%%%            trim_alt=out2;                          % alternate trim candidate
%%%        case 2
%%%            trim=out2;
%%%            trim_alt=out1;                          % alternate trim candidate
%%%        case 3
%%%            trim=out1;
%%%            trim_alt=out2;                          % alternate trim candidate
%%%        case 4
%%%            trim=out2;
%%%            trim_alt=out1;                          % alternate trim candidate
%%%    end
            
%    if max(dis1)>0 | max(dis2)>0                    % if there is red
%        if min(Jr1)>min(Jr2)                        % select closer red side
%            trim=out2;
%            trim_alt=out1;                          % alternate trim candidate
%        else
%            trim=out1;
%            trim_alt=out2;                          % alternate trim candidate
%        end




%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%else                                            % if there is red
%        if min(J1)>min(J2)
%            trim=out2;
%            trim_alt=out1;                          % alternate trim candidate
%        else
%            trim=out1;
%            trim_alt=out2;                          % alternate trim candidate
%        end


%        if length(dis1_1)>length(dis2_1)
%            trim=out2;
%        else
%            trim=out1;
%        end
    %end
%    end
end
return

function [trim, trim_alt, GFP_trim]=move_7(DIC_image,stats_g,red1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Modified on June 5, 2008
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_dist=5;                            % distance differencee between (green:red1) and (green:red2)
area=sum(sum(stats_g.Image));
dx=size(DIC_image,1);
dy=size(DIC_image,2);

if nargin==3
    red=(red1>0);
    boundary=DIC_image-imerode(DIC_image,strel('disk',1));
    [L n]=bwlabel(red);
    stats_red=regionprops(L, 'Orientation', 'Centroid', 'Area'); clear L
    ang=abs([stats_red.Orientation]-stats_g.Orientation);
    ind=find(ang==min(ang));
    
    %d=sum(((cat(1,stats_red.Centroid)-(stats_g.Centroid'*ones(1,n))').^2),2);
    %d1=find(d==min(d(:)));
    ori_r=[stats_g.Orientation];
    %ori_r=stats_red(ind).Orientation;
elseif nargin==2
    red=zeros(size(DIC_image));
    boundary=DIC_image-imerode(DIC_image,strel('disk',1));
    ori_g=[stats_g.Orientation];
    ori_r=[stats_g.Orientation];
end
%ori_temp=abs(ori_g-ori_r);
%if ori_temp>90
%    ori_temp=ori_temp-90;
%end
%if ori_temp>30
%    ori=ori_g;
%else
%    ori=ori_r;
%end

out1=zeros(size(boundary)); out_test1=zeros(size(boundary));
out2=zeros(size(boundary)); out_test2=zeros(size(boundary)); trim=zeros(size(boundary));
out_r1=out1; out_r2=out2;
GFP_trim=zeros(size(boundary));

%if nargin==3
%    red=out1;
%end

ori=ori_r;

if ori<=0
    direction=(ori+90)*pi/180;
    if ori+90>45
        x_first=1;
    else
        x_first=0;
    end
elseif ori>0
    direction=(ori-90)*pi/180;
    if ori-90<-45
        x_first=1;
    else
        x_first=0;
    end
end
ratio=tan(direction);


%if nargin==3
%    start_x=round(stats_r.BoundingBox(2));
%    end_x=start_x+stats_r.BoundingBox(4)-1;
%    start_y=round(stats_r.BoundingBox(1));
%    end_y=start_y+stats_r.BoundingBox(3)-1;
%    red(start_x:end_x,start_y:end_y)=stats_r.Image;
%end

template=stats_g.Image;
start_x=round(stats_g.BoundingBox(2));
end_x=start_x+stats_g.BoundingBox(4)-1;
start_y=round(stats_g.BoundingBox(1));
end_y=start_y+stats_g.BoundingBox(3)-1;
new_g=zeros(size(boundary));
new_g(start_x:end_x,start_y:end_y)=stats_g.Image;
out=boundary & new_g;

test=1;
count=0;
while test==1
    new_g=zeros(size(out1));
    count=count+1;
    if x_first==0
        st_y=start_y+count;
        ed_y=end_y+count;
        st_x=start_x-round(ratio*count);
        ed_x=end_x-round(ratio*count);
    else
        st_x=start_x+count;
        ed_x=end_x+count;
        st_y=start_y-round(count/ratio);
        ed_y=end_y-round(count/ratio);
    end
    if st_x>size(boundary,1) | ed_x>size(boundary,1) | ed_y>size(boundary,2) | st_y>size(boundary,2) | st_x<1 | st_y<1 | ed_x<1 | ed_y<1
        if count==1
            dis1_1(count)=0;
            dis1=[0 0];
            if nargin==3
                dis1=0;
            end
            break
        end
        break
    end
    new_g(st_x:ed_x,st_y:ed_y)=template;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Modified on June 5, 2008
%
    %if nargin==2
        new_g=new_g & DIC_image;
        template=new_g(st_x:ed_x,st_y:ed_y);
    %elseif nargin==3
    %    new_g=new_g & DIC_image & red;
    %end
    if (x_first==1 & count>round(dx/3)) |(x_first==0 & count>round(dy/3)) | sum(sum(new_g))/area==0 %<0.30
        if count==1
            dis1_1(count)=0;
            dis1=[0 0];
            if nargin==2
                dis1=0;
            end
        end
        break
    end
    %if sum(sum(new_g))/area<0.20 %isempty(find(new_g))==1
    %    if count==1
    %        dis1_1(count)=0;
    %        if nargin==3
    %            dis1=0;
    %        end
    %    end
    %    break
    %end
%    out1=min(out1,count*(new_g & boundary));
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    out1=out1+(1-double(out1>0)).*(count*double(new_g & boundary));
    out_r1=out_r1+(1-double(out_r1>0)).*(count*double(new_g & red));
    if nargin==3
        out_test1=out_test1+(new_g & red);
        dis1(count)=size(find((new_g & red)==1),1);
        dis1_1(count)=size(find((new_g & boundary)==1),1);
    elseif nargin==2
        %out_test1=out_test1+(new_g);
        dis1(count)=size(find((new_g & red)==1),1);
        dis1_1(count)=size(find((new_g & boundary)==1),1);
    end
end

template=stats_g.Image;
count=0;
while test==1
    new_g=zeros(size(out1));
    count=count-1;
    if x_first==0
        st_y=start_y+count;
        ed_y=end_y+count;
        st_x=start_x-round(ratio*count);
        ed_x=end_x-round(ratio*count);
    else
        st_x=start_x+count;
        ed_x=end_x+count;
        st_y=start_y-round(count/ratio);
        ed_y=end_y-round(count/ratio);
    end
    if st_x>size(boundary,1) | ed_x>size(boundary,1) | ed_y>size(boundary,2) | st_y>size(boundary,2) | st_x<1 | st_y<1 | ed_x<1 | ed_y<1
        if -count==1
%            dis1_1(count)=0;
            dis2_1(-count)=0;
            dis2(1:2)=[0 0];
            if nargin==3
                dis2=0;
            end
            break
        end
        break
    end
    new_g(st_x:ed_x,st_y:ed_y)=template;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Modified on June 5, 2008
%
    %if nargin==2
        new_g=new_g & DIC_image;
        template=new_g(st_x:ed_x,st_y:ed_y);
    %elseif nargin==3
    %    new_g=new_g & DIC_image & red;
    %end
    if (x_first==1 & -count>round(dx/3)) |(x_first==0 & -count>round(dy/3)) | sum(sum(new_g))/area==0 %<0.30
        if count==-1
            dis2_1(-count)=0;
            dis2=[0 0];
            if nargin==2
                dis2=0;
            end
        end
        break
    end
    %if sum(sum(new_g))/area<0.05 %isempty(find(new_g))==1
    %    if count==-1
    %        dis2_1(-count)=0;
    %        if nargin==3
    %            dis2=0;
    %        end
    %    end
    %    break
    %end
    
%    out2=min(out2,-count*(new_g & boundary));
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

out2=out2+(1-double(out2>0)).*(-count*double(new_g & boundary));
out_r2=out_r2+(1-double(out_r2>0)).*(-count*double(new_g & red));
    if nargin==3
        out_test2=out_test2+(new_g & red);
        dis2(-count)=size(find((new_g & red)==1),1);
        dis2_1(-count)=size(find((new_g & boundary)==1),1);
    elseif nargin==2
        %out_test2=out_test2+(new_g);
        dis2(-count)=size(find((new_g & red)==1),1);
        dis2_1(-count)=size(find((new_g & boundary)==1),1);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Modified on June 5, 2008
%
if nargin==4                % originally it was 'if nargin==3'
%
%                           % bypass this routine
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    common_1=out1 & red1;
    test_out_1=double(common_1).*out1;
    test_red_1=double(common_1).*red1;
%    if isempty(find(test_out_1==1))==1              % Matlab 6.5
    if isempty(find(test_out_1,1))==1              % Matlab 7.1
        ave_dis_out_1=sum(test_out_1(:));
    else
        ave_dis_out_1=sum(test_out_1(:))/length(find(test_out_1));
    end
%    if isempty(find(test_red_1==1))==1               % Matlab 6.5
    if isempty(find(test_red_1,1))==1               % Matlab 7.1
        ave_dis_red_1=sum(test_red_1(:));
    else
        ave_dis_red_1=sum(test_red_1(:))/length(find(test_red_1));
    end
    common_2=out2 & red1;
    test_out_2=double(common_2).*out2;
    test_red_2=double(common_2).*red1;
%    if isempty(find(test_out_2==1))==1               % Matlab 6.5
    if isempty(find(test_out_2,1))==1               % Matlab 7.1
        ave_dis_out_2=sum(test_out_2(:));
    else
        ave_dis_out_2=sum(test_out_2(:))/length(find(test_out_2));
    end
%    if isempty(find(test_red_2==1))==1               % Matlab 6.5
    if isempty(find(test_red_2,1))==1               % Matlab 7.1
        ave_dis_red_2=sum(test_red_2(:));
    else
        ave_dis_red_2=sum(test_red_2(:))/length(find(test_red_2));
    end
    
    if ave_dis_out_1~=0 & sum(test_out_1(:))>sum(test_out_2(:))
        if ave_dis_out_1>ave_dis_red_1
            trim=out1;
            trim_alt=out2;                          % alternate trim candidate
        else
            trim=zeros(size(out1));
            trim_alt=zeros(size(out1));
            GFP_trim=out1;
        end
    elseif ave_dis_out_2~=0 & sum(test_out_2(:))>sum(test_out_1(:))
        if ave_dis_out_2>ave_dis_red_2
            trim=out2;
            trim_alt=out1;                          % alternate trim candidate
        else
            trim=zeros(size(out1));
            trim_alt=zeros(size(out1));
            GFP_trim=out2;
        end
    else
        trim=zeros(size(out1));
        trim_alt=zeros(size(out1));
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Modified on June 5, 2008
%
elseif nargin==2 | nargin==3    % originally it was 'elseif nargin==2'
%
%                               % all process goes this routine
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [I1 J1]=find(dis1_1==max(dis1_1));
    [I2 J2]=find(dis2_1==max(dis2_1));
    if length(dis1)<2
        dis1=[dis1, zeros(2-length(dis1))];
    end
    if length(dis2)<2
        dis2=[dis2, zeros(2-length(dis2))];
    end
            
    [Ir1 Jr1]=find(dis1(1:end)==max(dis1(1:end)));    %   added on June 5, 2008
    [Ir2 Jr2]=find(dis2(1:end)==max(dis2(1:end)));    %   added on June 5, 2008
    if max(dis1(1:end))>0 & max(dis2(1:end))>0
        if max(dis2(1:end))/max(dis1(1:end))>1.5
            Jr1=0;
        elseif max(dis1(1:end))/max(dis2(1:end))>1.5
            Jr2=0;
        end
    end
    if max(dis1(1:end))==0
        Jr1=0;
    elseif max(dis2(1:end))==0
        Jr2=0;
    end
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Jun3 12, 2008
%
%       Override if shape of green and bone surface similar
%
    if max(dis1_1)>0 & max(dis2_1)>0
        if max(dis2_1)/max(dis1_1)>1.5 & min(J1)/min(J2)>1.5 & min(min(J1), min(J2))>5
            J1=10000;
        elseif max(dis1_1)/max(dis2_1)>1.5 & min(J2)/min(J1)>1.5 & min(min(J1), min(J2))>5
            J2=10000;
        end
    end
%
    %if find(cumsum(dis1_1),1)/find(cumsum(dis2_1),1)>3
    %    J1=10000;
    %elseif find(cumsum(dis2_1),1)/find(cumsum(dis1_1),1)>3
    %    J2=10000;
    %end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       2/3/2016
%
%     if min(Jr1)>0 & min(Jr2)>0
%         if min(Jr1)>min(Jr2)
%             trim=out2;
%             trim_alt=out1;
%         else
%             trim=out1;
%             trim_alt=out2;
%         end
%     elseif min(Jr1)>0 & min(Jr2)==0
%         if min(Jr1)-min(J2)>test_dist
%             trim=out2;
%             trim_alt=out1;
%         else
%             trim=out1;
%             trim_alt=out2;
%         end
%     elseif min(Jr2)>0 & min(Jr1)==0
%         if min(Jr2)-min(J1)>test_dist
%             trim=out1;
%             trim_alt=out2;
%         else
%             trim=out2;
%             trim_alt=out1;
%         end
%     else
%         if min(J2)>min(J1)
%             trim=out1;
%             trim_alt=out2;
%         else
%             trim=out2;
%             trim_alt=out1;
%         end
%     end
%
        [tI1, tJ1, tVr1]=find(out_r1);
        [tI2, tJ2, tVr2]=find(out_r2);
        [tI1, tJ1, tV1]=find(out1);
        [tI2, tJ2, tV2]=find(out2);
        if isempty(tVr1)
            Mean(3)=1000; Std(3)=1000;
            Md(3)=1000;
        else
%             Mean(3)=mean(tVr1); Std(3)=max(std(tVr1/Mean(3)),0.001);       % normalized standard deviation
            Mean(3)=mean(tVr1); Std(3)=max(std(tVr1),0.001);       % normalized standard deviation
            Md(3)=Mean(3)*Std(3)/sqrt(length(tVr1));
        end
        if isempty(tVr2)
            Mean(4)=1000; Std(4)=1000;
            Md(4)=1000;
        else
%             Mean(4)=mean(tVr2); Std(4)=max(std(tVr2/Mean(4)),0.001);       % normalized standard deviation
            Mean(4)=mean(tVr2); Std(4)=max(std(tVr2),0.001);       % normalized standard deviation
            Md(4)=Mean(4)*Std(4)/sqrt(length(tVr2));
        end
        if isempty(tV1)
            Mean(1)=1000; Std(1)=1000;
            Md(1)=1000;
        else
%             Mean(1)=mean(tV1); Std(1)=max(std(tV1/Mean(1)),0.001);       % normalized standard deviation
            Mean(1)=mean(tV1); Std(1)=max(std(tV1),0.001);       % normalized standard deviation
            Md(1)=Mean(1)*Std(1)/sqrt(length(tV1));
        end
        if isempty(tV2)
            Mean(2)=1000; Std(2)=1000;
            Md(2)=1000;
        else
%             Mean(2)=mean(tV2); Std(2)=max(std(tV2/Mean(2)),0.001);       % normalized standard deviation
            Mean(2)=mean(tV2); Std(2)=max(std(tV2),0.001);       % normalized standard deviation
            Md(2)=Mean(2)*Std(2)/sqrt(length(tV2));
        end
        
%         Md=[Md1 Md2 Mdr1 Mdr2];
        ind=find(Md==min(Md));
        if ~isempty(ind)
            if length(ind)>1
                Mdd=Mean(ind);
                ind1=find(Mdd==min(Mdd));
                ind=ind(ind1);
                if length(ind)>1
                    ind=ind(1);
                end
            end
            switch ind
                case 1
                    trim=out1;
                    trim_alt=out2;
                case 2
                    trim=out2;
                    trim_alt=out1;                
                case 3
                    trim=out1;
                    trim_alt=out2;
                case 4
                    trim=out2;
                    trim_alt=out1; 
            end
        else
            trim=out1;
            trim_alt=out2;
        end
%
%       2/3/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Jun3 12, 2008
%
%       Override if shape of green and bone surface similar
%
    if J1==10000
        trim=out2;
        trim_alt=out1;
    elseif J2==10000
        trim=out1;
        trim_alt=out2;
    end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




    
    %if abs(min(J1)-min(J2))<round((length(dis2_1)+length(dis1_1))/2)
    %    if max(dis1_1)>max(dis2_1)
    %        trim=out1;
    %    else
    %        trim=out2;
    %    end
    %else
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Modified on June 5, 2008
%
%                               compare distance between green and red if red is there
%
%%%    measure=[min(J1), min(J2), (1-isempty(find(dis1)))*min(Jr1), (1-isempty(find(dis2)))*min(Jr2)];
%%%    [I J]=min(measure(find(measure)));
%%%    switch J
%%%        case 1
%%%            trim=out1;
%%%            trim_alt=out2;                          % alternate trim candidate
%%%        case 2
%%%            trim=out2;
%%%            trim_alt=out1;                          % alternate trim candidate
%%%        case 3
%%%            trim=out1;
%%%            trim_alt=out2;                          % alternate trim candidate
%%%        case 4
%%%            trim=out2;
%%%            trim_alt=out1;                          % alternate trim candidate
%%%    end
            
%    if max(dis1)>0 | max(dis2)>0                    % if there is red
%        if min(Jr1)>min(Jr2)                        % select closer red side
%            trim=out2;
%            trim_alt=out1;                          % alternate trim candidate
%        else
%            trim=out1;
%            trim_alt=out2;                          % alternate trim candidate
%        end




%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%else                                            % if there is red
%        if min(J1)>min(J2)
%            trim=out2;
%            trim_alt=out1;                          % alternate trim candidate
%        else
%            trim=out1;
%            trim_alt=out2;                          % alternate trim candidate
%        end


%        if length(dis1_1)>length(dis2_1)
%            trim=out2;
%        else
%            trim=out1;
%        end
    %end
%    end
end
return

function [trim, trim_alt, GFP_trim]=move_5(DIC_image,stats_g,red1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Modified on June 5, 2008
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

area=sum(sum(stats_g.Image));
dx=size(DIC_image,1);
dy=size(DIC_image,2);

if nargin==3
    red=(red1>0);
    boundary=DIC_image-imerode(DIC_image,strel('disk',1));
    [L n]=bwlabel(red);
    stats_red=regionprops(L, 'Orientation', 'Centroid', 'Area'); clear L
    ang=abs([stats_red.Orientation]-stats_g.Orientation);
    ind=find(ang==min(ang));
    
    %d=sum(((cat(1,stats_red.Centroid)-(stats_g.Centroid'*ones(1,n))').^2),2);
    %d1=find(d==min(d(:)));
    ori_r=[stats_g.Orientation];
    %ori_r=stats_red(ind).Orientation;
elseif nargin==2
    red=zeros(size(DIC_image));
    boundary=DIC_image-imerode(DIC_image,strel('disk',1));
    ori_g=[stats_g.Orientation];
    ori_r=[stats_g.Orientation];
end
%ori_temp=abs(ori_g-ori_r);
%if ori_temp>90
%    ori_temp=ori_temp-90;
%end
%if ori_temp>30
%    ori=ori_g;
%else
%    ori=ori_r;
%end

out1=zeros(size(boundary)); out_test1=zeros(size(boundary));
out2=zeros(size(boundary)); out_test2=zeros(size(boundary)); trim=zeros(size(boundary));
GFP_trim=zeros(size(boundary));

%if nargin==3
%    red=out1;
%end

ori=ori_r;

if ori<=0
    direction=(ori+90)*pi/180;
    if ori+90>45
        x_first=1;
    else
        x_first=0;
    end
elseif ori>0
    direction=(ori-90)*pi/180;
    if ori-90<-45
        x_first=1;
    else
        x_first=0;
    end
end
ratio=tan(direction);


%if nargin==3
%    start_x=round(stats_r.BoundingBox(2));
%    end_x=start_x+stats_r.BoundingBox(4)-1;
%    start_y=round(stats_r.BoundingBox(1));
%    end_y=start_y+stats_r.BoundingBox(3)-1;
%    red(start_x:end_x,start_y:end_y)=stats_r.Image;
%end

start_x=round(stats_g.BoundingBox(2));
end_x=start_x+stats_g.BoundingBox(4)-1;
start_y=round(stats_g.BoundingBox(1));
end_y=start_y+stats_g.BoundingBox(3)-1;

test=1;
count=0;
while test==1
    new_g=zeros(size(out1));
    count=count+1;
    if x_first==0
        st_y=start_y+count;
        ed_y=end_y+count;
        st_x=start_x-round(ratio*count);
        ed_x=end_x-round(ratio*count);
    else
        st_x=start_x+count;
        ed_x=end_x+count;
        st_y=start_y-round(count/ratio);
        ed_y=end_y-round(count/ratio);
    end
    if st_x>size(boundary,1) | ed_x>size(boundary,1) | ed_y>size(boundary,2) | st_y>size(boundary,2) | st_x<1 | st_y<1 | ed_x<1 | ed_y<1
        if count==1
            dis1_1(count)=0;
            dis1(count)=0;
            if nargin==3
                dis1=0;
            end
            break
        end
        break
    end
    new_g(st_x:ed_x,st_y:ed_y)=stats_g.Image;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Modified on June 5, 2008
%
    %if nargin==2
        new_g=new_g & DIC_image;
    %elseif nargin==3
    %    new_g=new_g & DIC_image & red;
    %end
    if (x_first==1 & count>round(dx/3)) |(x_first==0 & count>round(dy/3)) | sum(sum(new_g))/area==0
        break
    end
    %if sum(sum(new_g))/area<0.20 %isempty(find(new_g))==1
    %    if count==1
    %        dis1_1(count)=0;
    %        if nargin==3
    %            dis1=0;
    %        end
    %    end
    %    break
    %end
%    out1=min(out1,count*(new_g & boundary));
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    out1=out1+(1-double(out1>0)).*(count*double(new_g & boundary));
    if nargin==3
        out_test1=out_test1+(new_g & red);
        dis1(count)=size(find((new_g & red)==1),1);
        dis1_1(count)=size(find((new_g & boundary)==1),1);
    elseif nargin==2
        %out_test1=out_test1+(new_g);
        dis1(count)=size(find((new_g & red)==1),1);
        dis1_1(count)=size(find((new_g & boundary)==1),1);
    end
end

count=0;
while test==1
    new_g=zeros(size(out1));
    count=count-1;
    if x_first==0
        st_y=start_y+count;
        ed_y=end_y+count;
        st_x=start_x-round(ratio*count);
        ed_x=end_x-round(ratio*count);
    else
        st_x=start_x+count;
        ed_x=end_x+count;
        st_y=start_y-round(count/ratio);
        ed_y=end_y-round(count/ratio);
    end
    if st_x>size(boundary,1) | ed_x>size(boundary,1) | ed_y>size(boundary,2) | st_y>size(boundary,2) | st_x<1 | st_y<1 | ed_x<1 | ed_y<1
        if -count==1
%            dis1_1(count)=0;
            dis2_1(-count)=0;
            dis2(-count)=0;
            if nargin==3
                dis2=0;
            end
            break
        end
        break
    end
    new_g(st_x:ed_x,st_y:ed_y)=stats_g.Image;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Modified on June 5, 2008
%
    %if nargin==2
        new_g=new_g & DIC_image;
    %elseif nargin==3
    %    new_g=new_g & DIC_image & red;
    %end
    if (x_first==1 & count>round(dx/3)) |(x_first==0 & count>round(dy/3)) | sum(sum(new_g))/area==0
        break
    end
    %if sum(sum(new_g))/area<0.05 %isempty(find(new_g))==1
    %    if count==-1
    %        dis2_1(-count)=0;
    %        if nargin==3
    %            dis2=0;
    %        end
    %    end
    %    break
    %end
    
%    out2=min(out2,-count*(new_g & boundary));
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

out2=out2+(1-double(out2>0)).*(-count*double(new_g & boundary));
    if nargin==3
        out_test2=out_test2+(new_g & red);
        dis2(-count)=size(find((new_g & red)==1),1);
        dis2_1(-count)=size(find((new_g & boundary)==1),1);
    elseif nargin==2
        %out_test2=out_test2+(new_g);
        dis2(-count)=size(find((new_g & red)==1),1);
        dis2_1(-count)=size(find((new_g & boundary)==1),1);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Added February 28 2013
if ~exist('dis2_1')
    trim=out1;
    trim_alt=zeros(size(out1));
    return
end
if ~exist('dis1_1')
    trim=out2;
    trim_alt=zeros(size(out2));
    return
end
%
%   Added February 28 2013
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Modified on June 5, 2008
%
if nargin==4                % originally it was 'if nargin==3'
%
%                           % bypass this routine
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    common_1=out1 & red1;
    test_out_1=double(common_1).*out1;
    test_red_1=double(common_1).*red1;
%    if isempty(find(test_out_1==1))==1              % Matlab 6.5
    if isempty(find(test_out_1,1))==1              % Matlab 7.1
        ave_dis_out_1=sum(test_out_1(:));
    else
        ave_dis_out_1=sum(test_out_1(:))/length(find(test_out_1));
    end
%    if isempty(find(test_red_1==1))==1               % Matlab 6.5
    if isempty(find(test_red_1,1))==1               % Matlab 7.1
        ave_dis_red_1=sum(test_red_1(:));
    else
        ave_dis_red_1=sum(test_red_1(:))/length(find(test_red_1));
    end
    common_2=out2 & red1;
    test_out_2=double(common_2).*out2;
    test_red_2=double(common_2).*red1;
%    if isempty(find(test_out_2==1))==1               % Matlab 6.5
    if isempty(find(test_out_2,1))==1               % Matlab 7.1
        ave_dis_out_2=sum(test_out_2(:));
    else
        ave_dis_out_2=sum(test_out_2(:))/length(find(test_out_2));
    end
%    if isempty(find(test_red_2==1))==1               % Matlab 6.5
    if isempty(find(test_red_2,1))==1               % Matlab 7.1
        ave_dis_red_2=sum(test_red_2(:));
    else
        ave_dis_red_2=sum(test_red_2(:))/length(find(test_red_2));
    end
    
    if ave_dis_out_1~=0 & sum(test_out_1(:))>sum(test_out_2(:))
        if ave_dis_out_1>ave_dis_red_1
            trim=out1;
            trim_alt=out2;                          % alternate trim candidate
        else
            trim=zeros(size(out1));
            trim_alt=zeros(size(out1));
            GFP_trim=out1;
        end
    elseif ave_dis_out_2~=0 & sum(test_out_2(:))>sum(test_out_1(:))
        if ave_dis_out_2>ave_dis_red_2
            trim=out2;
            trim_alt=out1;                          % alternate trim candidate
        else
            trim=zeros(size(out1));
            trim_alt=zeros(size(out1));
            GFP_trim=out2;
        end
    else
        trim=zeros(size(out1));
        trim_alt=zeros(size(out1));
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Modified on June 5, 2008
%
elseif nargin==2 | nargin==3    % originally it was 'elseif nargin==2'
%
%                               % all process goes this routine
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [I1 J1]=find(dis1_1==max(dis1_1));
    [I2 J2]=find(dis2_1==max(dis2_1));
    [Ir1 Jr1]=find(dis1(2:end)==max(dis1(2:end)));    %   added on June 5, 2008
    [Ir2 Jr2]=find(dis2(2:end)==max(dis2(2:end)));    %   added on June 5, 2008
    
    if max(dis1(2:end))>0 & max(dis2(2:end))>0
        if max(dis2(2:end))/max(dis1(2:end))>2
            Jr1=0;
        elseif max(dis1(2:end))/max(dis2(2:end))>2
            Jr2=0;
        end
    end

    %if abs(min(J1)-min(J2))<round((length(dis2_1)+length(dis1_1))/2)
    %    if max(dis1_1)>max(dis2_1)
    %        trim=out1;
    %    else
    %        trim=out2;
    %    end
    %else

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%   Modified on February 28, 2013    
if length(dis1_1)>length(dis2_1) &  max(dis1_1)<max(dis2_1)
    trim=out2;
    trim_alt=out1;
elseif length(dis1_1)<length(dis2_1) & max(dis1_1)>max(dis2_1)
        trim=out1;
        trim_alt=out2;
else
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% %   Modified on June 5, 2008
% %
% %                               compare distance between green and red if red is there
% %
    measure=[min(J1), min(J2), (1-isempty(find(dis1)))*min(Jr1), (1-isempty(find(dis2)))*min(Jr2)];
    [I J]=min(measure(find(measure)));
    switch J
        case 1
            trim=out1;
            trim_alt=out2;                          % alternate trim candidate
        case 2
            trim=out2;
            trim_alt=out1;                          % alternate trim candidate
        case 3
            trim=out1;
            trim_alt=out2;                          % alternate trim candidate
        case 4
            trim=out2;
            trim_alt=out1;                          % alternate trim candidate
    end
end            
% %    if max(dis1)>0 | max(dis2)>0                    % if there is red
% %        if min(Jr1)>min(Jr2)                        % select closer red side
% %            trim=out2;
% %            trim_alt=out1;                          % alternate trim candidate
% %        else
% %            trim=out1;
% %            trim_alt=out2;                          % alternate trim candidate
% %        end
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        
        
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%else                                            % if there is red
%        if min(J1)>min(J2)
%            trim=out2;
%            trim_alt=out1;                          % alternate trim candidate
%        else
%            trim=out1;
%            trim_alt=out2;                          % alternate trim candidate
%        end


%        if length(dis1_1)>length(dis2_1)
%            trim=out2;
%        else
%            trim=out1;
%        end
    %end
%    end
end
return

function [trim1]=move1_2(green1, red1, MARn)

col=green1.trim.col;
row=green1.trim.row;
val=green1.trim.val;
ver=green1.trim.size_trim_ver;
hor=green1.trim.size_trim_hor;
out1=full(sparse(col,row,val,ver,hor));

col=green1.trim_alt.col;
row=green1.trim_alt.row;
val=green1.trim_alt.val;
ver=green1.trim_alt.size_trim_ver;
hor=green1.trim_alt.size_trim_hor;
out2=full(sparse(col,row,val,ver,hor));
%out1=green1;
[I1 J1 V1]=find(out1);
[I2 J2 V2]=find(out2);

if nargin==3
    common_1=out1 & red1;
    test_out_1=double(common_1).*out1;
    test_red_1=double(common_1).*red1;
%    if isempty(find(test_out_1==1))==1              % Matlab 6.5
    if isempty(find(test_out_1,1))==1              % Matlab 7.1
        ave_dis_out_1=sum(test_out_1(:));
    else
        ave_dis_out_1=sum(test_out_1(:))/length(find(test_out_1));
    end
%    if isempty(find(test_red_1==1))==1               % Matlab 6.5
    if isempty(find(test_red_1,1))==1               % Matlab 7.1
        ave_dis_red_1=sum(test_red_1(:));
    else
        ave_dis_red_1=sum(test_red_1(:))/length(find(test_red_1));
    end
    common_2=out2 & red1;
    test_out_2=double(common_2).*out2;
    test_red_2=double(common_2).*red1;
%    if isempty(find(test_out_2==1))==1               % Matlab 6.5
    if isempty(find(test_out_2,1))==1               % Matlab 7.1
        ave_dis_out_2=sum(test_out_2(:));
    else
        ave_dis_out_2=sum(test_out_2(:))/length(find(test_out_2));
    end
%    if isempty(find(test_red_2==1))==1               % Matlab 6.5
    if isempty(find(test_red_2,1))==1               % Matlab 7.1
        ave_dis_red_2=sum(test_red_2(:));
    else
        ave_dis_red_2=sum(test_red_2(:))/length(find(test_red_2));
    end
    
    if ave_dis_red_1>0 & ave_dis_red_2>0
        %mar1=ave_dis_out_1-ave_dis_red_1;
        %mar2=ave_dis_out_2-ave_dis_red_2;
        %if mar1>MARn
        %    if mar2>MARn
        %        trim=out1;
        %    else
        %        trim=out2;
        %    end
        %elseif mar1<MARn
        %    trim=out1;
        %end
        if min(V1)>MARn
            if min(V2)>MARn
                trim1=out1;
            else
                trim1=out2;
            end
        %elseif min(V1)<MARn
        else
            trim1=out1;
        end
    elseif ave_dis_red_1>0 & ave_dis_red_2==0
        %mar1=ave_dis_out_1-ave_dis_red_1;
        %if mar1>MARn
        %    trim=out2;
        %else
        %    trim=out1;
        %end
        if min(V1)>MARn
            trim1=out2;
        else
            trim1=out1;
        end
    %elseif ave_dis_red_2>0 & ave_dis_red_1==0
    %    mar2=ave_dis_out_2-ave_dis_red_2;
    %    if mar2>MARn
    %        trim=out1;
    %    else
    %        trim=out2;
    %    end
    else%if ave_dis_red_1==0 & ave_dis_red_2==0
        trim1=out1;
    end
        
    %if ave_dis_out_1~=0 & ave_dis_out1>ave_out_dis2  % sum(test_out_1(:))>sum(test_out_2(:))   June, 10,2008
    %    if ave_dis_out_1>ave_dis_red_1 & ave_dis_out_1<MARn
    %        trim=out1;
    %    else
    %        if ave_dis_out_2>ave_dis_out_1
    %            trim=out1;
    %        else
    %            trim=out2;
    %        end
    %    end
    %elseif ave_dis_out_2~=0 & ave_dis_out2>ave_out_dis1  % sum(test_out_2(:))>sum(test_out_1(:))   June, 10,2008
    %    if ave_dis_out_2>ave_dis_red_2 & ave_dis_out_2<MARn
    %        trim=out2;
    %    else
    %        if ave_dis_out_1>ave_dis_out_2
    %            trim=out2;
    %        else
    %            trim=out1;
    %        end
    %    end
    %else
%        trim=zeros(size(out1));
    %    trim=zeros(size(out2));
    %end
elseif nargin==1
    %%%%%%%%%%%%%%%%%%%%
    %
    %   3/25
    %
    if isempty(find(out1,1))==1              % Matlab 7.1
        ave_dis_out_1=sum(out1(:));
    else
        ave_dis_out_1=sum(out1(:))/length(find(out1));
    end
    if isempty(find(out2,1))==1              % Matlab 7.1
        ave_dis_out_2=sum(out2(:));
    else
        ave_dis_out_2=sum(out2(:))/length(find(out2));
    end

    if ave_dis_out_1 > ave_dis_out_2
        trim1=out2;
    else
        trim1=out1;
    end
    %
    %   3/25
    %
    %%%%%%%%%%%%%%%%%%%%


elseif nargin==2 %| nargin==3
    %['use ''[output]=move1_2(green_label1, green_alternative_label1, red_label, standard_deviation_of_MAR);''']
    trim1=out1;
end
return


function [new_1]=find_min(old, new)

template = double(old & new);

if isempty(find(template))==1
    new_1 = old + new;
    %new_1 = new_1 + new .* (1-template);
    return
else
    old_template = old .* template;
    new_template = new .* template;
    
    min_old_new = min(old_template, new_template);
    
    new_1 = (old .* (1-template)) + min_old_new;
    new_1 = new_1 + new .* (1-template);
    return
end

function [label]=inside_trim_label(ref, label, loop_count, exact_boundary)

temp=imerode(ref,strel('disk',loop_count));
    
for loop=loop_count:-1:1                            % from outside of the bone to inside
    if exact_boundary==0
        a=(ref)-imerode(ref,strel('disk',loop));% Mar/21
        label=label & a;% Mar/21
    else
        temp_1=imerode(ref,strel('disk',loop-1));
        a=imsubtract(temp,temp_1);
        temp3=label & a;
    end
    if loop>1
        if exact_boundary==0
            label=imdilate(label,strel('disk',1));% Mar/21
        else
            label=label | imdilate(temp3,strel('disk',1));
        end
    end
    if exact_boundary==1
        temp=temp_1;
    end
end
clear temp temp_1
clear a
if exact_boundary==1
    label=temp3;clear temp3;
end
return

function [label]=outside_trim_label(ref, label, loop_count, exact_boundary)

temp=imdilate(ref,strel('disk',loop_count));

for loop=loop_count:-1:0                            % from outside of the bone to green label
    if loop>0
        if exact_boundary==0
            a=imdilate(ref,strel('disk',loop));%-imdilate(a4,strel('disk',loop-1));% Mar/21
        else
            temp_1=imdilate(ref,strel('disk',loop-1));
            a=imsubtract(temp,temp_1);
        end
    else
        a=(ref)-imerode((ref),strel('disk',1));
    end
    if exact_boundary==0
        label=label & a;% Mar/21
    else
        temp3=label & a;
    end
    if loop>0
        if exact_boundary==0
            label=imdilate(label,strel('disk',1));% Mar/21
        else
            label=label | imdilate(temp3,strel('disk',1));
        end
    end
    if exact_boundary==1
        temp=temp_1;
    end
end
clear temp temp_1
if exact_boundary==1
    label=temp3; clear temp3;
end
return

function [x y min_d]=find_closest_points(in1, in2)

[x1 y1]=find(in1);
[x2 y2]=find(in2);
min_d=10000000;
for i=1:length(x1)
    for j=1:length(x2)
        d=(x1(i)-x2(j))^2+(y1(i)-y2(j))^2;
        if d<min_d
            min_d=d;
            x(1)=x1(i);
            x(2)=x2(j);
            y(1)=y1(i);
            y(2)=y2(j);
        end
    end
end
if isempty(x)
    x(1)=1;x(2)=1;y(1)=2;y(2)=1;
end
return

function [test5, test6]=get_rid_of_bubbles(a, threshold)

[L n]=bwlabel((imdilate((double(a>threshold)),strel('disk',3))));
st=regionprops(L,'Area','Image','FilledImage','BoundingBox','Orientation'); clear L
test=zeros(size(a));
for i=1:n
    start_x=round(st(i).BoundingBox(2));
    start_y=round(st(i).BoundingBox(1));
    end_x=start_x+st(i).BoundingBox(4)-1;
    end_y=start_y+st(i).BoundingBox(3)-1;
    if length(find(st(i).FilledImage))>400
        test(start_x:end_x,start_y:end_y)=test(start_x:end_x,start_y:end_y)+st(i).Image;
    end
end
%figure;imshow(test)

[L n]=bwlabel(test);
st1=regionprops(L,'Area','Image','FilledImage','BoundingBox','Orientation','Perimeter'); clear L
get_rid=[];
for i=1:n
    start_x=round(st1(i).BoundingBox(2));
    start_y=round(st1(i).BoundingBox(1));
    end_x=start_x+st1(i).BoundingBox(4)-1;
    end_y=start_y+st1(i).BoundingBox(3)-1;
    temp(i)=st1(i).BoundingBox(4)/st1(i).BoundingBox(3);
    temp1(i)=length(find(st1(i).Image))/prod(size(st1(i).Image)); %(st1(i).BoundingBox(3)*st1(i).BoundingBox(4));
    temp2(i)=max(st1(i).BoundingBox(4),st1(i).BoundingBox(3));
    tt(i)=length(find(st1(i).Image))/length(find(st1(i).FilledImage));
    test_circle(i)=length(find(st1(i).FilledImage))*4*pi/(st1(i).Perimeter)^2;
    if (tt(i)>0.9) 
        get_rid=[get_rid;i];
    end
end

get_rid1=[];
for i=1:length(get_rid)
    if (temp(get_rid(i))<2 & temp(get_rid(i))>1/2)      % find long shape
        get_rid1=[get_rid1; get_rid(i)];                % test - get_rid1
    end
end
get_rid2=[];
for i=1:length(get_rid1)
    if temp1(get_rid1(i))>0.2                           % find narrow shape in a big area
        get_rid2=[get_rid2; get_rid1(i)];               % test - get_rid2
    end
end
get_rid3=[];
for i=1:length(get_rid2)
    if test_circle(get_rid2(i))<0.5                     % find circular shape
        get_rid3=[get_rid3; get_rid2(i)];               % test - get_rid3
    end
end
get_rid4=[];
for i=1:length(get_rid3)
    if temp2(get_rid3(i))<50                     % find circular shape
        get_rid4=[get_rid4; get_rid3(i)];               % test - get_rid3
    end
end

test1=zeros(size(a));
for i=1:length(get_rid4)
    start_x=round(st1(get_rid4(i)).BoundingBox(2));
    start_y=round(st1(get_rid4(i)).BoundingBox(1));
    end_x=start_x+st1(get_rid4(i)).BoundingBox(4)-1;
    end_y=start_y+st1(get_rid4(i)).BoundingBox(3)-1;
    test1(start_x:end_x,start_y:end_y)=test1(start_x:end_x,start_y:end_y)+st1(get_rid4(i)).Image;
end

test2=test-test1;
%figure;imshow(test2)

test3=a.*double(imdilate(test2,strel('disk',1))-test2);
%figure;imshow(test3)
test4=imdilate(test3,strel('disk',10));
test5=a-a.*test2+double(test4.*test2)/1.5;
%figure;imshow(test5/1)
test6=imerode(test2, strel('disk',3));

return


function select_move_images(varargin)
%function select_move_images(direc, d, phr1, samples_per_bone, bone_n1, bone_n2, bone_n3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       Read individual images and Save combined image
%
%direc='C:\Users\shhong\Desktop\C-A labe\Analysis\';

if length(varargin)<5
    error('Not enough argments');
end

direct=varargin{1};
out_directory=varargin{2};
exp_name=varargin{3};
samples_per_bone=varargin{4};
for k = 5:length(varargin)
    eval(['bone_n',num2str(k-4),' = varargin{',num2str(k),'};']);
end

no_bone_type=length(varargin)-4;

eval(['a=dir(''', direct,''');'])

file_name=['_c_pseudo1_tr1.jpg       ';
           '_c_roi.jpg               ';
           '_c_shift.jpg             ';
           '_c_trim12_tr1.jpg        ';
           '_c_pseudo1_tr1_inside.jpg';
           '_c_roi_inside.jpg        ';
           '_c_trim12_tr1_inside.jpg '];

for i=1:no_bone_type  
    eval(['bone_number=bone_n',num2str(i),';']);
    for j=1:length(bone_number)
        phr=[varargin{3}(1,:),num2str(bone_number(j)),'FL'];
        dd=[out_directory,'\',phr];
        if isdir(dd)==0
            mkdir(dd)
        end

        %eval(['copyfile(','''',direct,'analysis1_tr1.xls',''', ''',out_directory,'\analysis1_tr1.xls','''',')'])
        
        for s=1:varargin{4}
            if varargin{4}==1
                phr=[varargin{3}(i,:), num2str(bone_number(j)), exp_type,'.jpg_Files'];
            else
                phr=[varargin{3}(i,:), num2str(bone_number(j)), exp_type,'_',num2str(s),'.jpg_Files'];
            end
            phr2=[varargin{3}(i,:), num2str(bone_number(j)),exp_type]
            direc=[direct, phr,'\'];
            a=dir(direc);
            for file_number=1:7
                if varargin{4}==1
                    in_out_file=[phr2,file_name(file_number,:)];
                else
                    in_out_file=[phr2,num2str(s),file_name(file_number,:)];
                end
                in_file_name=[direc,in_out_file];
                %if isempty(regexpi([a.name],[phr,'_',num2str(s),file_name(file_number,:)]))==1
                if isempty(regexpi(in_file_name,in_out_file))==1
                %if isempty(regexpi(in_file_name,[phr2,'_',num2str(s),file_name(file_number,:)]))==1
                    break
                else
                    eval(['copyfile(','''',in_file_name,''', ''',dd,'\',in_out_file,'''',')'])
                end
                
            end
        end
    end
end
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