function find_thresholding_values_KOMP_test_for_BioHead_1_batch(direct, exp_name, exp_type, labels, GFP, bt, bn, section, delimeter)

% direct='Z:\Histomorphometry\His31   (sent 01-10-14)\01_Submitted\Layers\';
% exp_name={'His31_C';'His31_T'};
% thresh_growth=0;  % find threshold + find growth_plate (default)
%               1;  % find threshold only
%               2;  % find growth_plate only
% labels={'green','green';'red','red'};      % first label, second label
% exp_type='F';
% GFP=0;            % no GFP
%     1;            % there is GFP

close all
thresh_growth=0;

% direct=cell2mat(varargin(1));
if strcmp(direct(end),delimeter)~=1
    direct=[direct,delimeter];
end
direc=[direct,'Images',delimeter];
out_direc=[direct,'Threshold',delimeter];
if isempty(dir(out_direc))
    mkdir(out_direc)
end

eval(['load ''',direct, 'info'''])
sample_per_bone=0;
no_bone_type=length(exp_name);
for i=1:no_bone_type
    eval(['test=isempty(info',num2str(i),');'])
    if test==1
        continue
    end
    eval(['temp=info',num2str(i),'.im;'])
%     temp1=sum(temp,2);
%     temp1=sum(cell2mat(temp(:,2:end)),2);
%     bone_no=find(temp1);
    bone_no=temp(:,1);
    sample_per_bone(i)=max(max(sample_per_bone,size(temp,2)-1));
    eval(['bone_number',num2str(i),'=bone_no;'])
    if strcmp(exp_name{i}(end),'_')==0
        exp_name{i}=[cell2mat(exp_name(i))];
    end
end


for i=1:no_bone_type
    eval(['temp=exp_name{',num2str(i),'};'])
%     if strcmp(temp(end),'_')~=1
%         temp=[temp,'_'];
%     end
    eval(['p',num2str(i),'=temp;'])
end

if thresh_growth==0 | thresh_growth==1
%     matlabpool
    

    for i=bt
        eval(['information=info',num2str(i),';'])
        eval(['bone_number=bone_number',num2str(i),';']);
        eval(['phr1=p',num2str(i),';'])
        eval(['a=dir(''',direct, 'threshold_values',num2str(i),'.mat'');']);
        if isempty(a)==0
            eval(['load ''',direct, 'threshold_values',num2str(i),''''])
        else
%             eval(['threshold_values',num2str(i),'(',bone_number{end},',',num2str(sample_per_bone(i)),',9)=0;'])
            threshold_value=information.im;
            for column=1:length(threshold_value)
                for row=2:size(information.im,2)
                    threshold_value{column,row}=0;
                end
            end
            eval(['threshold_values',num2str(i),'=threshold_value;'])
        end
        eval(['threshold=threshold_values',num2str(i),';'])

        for j=bn
%                 parfor j=1:length(bone_number)
            for s=section
                
                'ImageToolbox License test start'
                [l m]=license('checkout','image_Toolbox');
                while (~l)
                   pause(1);
                [l m]=license('checkout','image_Toolbox');
                end
                'ImageToolbox License test end'
                
                k=1;
                for ii=1:length(threshold)

                    if strcmp(threshold{ii,1},bone_number{bn})==1
                        break
                    else
                        k=k+1;
                    end
                end
                tt=threshold{k,s+1};
                if tt>0
                    continue
                end
                for loop=3%[1 2 3 4 6 7]%channel            % loop : 1->green, 2->red, 3->DIC, 4->TRAP, 6->AP, 7->DAPI, 9-->GFP(CFP)
%                 for loop=[7]%channel            % loop : 1->green, 2->red, 3->DIC, 4->TRAP, 6->AP, 7->DAPI, 9-->GFP(CFP)
                    count_trap(j,s,i)=0;
                    count_ap(j,s,i)=0;
                    close all
%                     [bone_number{j},',  ',num2str(s)]

        %             phr3=[phr1, num2str(bone_number(j)), exp_type,'_',num2str(s),'(AP-DAPI).jpg_Files\'];
        %             phr2=[phr1, num2str(bone_number(j)), exp_type,'_',num2str(s),'(trap).jpg_Files\'];
        %             phr=[phr1, num2str(bone_number(j)), exp_type,'_',num2str(s),'.jpg_Files\'];


                   % direc=direct;
%                     nn=bone_number{j};

%                     c=phr1;

%                     eval(['tt=threshold(',bone_number{j},',',num2str(s),',',num2str(loop),');'])
%                    if strcmp(exp_type,'V')
%                        eval(['file_name=''',direc,phr1, bone_number{k},'_h',exp_type,'_s',num2str(s),'_',num2str(loop-1),'_shift3.jpg'';']);
%                    elseif strcmp(exp_type,'F')
                        eval(['file_name=''',direc,phr1, 'L',num2str(s),'_s',bone_number{k},'_',num2str(loop-1),'_shift2.jpg'';']);
%                    end
%                     eval(['file_name=''',direc,phr1, '0' ,num2str(bone_number(j)),'_',num2str(s),'_',num2str(loop-1),'_shift2.jpg'';']);
                    a=dir(file_name);
                    if isempty(a)==0

                        eval(['im=imread(''',direc,a.name,''');'])
                        im=im(:,:,1);
%                         if loop~=4 & length(find(im>20))/numel(im)<0.1
%                             continue
%                         end


%                         eval(['bone',num2str(i),'(',num2str(k),',',num2str(s),',',num2str(loop),')=1;']);
%                         phr=['bone : ',num2str(i), ', bone number : ',bone_number{j},', section : ', num2str(s),', ',cell2mat(colors(loop))];
%                         if loop==1 | loop==2 | loop==4 | loop==6 | loop==9
%                             dic=0;
%                         elseif loop==3
%                             dic=1;
%                         elseif loop==7
%                             dic=2;
%                         end
                        if loop==1 | loop==9                % EGFP, GFP
                            dic=0;
                        elseif loop==3                      % DIC
                            dic=1;
                        elseif loop==4                      % TRAP
                            dic=3;
                        elseif loop==7                      % DAPI
                            dic=2;
                        elseif loop==2                      % AC
                            dic=4;
                        elseif loop==6                      % AP
                            dic=5;
                        end                        
                        [d, thre] = double_step_threshold_4(im,dic);
%                         if loop==7
%                         if loop==8
%                             [d, thre] = double_otsu(im);                        % Jan 30, 2014 for DAPIloop~=4 & loop ~=7                                   % TRAP is thresholded already from registration
%                         else
% %                             [d, thre] = double_step_threshold_3(im,dic);        % Jul 31, 2012
%                             [d, thre, BackgroundRemoved] = double_step_threshold_4(im,dic);        % Jul 21, 2014
%                         end
%                         th=graythresh(im)*255;
%                         if th>=40
%                             th=round(th*0.4);
%                         elseif th>=30 & th<40
%                             th=round(th*0.4);
%                         else
%                             th=round(th*0.5);
%                         end
% % %                         im1=im(:);
% % %                         im2=im1(im1>0);
% % %                         th=mean(im2);
% % % %                         d=im>th-20;
% % %                         d=im>th+std(double(im2));%/5;
%                         eval(['threshold_values',num2str(i),'(',num2str(bone_number(j)),',',num2str(s),',',num2str(loop),')=GUI_thresholding(im,phr);']);
%                        if strcmp(exp_type,'V')
%                            eval(['out_file_name=''',out_direc,phr1, bone_number{k},'_h',exp_type,'_s',num2str(s),'_',num2str(loop-1),'_shift3.jpg'';']);
%                        elseif strcmp(exp_type,'F')
                            eval(['out_file_name=''',out_direc,phr1, 'L',num2str(s),'_s',bone_number{k},'_',num2str(loop-1),'_shift3.jpg'';']);
%                        end
%                         eval(['out_file_name=''',out_direc,phr1,'0', num2str(bone_number(j)),'_',num2str(s),'_',num2str(loop-1),'_shift2.jpg'';']);
                        if loop==4                                          % TRAP is thresholded already from registration
                                                                            % TRAP-DAPI
                                                                            % July 26, 2013
%                             imwrite(im,out_file_name,'jpg','quality',100)
                            imwrite(d,out_file_name,'jpg','quality',100)
                        else
                            imwrite(d,out_file_name,'jpg','quality',100)
                        end
                        if loop==2 | loop==6                % AC and AP
                            eval(['out_file_name=''',direc,phr1, 'L',num2str(s),'_s',bone_number{k},'_',num2str(loop-1),'_shift3_BR.jpg'';']);
                            imwrite(BackgroundRemoved,out_file_name,'jpg','quality',100)
                        end
                    else
                        continue
                    end
                end
            end
        end
%         eval(['save ''',out_direc, 'threshold_values',num2str(i),''' threshold_values',num2str(i)'])
    end
    close
%     matlabpool close
end
 