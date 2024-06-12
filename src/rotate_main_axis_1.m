function rotate_main_axis_1(direct,gene_exp, exp_type)

% % direct='Z:\KOMP\Bzw_E01F\01_Submitted\Layers\Threshold\';
% direct='Z:\KOMP\Bzw_E01F\01_Submitted\Layers\Images\';
% % o_dir='Z:\KOMP\Bzw_E01F\01_Submitted\Layers\Threshold_Rotated\';
% gene_exp='Bzw2_E01_hF';
% exp_type='F';


if ~isempty(strfind(direct,'\'))
    delimeter='\';
elseif ~isempty(strfind(direct,'/'))
    delimeter='/';
end

% ed=strfind([direct, 'Images'],delimeter);
ed=strfind(direct,['Images',delimeter]);
o_dir=direct(1:ed-1)
% o_dir=[direct(1:ed-1),'Threshold'];

aa=dir(direct);
[rot_ang1]=init_rot_ang;
rot_ang2=rot_ang1;

% for i=1:length(aa)
%     if i==138
%         i;
%     end
%     n=aa(i).name;
              
'ImageToolbox License test start'
[a m]=license('checkout','image_Toolbox');
while (~a)
   pause(1);
   [a m]=license('checkout','image_Toolbox');
end
'ImageToolbox License test end'

for i=1:2
    switch i
        case 1
            gender='F';
        case 2
            gender='M';
    end
    for section=1:3
        for bn=1:8
%             n=[direct,gene_exp,'_',gender,'L',num2str(section),'_s',num2str(bn),'_6_shift2.jpg'];
            n=[direct,gene_exp,'_',gender,'L',num2str(section),'_s',num2str(bn),'_shift2.jpg'];
            if isempty(dir(n))
                n=[direct,gene_exp,'_',gender,'L',num2str(section),'_s',num2str(bn),'_2_shift2.jpg'];
                if isempty(dir(n))   
                    continue
                end
            end
            n
            if exp_type=='F'
                a=imread(n);
%                 a=a(:,:,1);
%                 a=a>127;
                a=sum(a,3);
                a_1=a>graythresh(a)*max(a(:))*0.3;
%                 temp=sum(a,1);
%                 ind=find(temp==max(temp));
%                 ind_temp=abs(ind-size(a,2));
%                 ind=ind_temp(find(ind_temp==min(ind_temp)));
%                 ind=ind(1);
%                 a_1=a;
%                 a_1(:,ind)=true(size(a,1),1);
                [ang] = find_rot_angle_of_main_axis_using_DAPI(a_1)
            elseif exp_type=='V'
                ang=90;
            end

            switch gender
                case 'F'
                    rot_ang1{bn,section+1}=ang;
                case 'M'
                   rot_ang2{bn,section+1}=ang;
            end
        end
    end
end
% eval(['save ''',direct,'rot_ang'' rot_ang*'])
eval(['save ''',o_dir,'rot_ang'' rot_ang*'])

function [rot_ang]=init_rot_ang
for i=1:8
    rot_ang{i,1}=num2str(i);
end
rot_ang{8,4}=[];
