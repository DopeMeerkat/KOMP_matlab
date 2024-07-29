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
% info1

%     if str2num(channel)==1
%         bn=[bn, str2num(bone_number)];
%         row=row+1;
%         if length(find(bn==str2num(bone_number)))>1
%             bn=str2num(bone_number);
%             row=1;
%         end
%     end
%     


%     for ii=1:length(bn)
%         if bone_number==bn(ii)
%             row=ii;
%             flag=0;
%             break
%         end
%     end
%     if flag==1
%         row=row+1;
%     end
    
    
    
%     if gender=='M' 
%         if first
%             row=1;
%             info1=info;
%             info=info_M;
%             first=0;
%         end
%     elseif gender=='F'
% 
%     end
%     if filter_type=='M' 
%         if first_filter_m
%             row=1;
%             first_filter_m=0;
%         end
%     elseif filter_type=='T' 
%         if first_filter_t
%             row=1;
%             first_filter_t=0;
%         end
%     end
%     info.im{row,1}=bone_number;
%     info.ap{row,1}=bone_number;
%     info.tr{row,1}=bone_number;


%     im{row,1}=bone_number;
%     ap{row,1}=bone_number;
%     tr{row,1}=bone_number;
%     
%     switch colors
%         case 'A1'
%             stain='DAPI';
%         case 'A2'
%             stain='G_beads_A';
%         case 'A3'
%             stain='AP';
% %             info.ap{row,str2num(section_number)+1}=1;
%             ap{row,str2num(section_number)+1}=1;
%         case 'T1'
%             stain='G_beads_T';
%         case 'T2'
%             stain='TRAP';
% %             info.tr{row,str2num(section_number)+1}=1;
%             tr{row,str2num(section_number)+1}=1;
%         case 'T3'
%             stain='R_beads_T';
% %         case 'M1'
% %             stain='RAC';
% %         case 'M2'
% %             stain='Mineral';
% %             info.im{row,str2num(section_number)+1}=1;
% %         case 'M3'
% %             stain='G_beads_M';
% %         case 'M4'
% %             stain='R_beads_M';
%         case 'M1'
%             stain='Mineral';
% %             info.im{row,str2num(section_number)+1}=1;
%             im{row,str2num(section_number)+1}=1;
%         case 'M2'
%             stain='G_beads_M';
%         case 'M3'
%             stain='R_beads_M';
%     end
%     switch filter_type
%         case 'A'
%             info.ap=ap;
%         case 'M'
%             info.im=im;
%         case 'T'
%             info.tr=tr;
%     end
%     switch gender
%         case 'F'
%             info1=info;
%         case 'M'
%             info2=info;
%     end
% end
% info2=info; clear info

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


