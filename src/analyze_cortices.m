function [average_height, average_width, m_cortical, n, osteocytes]=analyze_cortices(tr, cortical, DAPI, dis_10x)
[x y]=find(tr);
average_height=max(x)-min(x);
width=zeros(1,average_height);
for i=min(x):max(x)
    temp=tr(i,:);
    width(i)=length(find(temp));
end
average_width=mean(width);

m_cortical=cortical & ~tr(dis_10x:end,:);
temp1=m_cortical(1:max(x)-dis_10x,:);
% temp=imdilate(m_cortical,strel('disk',20));
% m_cortical_frame=temp-imerode(temp,strel('disk',3));

% osteocytes=DAPI(start_point+dis_10x:start_point+size(m_cortical,1)-1+dis_10x,:) & m_cortical;
osteocytes=DAPI(dis_10x:size(temp1,1)-1+dis_10x,:) & temp1;
[L n]=bwlabel(osteocytes);
stat=regionprops(L, 'Area','Image','BoundingBox');clear L
n=length(find([stat.Area]<100));
% m_cortical=false(size(temp1));
% for i=1:n
%     if stat(i).Area<100                     % too big is not DAPI
%         st_x=round(stat(i).BoundingBox(2));
%         st_y=round(stat(i).BoundingBox(1));
%         ed_x=st_x+stat(i).BoundingBox(4)-1;
%         ed_y=st_y+stat(i).BoundingBox(3)-1;
%         m_cortical(st_x:ed_x,st_y:ed_y)=m_cortical(st_x:ed_x,st_y:ed_y)|stat(i).Image;
%     end
% end

% number_osteocytes_per_area=n/(length(find(m_cortical))*(400/dis_10x)^2);
% osteocytes_area_ratio=length(find(osteocytes))/length(find(m_cortical));
return