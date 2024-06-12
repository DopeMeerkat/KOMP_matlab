function [cortical] = find_cortical(test_test_im, result, roi)

[x y]=find(roi);
top=min(x);
bottom=max(x);

template=test_test_im-imerode(test_test_im,strel('disk',1));
template(1:top-1,:)=false(top-1,size(test_test_im,2));
template(bottom+1:end,:)=false(size(test_test_im,1)-bottom,size(test_test_im,2));

% test_test_im(1:top-1,:)=false(top-1,size(test_test_im,2));
% test_test_im(bottom+1:end,:)=false(size(test_test_im,1)-bottom,size(test_test_im,2));

[L n]=bwlabel(template);
stat=regionprops(L, 'BoundingBox', 'Image', 'Area'); clear L
% ar=[stat.Area];
% sar=sort(ar);
% ind1=find(ar==sar(end));
% ind2=find(ar==sar(end-1));
for i=1:n
    len(i)=stat(i).BoundingBox(4);
end
slen=sort(len);
ind=find(len==slen(end));
ind1=ind(1);
ind2=ind(2);
if stat(ind1).BoundingBox(1)>stat(ind2).BoundingBox(1)
    left=ind2;
    right=ind1;
else
    left=ind1;
    right=ind2;
end

cortical=false(size(roi));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   left cortex
temp=false(size(roi));
st_x=round(stat(left).BoundingBox(2));
ed_x=st_x+stat(left).BoundingBox(4)-1 ;
st_y=round(stat(left).BoundingBox(1));
ed_y=st_y+stat(left).BoundingBox(3)-1 ;
temp(st_x:ed_x,st_y:ed_y)=stat(left).Image;
temp_cortical=temp;

cont=1;
count=0;
s(1)=length(find(temp&result));
while cont
    count=count+1;
    s(count+1)=length(find(temp(:,1:end-count+1)&result(:,count:end)));
    if s(count+1)/s(1)<0.75 | count>500
        break
    end
    temp_cortical(:,count:end)=temp_cortical(:,count:end) | temp(:,1:end-count+1); 
end
cortical=cortical | (temp_cortical & result);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   right cortex
temp=false(size(roi));
st_x=round(stat(right).BoundingBox(2));
ed_x=st_x+stat(right).BoundingBox(4)-1 ;
st_y=round(stat(right).BoundingBox(1));
ed_y=st_y+stat(right).BoundingBox(3)-1 ;
temp(st_x:ed_x,st_y:ed_y)=stat(right).Image;
temp_cortical=temp;
clear s

cont=1;
count=0;
s(1)=length(find(temp&result));
while cont
    count=count+1;
    s(count+1)=length(find(temp(:,count:end)&result(:,1:end-count+1)));
    if s(count+1)/s(1)<0.75 | count>500
        break
    end
    temp_cortical(:,1:end-count+1)=temp_cortical(:,1:end-count+1) | temp(:,count:end); 
end
cortical=cortical | (temp_cortical & result);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%