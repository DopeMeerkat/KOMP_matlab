function [dic, dic2]=find_periosteum_4(dic)

dic(1,:)=false(1,size(dic,2));
dic(end,:)=false(1,size(dic,2));

dic1=imclose(dic,strel('disk',70));
[L n]=bwlabel(dic1); clear dic1
stat=regionprops(L, 'FilledImage','BoundingBox','Area'); clear L
temp=false(size(dic));
for i=1:n
    sx=round(stat(i).BoundingBox(2));
    sy=round(stat(i).BoundingBox(1));
    ex=sx+stat(i).BoundingBox(4)-1;
    ey=sy+stat(i).BoundingBox(3)-1;
    temp(sx:ex,sy:ey)=temp(sx:ex,sy:ey)|stat(i).FilledImage;
end
[L n]=bwlabel(temp);
stat=regionprops(L, 'FilledImage','BoundingBox','Area'); clear L
a=[stat.Area];
ind=find(a==max(a(:)));
dic2=false(size(dic));
sx=round(stat(ind).BoundingBox(2));
sy=round(stat(ind).BoundingBox(1));
ex=sx+stat(ind).BoundingBox(4)-1;
ey=sy+stat(ind).BoundingBox(3)-1;
dic2(sx:ex,sy:ey)=dic2(sx:ex,sy:ey)|stat(ind).FilledImage;

dic=dic2&dic;
return