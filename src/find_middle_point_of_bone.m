function [third]=find_middle_point_of_bone(a)
[L n]=bwlabel(a);
stat=regionprops(L, 'BoundingBox','Image','Area'); clear L
ar=[stat.Area];
test=ar/sum(ar);
ind=find(test>0.1);
test=false(size(a));
for i=ind
    st_x=round(stat(i).BoundingBox(2));
    st_y=round(stat(i).BoundingBox(1));
    ed_x=st_x+stat(i).BoundingBox(4)-1;
    ed_y=st_y+stat(i).BoundingBox(3)-1;
    test(st_x:ed_x,st_y:ed_y)=test(st_x:ed_x,st_y:ed_y)|stat(i).Image;
end
clear stat
count=1;
for i=1:size(test,1)
    temp=test(i,:);
    x=find(temp);
    if ~isempty(x)
        left(count)=min(x);
        right(count)=max(x);
        count=count+1;
    end
end
third=round((median(left)+median(right))/2);
return