function [d, thre]=double_otsu(a)

temp=a>graythresh(a)*max(a(:));
if length(find(temp))/numel(temp)<0.10
    d=temp;
    thre=graythresh(a)*max(a(:));
    return
end
temp1=immultiply(a,uint8(temp));
t=sort(temp1(:));
t1=(t(find(t,1,'first'):end));

% figure;imshow(a>graythresh(t1)*max(t(:)))
thre=graythresh(t1)*max(t(:));
d=uint8(a>thre)*255;