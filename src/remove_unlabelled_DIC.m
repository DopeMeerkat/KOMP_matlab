function [out]=remove_unlabelled_DIC(dic,rg)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   October 24, 2011
%
%   input
%       dic     : dic image (logical)
%       rg      : label image (logical) --> red and green
%   output
%       out     : unlebelled DIC removed DIC image (logical)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

size_dilation=20;
[L n]=bwlabel(dic);
stat=regionprops(L, 'Image', 'BoundingBox', 'Area'); clear L
out=false(size(dic));

for i=1:n
    start_x=round(stat(i).BoundingBox(2));
    start_y=round(stat(i).BoundingBox(1));
    end_x=start_x+stat(i).BoundingBox(4)-1;
    end_y=start_y+stat(i).BoundingBox(3)-1;
    temp=rg(start_x:end_x,start_y:end_y) & stat(i).Image;
    if isempty(find(temp(:),1))==1
        continue
    else
        temp1=false(size(temp)+size_dilation*2);
        temp1(size_dilation+1:end-size_dilation,size_dilation+1:end-size_dilation)=stat(i).Image;
        temp1=imclose(temp1, strel('disk',size_dilation));

        out(start_x:end_x,start_y:end_y)=out(start_x:end_x,start_y:end_y)...
            | temp1(size_dilation+1:end-size_dilation,size_dilation+1:end-size_dilation);
    end
end