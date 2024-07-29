function [ang] = find_rot_angle_of_main_axis_using_DAPI(a)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   a   : Thresholded DAPI
%
%   ang : rotated angle
%
%   March 12, 2014
%   Sean Hong
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if length(find(a))/numel(a)<0.05
    a1=imclose(a,strel('disk',30));
else
    a1=imclose(a,strel('disk',20));
end
[L n]=bwlabel(a1);
stat=regionprops(L, 'BoundingBox','Image','Area');clear L
ar=[stat.Area];
sar=sort(ar);
ind=(ar==sar(end));
[L n]=bwlabel(stat(ind).Image);
stat=regionprops(L, 'BoundingBox','Image','Area','Orientation');clear L
ang=stat.Orientation;
% if ang<0
%     a1=imrotate(a,-(90+ang));
% elseif ang>0
%     a1=imrotate(a,(90-ang));
% end