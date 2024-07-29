function [new_points]=find_points_after_rotate(ori_size, size_im, points, ang)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Find new growth_plate poiint and end point after rotate
%       to make the image straight up
%   
%   size_im     : input image size before rotate (vertical size only)
%   points      : growth plate point, end point point (vertical axis only)
%   ang         : angle to rotate (in degree)
%
%   new_points	: growth plate point, end point point after rotation
%               (vertical axis only)
%   May 6, 2014
%   Sean Hong
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ang==0
    new_points=points;
    return
end
for i=1:length(points)
    middle=round(ori_size/2);
    ang_rad=ang*pi/180/2;
    if points(i)<middle
        len=middle-points(i);
%         new_points(i)=round(size_im/2)-round(len*sqrt(1/(1+(1/tan(pi/2-ang_rad))^2)));
%         new_points(i)=len*(1-tan(ang_rad)^2)/(1+tan(ang_rad)^2);
        new_points(i)=(round(size_im/2)-middle)+points(i)+(len-len*(1-tan(ang_rad)^2)/(1+tan(ang_rad)^2));
    elseif points(i)>middle
        len=points(i)-middle;
%         new_points(i)=round(size_im/2)+round(len*sqrt(1/(1+(1/tan(pi/2-ang_rad))^2)));
%         new_points(i)=size_im-len*(1-tan(ang_rad)^2)/(1+tan(ang_rad)^2);
        new_points(i)=round(size_im/2)+len*(1-tan(ang_rad)^2)/(1+tan(ang_rad)^2);
    end
end
