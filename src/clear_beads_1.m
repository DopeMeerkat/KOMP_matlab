function [r1, g1, ap1]=clear_beads_1(file_name, d1,r1,g1,ap1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Remove beads (circular shaped objects) in red, green, AP channels
%
%   input   : file_name including path (cell structure) 
%               file_name={'C:\Users\shhong\Desktop\His2\His2_FL_FT & FC\01_Submitted\Images\His2_FC_1_1_2_shift2.jpg';
%                           'C:\Users\shhong\Desktop\His2\His2_FL_FT & FC\01_Submitted\Images\His2_FC_1_1_1_shift2.jpg';
%                           'C:\Users\shhong\Desktop\His2\His2_FL_FT & FC\01_Submitted\Images\His2_FC_1_1_0_shift2.jpg';
%                           'C:\Users\shhong\Desktop\His2\His2_FL_FT & FC\01_Submitted\Images\His2_FC_1_1_5_shift2.jpg'}
%             threshold (threshold values of each image)
%             mouse : 1 --> mouse (default)
%                     0 --> rat
%
%   output  : r1 (red channel without beads in uint8)               
%             g1 (green channel without beads in uint8)
%             ap1 (AP channel without beads in uint8)
%
%   Author  : Sean Hong
%   Date    : March 30, 2011
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if mouse
    min_size=50;
    max_size=500;
% else
%     min_size=10;
%     max_size=190;
% end

% switch length(file_name)
%     case 1
%         d=imread(file_name{1,:}, 'jpg');
%         d=d(:,:,1);
%         r=uint8(zeros(size(d)));
%         g=uint8(zeros(size(d)));
%         ap=uint8(zeros(size(d)));
%         threshold(2:4)=[255, 255, 255];
%     case 2
%         d=imread(file_name{1,:}, 'jpg');
%         d=d(:,:,1);
%         r=imread(file_name{2,:}, 'jpg');
%         g=uint8(zeros(size(d)));
%         ap=uint8(zeros(size(d)));
%         threshold(3:4)=[255, 255];
%     case 3
%         d=imread(file_name{1,:}, 'jpg');
%         d=d(:,:,1);
%         r=imread(file_name{2,:}, 'jpg');
%         g=imread(file_name{3,:}, 'jpg');
%         ap=uint8(zeros(size(d)));
%         threshold(4)=[255];
%     case 4
%         d=imread(file_name{1,:}, 'jpg');
%         d=d(:,:,1);
%         r=imread(file_name{2,:}, 'jpg');
%         g=imread(file_name{3,:}, 'jpg');
%         ap=imread(file_name{4,:}, 'jpg');
% end
% 
% ap=sum(ap,3); ap=uint8(ap/max(ap(:))*255);
% d=sum(d,3); d=uint8(d/max(d(:))*255);
% r=sum(r,3); r=uint8(r/max(r(:))*255);
% g=sum(g,3); g=uint8(g/max(g(:))*255);

temp_name=file_name{1,:};
underscore=strfind(temp_name,'_');
o_file_no_bead=[temp_name(1:underscore(end-1)),temp_name(underscore(end)+1:end-4),'_NoBeads.jpg'];
o_file_bead=[temp_name(1:underscore(end-1)),temp_name(underscore(end)+1:end-4),'_Beads.jpg'];

% d1=(d>threshold(1));
% g1=(g>threshold(3));
% r1=(r>threshold(2));
% ap1=(ap>threshold(4));
if nargin<5
    ap1=uint8(false(size(d1)));
end
a=imadd(uint8(d1)*56,imadd(uint8(r1)*255, uint8(ap1)*255));
a(:,:,2)=imadd(uint8(d1)*56,imadd(uint8(g1)*255, uint8(ap1)*128));
a(:,:,3)=uint8(d1)*56;
imwrite(a, o_file_bead, 'jpg', 'quality', 100)

[g2]=find_circles(g1,min_size,max_size);
g1=xor(g1,g2) & (g1 & ~g2);
[r2]=find_circles(r1,min_size,max_size);
r1=xor(r1,r2) & (r1 & ~r2);
if nargin==5
    [ap2]=find_circles(ap1,min_size,max_size);
    ap1=xor(ap1,ap2) & (ap1 & ~ap2);
end
a=imadd(uint8(d1)*56,imadd(uint8(r1)*255, uint8(ap1)*255));
a(:,:,2)=imadd(uint8(d1)*56,imadd(uint8(g1)*255, uint8(ap1)*128));
a(:,:,3)=uint8(d1)*56;
%figure;imshow(a)
imwrite(a, o_file_no_bead, 'jpg', 'quality', 100)

return


function [g1]=find_circles(g, min_size, max_size)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
%
%       input : g (input with beads)
%               min_size( minimum size of beads)
%               max_size( maximum size of beads)
%
%       output : g1 (image with found beads)
%
%       date  : January 6, 2010
%       author : Sean Hong
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55

[L n]=bwlabel(g);
stats=regionprops(L,'BoundingBox', 'Image', 'Area'); clear L
a=[stats.Area];
x=find(a>min_size & a<max_size);
g1=false(size(g));
for i=1:length(x)
    sx=round(stats(x(i)).BoundingBox(2));
    sy=round(stats(x(i)).BoundingBox(1));
    ex=sx+stats(x(i)).BoundingBox(4)-1;
    ey=sy+stats(x(i)).BoundingBox(3)-1;
    test=a(x(i))/(((max(stats(x(i)).BoundingBox(4), stats(x(i)).BoundingBox(3)))/2)^2*pi);
    radius_ratio=stats(x(i)).BoundingBox(4)/stats(x(i)).BoundingBox(3);
%    if (test>0.75 & test <1.25) & (radius_ratio>0.75 & radius_ratio <1.25)
    if (test>0.8 & test <1.2) & (radius_ratio>0.85 & radius_ratio <1.15)
        g1(sx:ex,sy:ey)=stats(x(i)).Image;
    end
end
return
