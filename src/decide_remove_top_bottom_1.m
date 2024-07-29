function [remove_top, remove_bottom, min_y, max_y]=decide_remove_top_bottom_1(a, min_y, max_y)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   input   : a (uint8 format)
%           : min_y --> vertical start point (remove above growth plate)
%           : max_y --> vertical end point (remove below growth plate)
%   output  : remove_top (1 --> remove, 0 --> do not remove)
%             remove_bottom (1 --> remove, 0 --> do not remove)
%   January 11, 2011
%
%   Sean Hong
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a(min_y-1,:)=zeros(1,size(a,2));
a(max_y+1,:)=zeros(1,size(a,2));

[L n]=bwlabel(a);
stat=regionprops(L,'Image','BoundingBox','Area');clear L
a=zeros(size(a));
for i=1:n
    sx=round(stat(i).BoundingBox(2));
    sy=round(stat(i).BoundingBox(1));
    ex=sx+stat(i).BoundingBox(4)-1;
    ey=sy+stat(i).BoundingBox(3)-1;
    if (sx<min_y & max_y-sx>6000) | (ex>max_y & ex-min_y>6000)
        continue
    else
        a(sx:ex,sy:ey)=a(sx:ex,sy:ey)+stat(i).Image;
    end
end

p=sum(a,2);
if find(p,1,'last')<=max_y & find(p,1,'first')>=min_y
    remove_top=1;
    remove_bottom=1;
    return
end
p=medfilt2(p,[51,1]);

d=mean(p)/3;

pp=diff(p);

if min_y==1
    val=d-1;
else
    pp1=pp(1:min_y);
    count=1;

    if pp1(1)>0
        init=1;
    elseif pp1(1)<0
        init=-1;
    else
        init=0;
    end
    before=init;
    val(count)=pp1(1);
    for i=2:length(pp1)
        if pp1(i)~=0
            sign=pp1(i)/abs(pp1(i));
            if sign==before
                val(count)=val(count)+pp1(i);
            else
                count=count+1;
                val(count)=pp1(i);
                cordinate(count)=i;
            end
            before=sign;
        end
    end
end

new_ind=find((val>mean(val(find(val>0)))),1,'first');
if isempty(new_ind)==0
    min_y=cordinate(new_ind);
end

if max(abs(val))>d & isempty(new_ind)==0
    remove_top=1;
else
    remove_top=0;    
end

clear cordinate
if max_y>=length(pp)
    val=d-1;
else
    pp1=pp(max_y:end);
    count=1;
    clear val

    if pp1(1)>0
        init=1;
    elseif pp1(1)<0
        init=-1;
    else
        init=0;
    end
    before=init;
    val(count)=pp1(1);
    for i=2:length(pp1)
        if pp1(i)~=0
            sign=pp1(i)/abs(pp1(i));
            if sign==before
                val(count)=val(count)+pp1(i);
            else
                count=count+1;
                val(count)=pp1(i);
                cordinate(count)=i;
            end
            before=sign;
        end
    end
end

new_ind=find((val<mean(val(find(val<0)))),1,'last');
if isempty(new_ind)==0
    max_y=max_y+cordinate(new_ind);
end

if max(abs(val))>d & isempty(new_ind)==0
    remove_bottom=1;
else
    remove_bottom=0;
end

