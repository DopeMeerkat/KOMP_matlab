function [b] = find_dips_and_remove_1(b, a, ap)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   find left dip
%
%   October 05, 2011
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(find(ap,1,'first'))
    return
end
t=sum(b,1);
ind=find(t,1,'first');
t=b(:,ind:ind+1000);
[x1 y1]=find(t);
for h=1:max(y1)
    profile(h)=min(x1(y1==h));
end
profile(end)=0;
d=diff(profile); clear profile

d1=find(d~=0);
dd=diff(d1);
ddd=diff(dd);
id=find(ddd==max(ddd));
st=d1(id);
ed=dd(id+1)+st;

% st=max(find(d>0));
% % ed=st+min(find(d(st:end)<0));
% d1=d(st:end);
% ed=min(find(d1==(min(d1(d1<0)))))+st-1;

if isempty(st)==0 & isempty(ed)==0
    st=st+ind;
    ed=ed+ind;
    y_val=min(find(b(:,st)));
    s=ed-st+1;
    tt=a(y_val:y_val+2*s-1,st:ed);
    for x=1:size(tt,2)
       test=tt(:,x);
       p=min(find(test));
       if isempty(p)==0
           tt(p:end,x)=ones(size(tt,1)-p+1,1);
       end
    end
    tt=sum(tt,2);
%     before=2*s;
%     for y=y_val:y_val+2*s-1
%         number=length(find(a(y,st:ed)));
%         if number<=before
%             b(y,st:ed)=zeros(1,s);
%             before=number;
%         else
%             break
%         end
%     end
    ttt=diff(tt);
    tttt=find(ttt~=0);
    ttttt=diff(tttt);
    stop_point=tttt(find(ttttt==max(ttttt)))+ttt(find(ttttt==max(ttttt)));
%     pos=find(ttt>=0);
%     neg=find(ttt<0);
%     for i=2:length(neg)
%         test=neg(i)-min(pos(pos<neg(i)&pos>neg(i-1)));
%         if isempty(test)==0
%             temp(i)=test;
%         end
%     end
%     stop_point=neg(find(temp==max(temp)));


%     b(y_val:y_val+stop_point,st:ed)=zeros(stop_point+1,s);
%%%%%%%%%%%%%%%
%
%   added on November 10, 2011
%

    for i=st:ed
%         st_a(i)=find(a(:,i),1,'first');
        temp=find(a(:,i),1,'first');
        if isempty(temp)
            st_a(i)=NaN;
            continue
        else
            st_a(i)=temp;
        end
%         st_b(i)=find(b(:,i),1,'first');
        temp=find(b(:,i),1,'first');
        if isempty(temp)
            st_b(i)=NaN;
            continue
        else
            st_b(i)=temp;
        end
        temp_ap=ap(st_a(i)-100:st_a(i)+1000,i);
        temp=find(temp_ap,1,'first');
        if isempty(temp)
            point_ap(i)=10000;
            continue
        end
        point_ap(i)=st_a(i)+temp-100;
    end
    minimum=min(point_ap(st:ed), st_a(st:ed));
    for i=st:ed
%         b(st_b(i):minimum(i-st+1),i)=zeros(minimum(i-st+1)-st_b(i)+1,1);
        if isnan(st_b(i)) | isnan(minimum(i-st+1)) | st_b(i)==0 | minimum(i-st+1)==0
            continue
        else
            b(st_b(i):minimum(i-st+1),i)=zeros(minimum(i-st+1)-st_b(i)+1,1);
        end
    end
%     for i=st:ed
%         temp_ap=ap(y_val+stop_point+1:y_val+stop_point+1+1000,i);
%         temp=find(temp_ap,1,'first');
%         if isempty(temp)
%             continue
%         end
%         point_ap=y_val+stop_point+temp;
%         b(y_val+stop_point+1:point_ap,i)=zeros(point_ap-y_val-stop_point,1);
%     end
%
%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   find right dip
%
%   October 05, 2011
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t=sum(b,1);
ind=find(t,1,'last');
t=b(:,ind-1000:ind);
[x1 y1]=find(t);
for h=1:max(y1)
    profile(h)=min(x1(y1==h));
end
profile(1)=0;
d=diff(profile); clear profile

d1=find(d~=0);
dd=diff(d1);
ddd=diff(dd);
id=find(ddd==max(ddd));
st=d1(id);
ed=dd(id+1)+st;

% % ed=min(find(d<0));
% % st=max(find(d(1:ed)>0));
% ed=min(find(d==(min(d(d<0)))));
% d1=d(1:ed);
% st=max(find(d1==(max(d1(d1>0)))));
% if isempty(st)
%     st=find(d==max(d))-1;
% %     ed=min(find(d(st:end)<0))+st-1;
%     d1=d(st:end);
%     ed=min(find(d1==(min(d1(d1<0)))))+st-1;
% end
ind1=ind-1000;
st=st+ind1;
ed=ed+ind1;
if isempty(st)==0 & isempty(ed)==0
    y_val=min(find(b(:,st)));
    s=ed-st+1;
%     before=2*s;
%     for y=y_val:y_val+2*s-1
%         number=length(find(a(y,st:ed)));
%         if number<=before
%             b(y,st:ed)=zeros(1,s);
%             before=number;
%         else
%             break
%         end
%     end
% %     b(y_val:y_val+2*s-1,st:ed)=zeros(2*s,s);
%     before=2*s;
%     for y=y_val:y_val+2*s-1
%         number=length(find(a(y,st:ed)));
%         if number<=before
%             b(y,st:ed)=zeros(1,s);
%             before=number;
%         else
%             break
%         end
%     end
    tt=a(y_val:y_val+2*s-1,st:ed);
    for x=1:size(tt,2)
       test=tt(:,x);
       p=min(find(test));
       if isempty(p)==0
           tt(p:end,x)=ones(size(tt,1)-p+1,1);
       end
    end
    tt=sum(tt,2);
    ttt=diff(tt);
    tttt=find(ttt~=0);
    ttttt=diff(tttt);
    stop_point=tttt(find(ttttt==max(ttttt)))+ttt(find(ttttt==max(ttttt)));
%     pos=find(ttt>=0);
%     neg=find(ttt<0);
%     for i=2:length(neg)
%         test=neg(i)-min(pos(pos<neg(i)&pos>neg(i-1)));
%         if isempty(test)==0
%             temp(i)=test;
%         end
%     end
%     stop_point=neg(find(temp==max(temp)));



%     b(y_val:y_val+stop_point,st:ed)=zeros(stop_point+1,s);
%%%%%%%%%%%%%%%
%
%   added on November 10, 2011
%
    clear st_a st_b minimum
    for i=st:ed
        temp=find(a(:,i),1,'first');
        if isempty(temp)
            st_a(i)=NaN;
            continue
        else
            st_a(i)=temp;
        end
        temp=find(b(:,i),1,'first');
        if isempty(temp)
            st_b(i)=NaN;
            continue
        else
            st_b(i)=temp;
        end
        temp_ap=ap(max(1,st_a(i)-100):min(size(ap,1),st_a(i)+1000),i);
        temp=find(temp_ap,1,'first');
        if isempty(temp)
            point_ap(i)=10000;
            continue
        end
        point_ap(i)=st_a(i)+temp-100;
    end
    minimum=min(point_ap(st:ed), st_a(st:ed));
    for i=st:ed
%         b(st_b(i):minimum(i-st+1),i)=zeros(minimum(i-st+1)-st_b(i)+1,1);
        if isnan(st_b(i)) | isnan(minimum(i-st+1))
            continue
        else
            b(st_b(i):minimum(i-st+1),i)=zeros(minimum(i-st+1)-st_b(i)+1,1);
        end
    end
%     for i=st:ed
%         temp_ap=ap(y_val+stop_point+1:y_val+stop_point+1+1000,i);
%         temp=find(temp_ap,1,'first');
%         if isempty(temp)
%             continue
%         end
%         point_ap=y_val+stop_point+temp;
%         b(y_val+stop_point+1:point_ap,i)=zeros(point_ap-y_val-stop_point,1);
%     end
%
%%%%%%%%%%%%%%%
end

[L n]=bwlabel(b|a);
stat=regionprops(L, 'FilledImage','BoundingBox', 'Area');
ar=[stat.Area];
ind=find(ar==max(ar));
b=zeros(size(b));
st_x=round(stat(ind).BoundingBox(2));
st_y=round(stat(ind).BoundingBox(1));
ed_x=st_x+stat(ind).BoundingBox(4)-1;
ed_y=st_y+stat(ind).BoundingBox(3)-1;
b(st_x:ed_x,st_y:ed_y)=stat(ind).FilledImage;

return