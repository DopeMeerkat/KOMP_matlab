function [d, thre] = double_step_threshold_4_batch_CCCE09(a,dic)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   input   : a (uint8 format)
%             dic : 0 --> EGFP, GFP
%                   1 --> DIC
%                   2 --> DAPI
%                   3 --> TRAP
%                   4 --> AC
%                   5 --> AP
%   output  : d (logical format)
%
%   October 27, 2011
%
%   Sean Hong
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   OCtober 13, 2011
%
% thre=graythresh(a)*255;                   % removed on October 13, 2011

% ang=-10;
% if dic==1 | dic==4 | dic==5         % DIC, AC, AP
%     ang=-10;
% else
%     ang=-1;
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   July 21, 2014
if dic==1                           % DIC
%     ang=-10;
    ang=-10;
elseif dic==0                       % EGFP, GFP
    ang=-10;
elseif dic==3                       % TRAP
    ang=-1;
elseif dic==2                       % DAPI
%     ang=-5;
    ang=-50;
elseif dic==4                       % AC
%     ang=-60;
    ang=-85;
elseif dic==5                       % AP
%     ang=-60;
    ang=-85;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

thresh_offset=0;
low=0;

cont=1;
count=1;
% while cont
%     if dic==1                       % DIC
%         dividend=5;
%     elseif dic==0                   % labels, TRAP, AP
%         dividend=100;
%     elseif dic==2                   % DAPI
%         dividend=7;
%     end

    switch dic                          % if not DIC, use medfilt2 for background removal
        case 0                          % EGFP, GFP
%             w_size=31;
            w_size=51;
        case 1                          % DIC
            w_size=1201;
        case 2                          % DAPI
            w_size=21;
        case 3                          % TRAP
            w_size=81;
        case 4                          % AC
%             w_size=31;            
            w_size=51;            
        case 5                          % AP
%             w_size=51;
            w_size=71;
    end

    if dic==0 | dic==3 | dic==4 | dic==5 | dic==2
        a12=medfilt2_resize(a, w_size);   % image resize scale is 0.5
        a=a-a12;
    end
    
    h=hist(double(a(:)),255);
    if isempty(find(h,1))
        d=a; thre=0;
        return
    end
    h(1)=0;
    h(255)=0;
%     h=medfilt2(h, [1,7]);        % November 7, 2011
    h=medfilt2(h, [1,4]);        % November 7, 2011

    h=h/max(h)*255;
    h;
    len=10;
    for i=1:255-len
        cur(i)=atan((h(i+len)-h(i))/len)*180/pi;
    end
    ind=find(h==max(h(:)),1,'last');

%     cur=medfilt2(cur,[1,7]);        % November 7, 2011
    cur=medfilt2(cur,[1,4]);        % November 7, 2011
    
    thre=ind+find(cur(ind:end)>ang,1,'first')+len/2+thresh_offset;
    if isempty(thre)
        d=a>5; thre=0;
        return
    end

    h1=h(thre:end);
    ind1=find(h1==max(h1(:)),1,'last');
    if h(ind1+thre)>h(ind)/2 & ind1<50
        if isempty(find(cur(ind1:end)>ang,1,'first'))==0
            thre=thre+ind1+find(cur(thre+ind1:end)>ang,1,'first')+len/2+thresh_offset;
        end
    end

    % thre=min(find(h(ind+1:end)<h(ind)/dividend))+ind;

%
%   OCtober 13, 2011
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if dic==2
        thre=graythresh(a)*max(a(:));
    end

    d=(a>thre);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Jul 31, 2012
%
    if low==1
        cont=0;
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%     if (dic==1 | dic==0 | dic==3) & length(find(d(:)))/numel(d)>0.3

%     if (dic==1 | dic==0 | dic==2 | dic==3) & length(find(d(:)))/numel(d)>0.25 %0.3

%     if (dic==1 | dic==2) & length(find(d(:)))/numel(d)>0.25 %0.3
    if (dic==1 & length(find(d(:)))/numel(d)>0.25)  | (dic==2 & length(find(d(:)))/numel(d)>0.5)
        if count>1
            cont=0;
%             continue
        end
        count=count+1;
%         a=a-medfilt2(a,[1001,1001]);

        ad=double(a);
        ad=ad/max(ad(:));
        a=ad.^.7;        
        a12=medfilt2_resize(a, w_size);   % image resize scale is 0.5
        a=a-a12;
        switch dic
            case 1
                ['Medfilt2 is applied with resize : too strong DIC']
            case 0
                ['Medfilt2 is applied with resize : too strong EGFP']
            case 2
                ['Medfilt2 is applied with resize : too strong DAPI']
            case 3
                ['Medfilt2 is applied with resize : too strong TRAP']
        end

%     elseif dic==3 & length(find(d(:)))/numel(d)<0.1

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Jul 31, 2012
%
    elseif dic==1 & length(find(d(:)))/numel(d)<0.15 
        if count>1
            cont=0;
%             continue
        end
        count=count+1;
%         a=a-medfilt2(a,[1001,1001]);

        ad=double(a);
        ad=ad/max(ad(:));
        a=ad.^.4;
        a1=uint8(a*255);
%         a12=medfilt2_resize(a, w_size);   % image resize scale is 0.5
        a12=medfilt2_resize(a1, w_size);
        a=a1-a12;
        st=1.5;
        ed=-0.5;
        step=0.1;
        count=round((st-ed)/step)+1;
        for c=st:-step:ed
            thre=graythresh(a(:))*max(a(:))-std(double(a(:)))*c;
            d=a>thre;
            [L,n]=bwlabel(d);
            stats=regionprops(L,'Area');
            ar=[stats.Area];

            dots_ratio(count-round((c-ed)/step))=length(find(ar<5))/n;
            white_ratio(count-round((c-ed)/step))=length(find(d(:)))/numel(d);
            if white_ratio(count-round((c-ed)/step))>0.09% | dots_ratio(count-round((c-ed)/step))>0.5
                continue
            else
                break
            end
        end
        [white_ratio;dots_ratio]'
        low=1;
        ['Medfilt2 is applied with resize : too weak DIC']
%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif dic==2 & length(find(d(:)))/numel(d)<0.005
        thresh_offset=thresh_offset-5;
%         continue;
    else
        cont=0;
    end
% end


function [a12]=medfilt2_resize(a1, window_size)

%%%%%%%%%%%%%%%%%%%%%%
%
%   a1          : input image
%   window_size : window_size of 2D median filter
%   
%   a12         : median filtered image (same size of a1)
%
%   July 21, 2014
%   
%%%%%%%%%%%%%%%%%%%%%%
    [size_x, size_y]=size(a1);
    x_flag=0; y_flag=0;
    if (size_x/2)~=round(size_x/2)
        a1(end+1,end)=0;
        x_flag=1;
    end
    if (size_y/2)~=round(size_y/2)
        a1(end,end+1)=0;
        y_flag=1;
    end
    a1=imresize(a1,0.5);
    a1m=medfilt2(a1,[window_size window_size]);
    a11=imresize(a1m,2);
    if x_flag==1
        a12=a11(1:end-1,:);
        if y_flag==1
            a12=a12(:,1:end-1);
        end
    elseif x_flag==0
        if  y_flag==1
            a12=a11(:,1:end-1);
        else
            a12=a11;
        end
    end
return


% if size(find(b),1)/numel(a)>0.06 & (dic==0 | dic==1)
%     c=immultiply(a,uint8(b));
%     [x y v]=find(c); clear x y
% 
% %     thre=graythresh(v)*255;
% %     d=(c>thre);
%     hh=hist(double(v(:)),255);
%     hh1=hh(hh>0);
%     mi=find(hh1<max(hh1)/dividend,1,'first');
%     if isempty(mi)
%         thre=ind;
%     else
%         thre=ind+mi;
%     end
%     d=(c>thre);
% elseif  size(find(b),1)/numel(a)>0.2 & dic==2
%     c=immultiply(a,uint8(b));
%     [x y v]=find(c); clear x y
% 
% %     thre=graythresh(v)*255;
% %     d=(c>thre);
%     hh=hist(double(v(:)),255);
%     hh1=hh(hh>0);
%     mi=find(hh1<max(hh1)/5,1,'first');
%     if isempty(mi)
%         thre=ind;
%     else
%         thre=ind+mi;
%     end
%     d=(c>thre);
% else
%     d=b;
% end
