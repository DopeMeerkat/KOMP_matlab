n=0;
for i=1:2
    for j=1:3
        for k=1:8
            n=n+1;
            bt(n)=i;
            bn(n)=k;
            section(n)=j;
        end
    end
end

% bt=[1 2];
% bn=[8 7];
% section=[1 2];
% n=length(bt);

% par
for i=1:n
%     if i==5
%         i
%     end
    control_threshold_1(bt(i), bn(i), section(i))
end