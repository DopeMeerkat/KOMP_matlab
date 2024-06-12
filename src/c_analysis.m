n=0;
for i=2
    for j=4:8
        for k=1:3
            n=n+1;
            bt(n)=i;
            bn(n)=j;
            section(n)=k;
        end
    end
end

% bt=[1 1 1 2 2 2 2 2 2 2];
% bn=[3 5 7 1 2 3 4 6 6 7];
% section=[2 3 1 2 2 2 3 1 2 1];
% n=length(bt);

% 1 3 2
% 1 5 3
% 1 7 1
% 2 1 2
% 2 2 2
% 2 3 2
% 2 4 3
% 2 6 1
% 2 6 2
% 2 7 1
% par
% parpool(3)
% par
for i=1:n
    control_analysis(bt(i), bn(i), section(i))
end