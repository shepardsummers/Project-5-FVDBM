clc, clear;

guh = [-2 -3 -4; ...
        2  4 -6; ...
       -1 -2  3; ...
        6  5  3];

guh_1 = [2 2 2; ...
         2 2 2; ...
         2 2 2; ...
         2 2 2];

guh_2 = [8 8 8; ...
         8 8 8; ...
         8 8 8; ...
         8 8 8];

out = zeros(4,3);

tic;
out1 = (max(guh, 0)~=0).*guh_1 + (min(guh, 0)~=0).*guh_2;
time1 = toc;

tic;
for j = 1:4
    for i = 1:3
        if guh(j,i) >= 0
            out(j,i) = guh_1(j,i);
        else
            out(j,i) = guh_2(j,i);
        end
    end
end
time2 = toc;

rat = time2/time1;
% first method is 67% faster


