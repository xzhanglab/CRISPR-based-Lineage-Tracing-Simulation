% This program finds the matches between two matrices
function [y,spy] = findm(v1,v2,spind)
%spind is the splitting node index
% y is total matched nodes, spy is total matched splitting nodes
s1 = size(v1);
s2 = size(v2);
sum1 = sum(v1')';
sum2 = sum(v2')';
nv1 = sortrows([sum1 v1],1);
permsp = sortrows([sum1 spind],1); % permutated splitting index
nv2 = sortrows([sum2 v2],1);
y = 0;
spy = 0;
for i = 1:s1(1)
    cursp = 1; % current starting point
    for j = cursp:s2(1)
        if (nv1(i,1)==nv2(j,1) && nnz(nv1(i,2:end)-nv2(j,2:end))==0) % found one match
            y = y + 1;
            cursp = j;
            if (permsp(i,2)==1) % a splitting node
                spy = spy+1;
            end
            break;
        elseif (nv1(i,1)<nv2(j,1))
            break;
        elseif (nv1(i,1)>nv2(j,1))
            cursp = j;
        else
        end
    end
end
