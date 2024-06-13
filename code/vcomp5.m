% This function compares the similarity between two vectors, including 0
% entries, need to trim leading and tail 0s
function [y, ost]=vcomp5(s,t) 
ls = length(s);
lt = length(t);
% empty,ACGT=01234, then add 1 to each to get 12345 as index, 0==0 has
% little reward
scm = [0.25 -1.5 -1.5 -1.5 -1.5
      -1.5  1 -2 -2 -2
      -1.5 -2 1 -2 -2
      -1.5  -2 -2 1 -2
      -1.5 -2 -2 -2 1];
mismat = scm(2,3); % penalty of mismatch, nonzeros
mismat0 = scm(1,2);% mismatch penalty involving 0, agree with the scm

rewm = zeros(ls+1,lt+1); % reward matrix, rewm(1,1) means no entry in either s or t
for i=1:lt % s is over, but t has gaps
    rewm(1,i+1)=i*mismat0;
end
for i=1:ls % t is over, but s has gaps
    rewm(i+1,1)=i*mismat0;
end

% DP forward
for i=2:ls+1
    for j=2:lt+1
        rewm(i,j)=max([rewm(i-1,j-1)+scm(s(i-1)+1,t(j-1)+1), rewm(i-1,j)+mismat, rewm(i,j-1)+mismat]);
    end
end
y = rewm(i,j); % highest matching score


% trace back to find the matching string
if (nargout >1)
    lst=ls+lt; % length of st

    st=zeros(lst,2);
    while ((i>1) && (j>1))
        if (rewm(i,j)==rewm(i-1,j-1)+scm(s(i-1)+1,t(j-1)+1))
            st(lst,1)=s(i-1);
            st(lst,2)=t(j-1);
            lst=lst-1;
            i = i-1;
            j = j-1;
        elseif (rewm(i,j)==rewm(i-1,j)+mismat)
            st(lst,1)=s(i-1);
            lst=lst-1;
            i=i-1;
        else
            st(lst,2)=t(j-1);
            lst=lst-1;
            j = j-1;
        end
    end
    % finish up

    if (i>1)
        while (i>1)
            st(lst,1)=s(i-1);
            i = i-1;
            lst=lst-1;
        end
    else
        while (j>1)
            st(lst,2)=t(j-1);
            j = j-1;
            lst=lst-1;
        end
    end
    ost=st(lst+1:end,:);
end
