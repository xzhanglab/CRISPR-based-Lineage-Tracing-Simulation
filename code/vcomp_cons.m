% This function compares the similarity between two vectors with reward for
% consecutive matches, and less penalty for consecutive mismatches
function [y, ost]=vcomp_cons(s,t)
ls = length(s);
lt = length(t);
scm = [1 -2 -2 -2
      -2  1 -2 -2
       -2 -2 1 -2
      -2  -2 -2 1];
mismat = -2; % penalty of mismatch
consmat = 0.2; % reward for consecutive matches
pc = 1; % percentage for further consecutive reward

rewm = zeros(ls+1,lt+1); % reward matrix, rewm(1,1) means no entry in either s or t
matchc = rewm; % match count
mismat_s = rewm; % mismatch in s vector, i.e., 0 counts
mismat_t = rewm; % mismatch in t vector
fr = 0.5; % fraction of penalty on consecutive mismatches

for i=1:lt % s is over, but t has gaps
    rewm(1,i+1)=rewm(1,i)+mismat*(fr)^(i-1);
end
for i=1:ls % t is over, but s has gaps
    rewm(i+1,1)=rewm(i,1)+mismat*(fr)^(i-1);
end

% DP forward
for i=2:ls+1
    for j=2:lt+1
        rewm(i,j)=max([rewm(i-1,j-1)+scm(s(i-1),t(j-1))+consmat*(pc)^matchc(i-1,j-1), rewm(i-1,j)+mismat*(fr)^mismat_t(i-1,j), rewm(i,j-1)+mismat*(fr)^mismat_s(i,j-1)]);
        if (rewm(i,j)==rewm(i-1,j-1)+scm(s(i-1),t(j-1))+consmat*(pc)^matchc(i-1,j-1)) % update consecutive matches count
            if (s(i-1)==t(j-1)) 
                matchc(i,j)=matchc(i-1,j-1)+1; % accumulative
%                 mismat_s(i,j)=0;
%                 mismat_t(i,j)=0;
                %matchc(i,j)=1;
            else
                matchc(i,j)=0;
                mismat_s(i,j)=0;
                mismat_t(i,j)=0;
            end
        elseif (rewm(i,j)==rewm(i-1,j)+mismat*(fr)^mismat_t(i-1,j)) % add one 0 to t
            matchc(i,j)=0;
            mismat_t(i,j)=mismat_t(i-1,j)+1;
            mismat_s(i,j)=0;
        else
            matchc(i,j)=0;
            mismat_s(i,j)=mismat_s(i,j-1)+1;
            mismat_t(i,j)=0;
        end
    end
end
y = rewm(i,j); % highest matching score

% trace back to find the matching string
if (nargout >1)
    lst=ls+lt; % length of st

    st=zeros(lst,2);
    while ((i>1) && (j>1))
        if (rewm(i,j)==rewm(i-1,j-1)+scm(s(i-1),t(j-1))+consmat*pc^matchc(i-1,j-1))
            st(lst,1)=s(i-1);
            st(lst,2)=t(j-1);
            lst=lst-1;
            i = i-1;
            j = j-1;
        elseif (rewm(i,j)==rewm(i-1,j)+mismat*(fr)^mismat_t(i-1,j))
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

