% this program builds the original tree
function [y,splt, nodct] = btree(celltag) % build tree
%splt records the nodes of splits
leaves = sort(celltag,1);
maxgen = max(celltag(:,2)); % current max generation number
leaind = leaves(:,1); % leaves indice, assume the same generation
m = length(leaind);
nodct = 0; % node count
for i = maxgen-1:-1:1
    flagind = leaind; % indicates if an entry has been taken
    cohort = [];    
    spltct = []; % split nodes count
    for j=1:m
        if (flagind(j)>0)
            res = mod(leaind(j),2^(maxgen-i));
            if (res==0) % happens to be the endpoint of a cohort
                tvec = zeros(1,2^(maxgen-i));
                nodct = nodct + 1;
                tvec(1) = leaind(j);
                cohort = [cohort; sort(tvec)];
                spltct = [spltct; 0]; % this node does not split
                flagind(j)=0;
            else %res>0
                tvec = zeros(1,2^(maxgen-i));
                nodct = nodct + 1;
                tvec(1) = leaind(j);
                flagind(j)=0;
                subentry = 1; % position in this node
                
                    %while (j+subentry<=m && res+subentry<= 2^(maxgen-i) && leaind(j+subentry)-leaind(j)<2^(maxgen-i))
                    while (j+subentry<=m && res+leaind(j+subentry)-leaind(j)<= 2^(maxgen-i))
                        % found a member that belongs to the same cohort
                           subentry = subentry + 1;
                           tvec(subentry) = leaind(j+subentry-1);     
                           flagind(j+subentry-1)=0;
                    end
                if (res>=2^(maxgen-i-1)+1)%no split
                    spltct = [spltct; 0];
                elseif (mod(tvec(subentry),2^(maxgen-i))==0)% the last is an end point of a cohort
                    spltct = [spltct; 1]; % a split
                elseif (mod(tvec(subentry),2^(maxgen-i))>=2^(maxgen-i-1)+1)
                    spltct = [spltct; 1]; % a split
                else % no split
                    spltct = [spltct; 0];
                end
                cohort = [cohort; sort(tvec)];
             end
        end
    end
    y{i}= cohort;
    splt{i}=spltct;
end