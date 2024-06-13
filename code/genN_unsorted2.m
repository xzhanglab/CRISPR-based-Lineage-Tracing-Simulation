% This function assigns generation number to each cell
function outy  = genN_unsorted2(simiM2, it,ss,checklive) % definitions as in main code
% ss is the sample size, 0--1 means 0% -- 100%, the percentage of live leaves
% sampled
% create a random sample
% not sorting the cells due to similarity to root
totcell = length(simiM2);

m = nnz(checklive(:,4));
if (m==0)
    error('All empty barcodes');
end
newlea = zeros(m,4);
k = 0;
for i = 1:totcell
    if (checklive(i,4)>0) % live cell
        k = k+1;
        newlea(k,:)=checklive(i,:);
    end
end

rs = randsample(m, ceil(ss*m)); % random sample
m = length(rs);

y = zeros(m,2);
for i=1:m
    y(i,1) = newlea(rs(i));
end



y(:,2)=it;


outy = y(:,1:2);