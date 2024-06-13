% This function assigns generation number to each cell
function outy  = genN_unsorted(simiM, simiM2, it,ss,checklive) % definitions as in main code
% ss is the sample size, 0--1 means 0% -- 100%, the percentage of live leaves
% sampled
% create a random sample
% not sorting the cells due to similarity to root
totcell = length(simiM2);

m = nnz(checklive(:,2));
if (m==0)
    error('All empty barcodes');
end
newlea = zeros(m,2);
k = 0;
for i = 1:totcell
    if (checklive(i,2)>0)
        k = k+1;
        newlea(k,:)=checklive(i,:);
    end
end

rs = randsample(m, ceil(ss*m)); % random sample
m = length(rs);


tsimiM2 = zeros(m,1);
tsimiM = zeros(m,m);
y = zeros(m,3);
for i=1:m
    tsimiM2(i) = simiM2(newlea(rs(i),1));
    y(i,1) = newlea(rs(i));
    y(i,3) = tsimiM2(i);
    for j = 1:m
        tsimiM(i,j) = simiM(newlea(rs(i),1),newlea(rs(j),1));
    end
end



y(:,2)=it;
msrel = max(max(tsimiM)); % find the highest relative score
ave = mean(tsimiM2); % find average similarity to the root.
stdsim = std(tsimiM2); % standard deviation of similarity to the root
cohs = ones(m,1); % record the cohort size of each cell that it belongs to
for i = 1:m-1
    for j = i+1:m
        if(tsimiM(i,j)>=0.8*msrel) % set a threshold to group cells
            cohs(i) = cohs(i)+1;
            cohs(j) = cohs(j)+1;
        end
    end
end



for i = 1:m
    if (cohs(i)<0.1*m) % small proportion
        if (tsimiM2(i)>ave+stdsim) % high similarity to the root
            y(i,2)=max(y(i,2)-1,1); % at least generation one
        end
    elseif (cohs(i)>0.5*m) % large proportion
        if(tsimiM2(i)<ave-stdsim) % low similarity to the root
            y(i,2)=y(i,2)+1;
        end
    else
    end
end

outy = y(:,1:2);