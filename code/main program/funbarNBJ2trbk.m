% This program simulates barcodes evolution with irregular cell division
% and cell death, and with varying crispr-cas9 speed. traceback special
% version
% no filtering by root
function [fy1,fy2] = funbarNBJ2trbk(n,it,propm,propmi,ss,mupb,ins_sub,lgdelprob,divp,clive, pulse, trbk)

if (trbk>it-1)
    error('wrong trbk number');
end


ldtl = ceil(1*n); % lead and tail space on the ends of barcode
% assign a large negative number to indicate empty barcodes
empbar = -4*(n+2*ldtl);

barM = zeros(n+2*ldtl,2^(it+1)-1); % barcode matrix
barM(ldtl+1:ldtl+n,1)=ceil(4*rand(n,1)); % initial barcode A,C,G,T = 1,2,3,4
sted = ones(2,2^(it+1)-1); % record start and end position of each barcode
                            % mostly due to insertion
sted(1,:)=ldtl+1;
sted(2,:) = ldtl+n;

%-------------------repeat for barM2-------------------------------------
barM2 = zeros(n+2*ldtl,2^(it+1)-1); % barcode matrix
barM2(ldtl+1:ldtl+n,1)=ceil(4*rand(n,1)); % initial barcode A,C,G,T = 1,2,3,4
sted2 = ones(2,2^(it+1)-1); % record start and end position of each barcode
                            % mostly due to insertion
sted2(1,:)=ldtl+1;
sted2(2,:) = ldtl+n;
%-----------------------------------------------------------------------
inmupb = mupb;

% start mutation, construct the evolution tree-----------------------------
for k=1:it    % cell division iteration = generation number
    simtree{k}=zeros(2^k,4); % (cell index in this generation, 1st barcode generation number, 2nd barcode generation number, live cell indicator)
    simtree{k}(:,1)=1:2^k;  % generation number=0 means dead barcode, if both 1st and 2nd gen number=0, dead cell-set live indicator 0
    
    % test pulse dox concentration
    
	if (pulse == 1)
     		if (mod(k,2)==0)
         		mupb = 0.005; % base cut rate, suppose no dox
     		else
         		mupb = inmupb; % with pulses of dox
     		end
	end

%     if (k>=it-trbk)
%         mupb = inmupb; % with constant dox application
%     else
%         mupb = 0.005; % base cut rate, suppose no dox
%     end
    
    if (k==1) % no need to check parent's life, assume the root is a live cell
            for j=1:2^(k-1) 
                % first child---------------- 
                    rliv = rand; %probability of a live child
                    if (rliv<clive) % live child, generate it
                        [barM(:,2^k+2*(j-1)), sted(1,2^k+2*(j-1)), sted(2,2^k+2*(j-1))]=genchild(n,ldtl,sted(1,2^(k-1)+j-1),sted(2,2^(k-1)+j-1),barM(:,2^(k-1)+j-1),mupb,lgdelprob,ins_sub);
                        if (nnz(barM(:,2^k+2*(j-1)))>0) % live cell
                            simtree{k}(2*j-1,2)=k;
                        end
                        simtree{k}(2*j-1,4)=1;
                    end
                    % 2nd child-------------------------------------------    
                    rliv = rand; %probability of a live child
                    if (rliv<clive) % live child, generate it
                        [barM(:,2^k+2*j-1), sted(1,2^k+2*j-1), sted(2,2^k+2*j-1)]=genchild(n,ldtl,sted(1,2^(k-1)+j-1),sted(2,2^(k-1)+j-1),barM(:,2^(k-1)+j-1),mupb,lgdelprob,ins_sub);
                        if (nnz(barM(:,2^k+2*j-1))>0) % live cell
                            simtree{k}(2*j,2)=k;
                        end
                        simtree{k}(2*j,4)=1;
                    end
            end % j for this generation
    else %k>1, need to check parent's life
        for j=1:2^(k-1) 
            if (simtree{k-1}(j,4)>0 && simtree{k-1}(j,2)>0) % cell is live and 1st barcode is live in parent generation
                rdiv = rand; % probability of division
                if (rdiv<divp) % division occurs
                    % first child---------------- 
                    rliv = rand; %probability of a live child
                    if (rliv<clive) % live child, generate it
                        [barM(:,2^k+2*(j-1)), sted(1,2^k+2*(j-1)), sted(2,2^k+2*(j-1))]=genchild(n,ldtl,sted(1,2^(k-1)+j-1),sted(2,2^(k-1)+j-1),barM(:,2^(k-1)+j-1),mupb,lgdelprob,ins_sub);
                        if (nnz(barM(:,2^k+2*(j-1)))>0) % live cell
                            simtree{k}(2*j-1,2)=k;
                        end
                        simtree{k}(2*j-1,4)=1;
                    end
                    % 2nd child-------------------------------------------    
                    rliv = rand; %probability of a live child
                    if (rliv<clive) % live child, generate it
                        [barM(:,2^k+2*j-1), sted(1,2^k+2*j-1), sted(2,2^k+2*j-1)]=genchild(n,ldtl,sted(1,2^(k-1)+j-1),sted(2,2^(k-1)+j-1),barM(:,2^(k-1)+j-1),mupb,lgdelprob,ins_sub);
                        if (nnz(barM(:,2^k+2*j-1))>0) % live cell
                            simtree{k}(2*j,2)=k;
                        end
                        simtree{k}(2*j,4)=1;
                    end
                else % single child
                     % only child---------------- 
                    rliv = rand; %probability of a live child
                    if (rliv<clive) % live child, generate it
                        [barM(:,2^k+2*(j-1)), sted(1,2^k+2*(j-1)), sted(2,2^k+2*(j-1))]=genchild(n,ldtl,sted(1,2^(k-1)+j-1),sted(2,2^(k-1)+j-1),barM(:,2^(k-1)+j-1),mupb,lgdelprob,ins_sub);
                        if (nnz(barM(:,2^k+2*(j-1)))>0) % live cell
                            simtree{k}(2*j-1,2)=k;
                        end
                        simtree{k}(2*j-1,4)=1;
                    end
                end
                
            end
        end
    end% check parent life
end % i for division round

% draw barcodes of the entire simulation

% figure
% heatmap(barM);
% xlabel('Cell Label');

%--------------------repeat for
%barM2---------------------------------------------------------------------
% start mutation, construct the evolution tree-----------------------------
for k=1:it    % cell division iteration = generation number
    
    % test pulse dox concentration
    
% 	if (pulse == 1)
%      		if (mod(k,2)==0)
%          		mupb = 0.005; % base cut rate, suppose no dox
%      		else
%          		mupb = inmupb; % with pulses of dox
%      		end
%     end

    if (k>=it-trbk)
        mupb = inmupb; % with constant dox application
    else
        mupb = 0.005; % base cut rate, suppose no dox
    end
    
    if (k==1) % no need to check parent's life, assume the root is a live cell
            for j=1:2^(k-1) 
                % first child---------------- 
                    
                    if (simtree{k}(2*j-1,4)>0) % live child, generate it
                        [barM2(:,2^k+2*(j-1)), sted2(1,2^k+2*(j-1)), sted2(2,2^k+2*(j-1))]=genchild(n,ldtl,sted2(1,2^(k-1)+j-1),sted2(2,2^(k-1)+j-1),barM2(:,2^(k-1)+j-1),mupb,lgdelprob,ins_sub);
                        if (nnz(barM2(:,2^k+2*(j-1)))>0) % live cell
                            simtree{k}(2*j-1,3)=k;
                        end
                    end
                    % 2nd child-------------------------------------------    
                    
                    if (simtree{k}(2*j,4)>0) % live child, generate it
                        [barM2(:,2^k+2*j-1), sted2(1,2^k+2*j-1), sted2(2,2^k+2*j-1)]=genchild(n,ldtl,sted2(1,2^(k-1)+j-1),sted2(2,2^(k-1)+j-1),barM2(:,2^(k-1)+j-1),mupb,lgdelprob,ins_sub);
                        if (nnz(barM2(:,2^k+2*j-1))>0) % live cell
                            simtree{k}(2*j,3)=k;
                        end
                    end
            end % j for this generation
    else %k>1, need to check parent's life
        for j=1:2^(k-1) 
            if (simtree{k-1}(j,4)>0 && simtree{k-1}(j,3)>0) % cell is live and 2nd barcode is live in parent generation
                

                    % first child---------------- 
                   
                    if (simtree{k}(2*j-1,4)>0) % live child, generate it
                        [barM2(:,2^k+2*(j-1)), sted2(1,2^k+2*(j-1)), sted2(2,2^k+2*(j-1))]=genchild(n,ldtl,sted2(1,2^(k-1)+j-1),sted2(2,2^(k-1)+j-1),barM2(:,2^(k-1)+j-1),mupb,lgdelprob,ins_sub);
                        if (nnz(barM2(:,2^k+2*(j-1)))>0) % live cell
                            simtree{k}(2*j-1,3)=k;
                        end
                    end
                    % 2nd child-------------------------------------------    
  
                    if (simtree{k}(2*j,4)>0) % live child, generate it
                        [barM2(:,2^k+2*j-1), sted2(1,2^k+2*j-1), sted2(2,2^k+2*j-1)]=genchild(n,ldtl,sted2(1,2^(k-1)+j-1),sted2(2,2^(k-1)+j-1),barM2(:,2^(k-1)+j-1),mupb,lgdelprob,ins_sub);
                        if (nnz(barM2(:,2^k+2*j-1))>0) % live cell
                            simtree{k}(2*j,3)=k;
                        end
                    end

                
            end
        end
    end% check parent life
end % i for division round
%---------------------------------------------------------------------------------

% draw barcodes of the entire simulation

% figure
% heatmap(barM);
% xlabel('Cell Label');
% 
% figure
% heatmap(barM2);
% xlabel('Cell Label');

%reverse alignment between root and leaves, consecutive matches preferred-----------------------
revtrac = zeros(n+2*ldtl, 2^it+1); % reverse trace
revtrac(:,1)=barM(:,1); % show the root/first cell barcode
simiM2 = empbar*ones(2^it,1);
tvec1 = barM(1+ldtl:n+ldtl,1);
colv1 = tvec1(tvec1~=0);
for i=1:2^it
    if(simtree{it}(i,4)>0) % this is a live cell
        % collapse columns        
        tvec2 = barM(sted(1,2^it+i-1):sted(2,2^it+i-1),2^it+i-1);        
        colv2 = tvec2(tvec2~=0);
        
        if (isempty(colv2))% 
            simtree{it}(i,2)=0; % 1st barcode is empty
            %warning('empty barcode');
            %i
            %simiM2(i) = empbar;
        else % compare
           [y, outv]=vcomp_cons(colv1,colv2);
           %revtrac(ldtl+1:ldtl+length(outv(:,2)),i+1)=outv(:,2); % this
           %could cause over length
           
           revtrac(1:length(outv(:,2)),i+1)=outv(:,2);
           
           simiM2(i) = y;
           %clear outv;
        end
       %clear tvec2 colv2;
    end
end

% ------------repeat for barM2----------------------------------
revtrac2 = zeros(n+2*ldtl, 2^it+1); % reverse trace
revtrac2(:,1)=barM2(:,1); % show the root/first cell barcode
simiM22 = empbar*ones(2^it,1); % for second barcode
tvec1 = barM2(1+ldtl:n+ldtl,1);
colv1 = tvec1(tvec1~=0);
for i=1:2^it
    if(simtree{it}(i,4)>0) % this is a live cell
        % collapse columns        
        tvec2 = barM2(sted(1,2^it+i-1):sted(2,2^it+i-1),2^it+i-1);        
        colv2 = tvec2(tvec2~=0);
        
        if (isempty(colv2))% 
            simtree{it}(i,3)=0; % 2nd barcode is empty
            %warning('empty barcode');
            %i
            %simiM2(i) = empbar;
        else % compare
           [y, outv]=vcomp_cons(colv1,colv2);
           %revtrac(ldtl+1:ldtl+length(outv(:,2)),i+1)=outv(:,2); % this
           %could cause over length
           
           revtrac2(1:length(outv(:,2)),i+1)=outv(:,2);
           
           simiM22(i) = y;
           %clear outv;
        end
        if (simtree{it}(i,2)==0 && simtree{it}(i,3)==0) % both barcodes are empty, mark dead
            simtree{it}(i,4)=0;
        end
       %clear tvec2 colv2;
    end
end


% % draw reconstructed leaves
% figure
% heatmap(revtrac);
% title('Reverse alignment consecutive matches preferred')
% 
% figure
% heatmap(revtrac2);
% title('Reverse alignment consecutive matches preferred')


%clear tvec1 colv1;

% figure
% plot(simiM2,'*');
% title('Matching scores: Root vs leaves');

% compare similarity between cells in the last generation, leaves------------
simiM = empbar*ones(2^it,2^it); 
for i=1:2^it
    if (simtree{it}(i,4)>0) % live cell
        for j=i+1:2^it
            if (simtree{it}(j,4)>0) % live cell
                tvec1 = revtrac(:,i+1);
                tvec2 = revtrac(:,j+1);
                % trim leading and tail 0s--------
                    colv1 = trim(tvec1,2*ldtl+n);
                        %----------------                
                    colv2 = trim(tvec2,2*ldtl+n);

                % compare
                simiM(i,j)=vcomp5(colv1,colv2);
                %clear tvec1 tvec2 colv1 colv2;
                simiM(j,i)=simiM(i,j);
            end
        end
    end
end

%--repeat for 2nd barcode---------------------
simiMM = empbar*ones(2^it,2^it); 
for i=1:2^it
    if (simtree{it}(i,4)>0) % live cell
        for j=i+1:2^it
            if (simtree{it}(j,4)>0) % live cell
                tvec1 = revtrac2(:,i+1);
                tvec2 = revtrac2(:,j+1);
                % trim leading and tail 0s--------
                    colv1 = trim(tvec1,2*ldtl+n);
                        %----------------                
                    colv2 = trim(tvec2,2*ldtl+n);

                % compare
                simiMM(i,j)=vcomp5(colv1,colv2);
                %clear tvec1 tvec2 colv1 colv2;
                simiMM(j,i)=simiMM(i,j);
            end
        end
    end
end


% make a symmetric matrix
% simiM=simiM+simiM'-eye(size(simiM)).*diag(simiM);
% simiM=simiM+simiM';
% for i=1:2^it
%     simiM(i,i)=empbar;
% end


% figure
% heatmap(simiM);
% title('Similarity among leaves');
% 
% figure
% heatmap(simiMM);
% title('Similarity among leaves');


% reconstruct lineage tree-------------------------------
celltag = genN_unsorted2(simiM2, it,ss,simtree{it}); % generate a random sample
indcell = celltag(:,1); % contains all sampled live cells
%--------------------------------------------------------------------
valflag = valcheck(celltag);
if (valflag == 0) % not valid assignment
    error('Not valid assignment of generation numbers');
end

%---------------------build the original tree-----------
[otree, splt, ~] = btree(celltag); % original tree and number of nodes
%------------------------------------------------------

maxgen = max(celltag(:,2)); % current max generation number
m = length(indcell);

% record the barcodes of leaves
tbarC = zeros(n+2*ldtl,m);
for i = 1:m
    tbarC(:,i) = revtrac(:,indcell(i)+1); % record the sampled (not necessarily all) barcodes of the 'last' generation, i.e., available information at the beginning
                            % reorder the cells in tbarC
end

%-----------repeat for 2nd barcode-------------

tbarC2 = zeros(n+2*ldtl,m);
for i = 1:m
    tbarC2(:,i) = revtrac2(:,indcell(i)+1); % record the sampled (not necessarily all) barcodes of the 'last' generation, i.e., available information at the beginning
                            % reorder the cells in tbarC
end

% 
% figure
% heatmap(tbarC);
% 
% figure
% heatmap(tbarC2);

%lineage = zeros(length(indcell),2*(maxgen-1)); % both columns are nonempty=paired, left column nonempty = singleton, right column nonempty = unchanged.
C={}; % use cell array to store lineage
indice={}; % record cell orders
while maxgen > 1 % start rebuilding the binomial tree---------------------
    
%     figure
%     heatmap(tbarC);

    indice{maxgen-1}=indcell;
    
    m=length(indcell);
    
    nsimiM = empbar*ones(m,m);
    nsimiMM = nsimiM;

    for i = 1:m % construct pairwise similarity
        for j = i+1:m
            if (maxgen == it)
                nsimiM(i,j)=simiM(indcell(i),indcell(j));
                nsimiMM(i,j)=simiMM(indcell(i),indcell(j));
            else
                if (celltag(i,2)==maxgen && celltag(j,2)==maxgen)
                    tvec1 = trim(tbarC(:,i),2*ldtl+n);
                    colv1 = tvec1(tvec1~=0);
                    tvec2 = trim(tbarC(:,j),2*ldtl+n);
                    colv2 = tvec2(tvec2~=0);
                    if (~isempty(colv1) && ~isempty(colv2))
                        nsimiM(i,j) = vcomp5(colv1,colv2);
                    end
                    
                    %----repeat for 2nd barcode----------
                    tvec1 = trim(tbarC2(:,i),2*ldtl+n);
                    colv1 = tvec1(tvec1~=0);
                    tvec2 = trim(tbarC2(:,j),2*ldtl+n);
                    colv2 = tvec2(tvec2~=0);
                    
                    if (~isempty(colv1) && ~isempty(colv2))
                        nsimiMM(i,j) = vcomp5(colv1,colv2);
                    end
                end
            end
        end
    end

   TsimiM = nsimiM+nsimiMM; % total similarity of both matrices

if (maxgen < it)
    simiM = nsimiM; % record the original matrices, no need to recalculate later
    simiMM = nsimiMM;
end

    paired = [];
    lv = m; % length of vector
    ncelltag = []; % new cell tag
    ncount = 0; % new cell count
    singleton = [];
    unchanged = []; % unchanged barcodes due to lower generation


    
% ----start from the total similarity----------------
    [Mc,maxind] = max(TsimiM); % max of each column of TsimiM

    [M,maxcol] = max(Mc);
    while (M>0)  % set a threshold to stop searching in the matrix
        if (celltag(maxcol,2)==maxgen && celltag(maxind(maxcol),2)==maxgen) % both belong to the last generation
            tvec1 = trim(tbarC(:,maxcol),2*ldtl+n); % 1st barcode
            colv1 = tvec1(tvec1~=0);
            tvec2 = trim(tbarC(:,maxind(maxcol)),2*ldtl+n);
            colv2 = tvec2(tvec2~=0);

            tvec12 = trim(tbarC2(:,maxcol),2*ldtl+n);% 2nd barcode
            colv12 = tvec12(tvec12~=0);
            tvec22 = trim(tbarC2(:,maxind(maxcol)),2*ldtl+n);
            colv22 = tvec22(tvec22~=0);

            if (~isempty(colv1) && ~isempty(colv2) && ~isempty(colv12) && ~isempty(colv22)) % all 4 barcodes are nonempty barcodes
                if (M>2*max([propm*nnz(colv1),propm*nnz(colv2),propm*nnz(colv12),propm*nnz(colv22)])) % paired
                    paired = [paired; maxcol maxind(maxcol)];
                    ncount = ncount + 1;
                    ncelltag = [ncelltag; ncount maxgen-1];
                    celltag(maxcol,2)=0; % remove cell i
                    celltag(maxind(maxcol),2)=0; % remove cell curpar
                    nsimiM(:,maxcol)=empbar;
                    nsimiM(maxind(maxcol),:)=empbar;
                    nsimiMM(:,maxcol)=empbar;
                    nsimiMM(maxind(maxcol),:)=empbar;
                    TsimiM(:,maxcol)=2*empbar;
                    TsimiM(maxind(maxcol),:)=2*empbar;
                else % the matching score is lower than propm*nnz, so do not pair
                     TsimiM(maxind(maxcol),maxcol)=empbar; % here, update TsimiM only, not individual nsimiM or nsimiMM
                end
            else % there is at least one empty barcode
                TsimiM(maxind(maxcol),maxcol)=empbar; % here, update TsimiM only, not individual nsimiM or nsimiMM
            end
        else % the two barcodes do not belong to the same generation
            TsimiM(maxind(maxcol),maxcol)=empbar; % here, update TsimiM only, not individual nsimiM or nsimiMM
        end

        % update the new M
        [Mc,maxind] = max(TsimiM); % max of each column
        [M,maxcol] = max(Mc);
    end
   
    %-----------for individual barcode---------------------------
    [Mc1,maxind1] = max(nsimiM); % max of each column of nsimiM
    [M1,maxcol1] = max(Mc1);

    [Mc2,maxind2] = max(nsimiMM); % max of each column of nsimiMM
    [M2,maxcol2] = max(Mc2);
    M=max(M1,M2);
     while (M>0.1*empbar)  % set a threshold to stop searching in the matrix
         if (M1>M2)
                if (celltag(maxcol1,2)==maxgen && celltag(maxind1(maxcol1),2)==maxgen) % both belong to the last generation
                    tvec1 = trim(tbarC(:,maxcol1),2*ldtl+n);
                    colv1 = tvec1(tvec1~=0);
                    tvec2 = trim(tbarC(:,maxind1(maxcol1)),2*ldtl+n);
                    colv2 = tvec2(tvec2~=0);
                    if (~isempty(colv1) && ~isempty(colv2)) % both are nonempty barcodes
                        if (M>max(propmi*nnz(colv1),propmi*nnz(colv2))) % paired
                            paired = [paired; maxcol1 maxind1(maxcol1)];
                            ncount = ncount + 1;
                            ncelltag = [ncelltag; ncount maxgen-1];
                            celltag(maxcol1,2)=0; % remove cell i
                            celltag(maxind1(maxcol1),2)=0; % remove cell curpar
                            nsimiM(:,maxcol1)=empbar;
                            nsimiM(maxind1(maxcol1),:)=empbar;
                            nsimiMM(:,maxcol1)=empbar;
                            nsimiMM(maxind1(maxcol1),:)=empbar;
                        else % the matching score is lower than propm*nnz, so do not pair
                             nsimiM(maxind1(maxcol1),maxcol1)=empbar;
                        end
                    else % there is at least one empty barcode
                        nsimiM(maxind1(maxcol1),maxcol1)=empbar; % update nsimiM only
                    end
                else % the two barcodes do not belong to the same generation
                    nsimiM(maxind1(maxcol1),maxcol1)=empbar;
                end
                [Mc1,maxind1] = max(nsimiM); % max of each column of TsimiM
                [M1,maxcol1] = max(Mc1);
         else % 2nd barcode is used
                if (celltag(maxcol2,2)==maxgen && celltag(maxind2(maxcol2),2)==maxgen) % both belong to the last generation
                    tvec1 = trim(tbarC2(:,maxcol2),2*ldtl+n);
                    colv1 = tvec1(tvec1~=0);
                    tvec2 = trim(tbarC2(:,maxind2(maxcol2)),2*ldtl+n);
                    colv2 = tvec2(tvec2~=0);
                    if (~isempty(colv1) && ~isempty(colv2)) % both are nonempty barcodes
                        if (M>max(propmi*nnz(colv1),propmi*nnz(colv2))) % paired
                            paired = [paired; maxcol2 maxind2(maxcol2)];
                            ncount = ncount + 1;
                            ncelltag = [ncelltag; ncount maxgen-1];
                            celltag(maxcol2,2)=0; % remove cell i
                            celltag(maxind2(maxcol2),2)=0; % remove cell curpar
                            nsimiM(:,maxcol2)=empbar;
                            nsimiM(maxind2(maxcol2),:)=empbar;
                            nsimiMM(:,maxcol2)=empbar;
                            nsimiMM(maxind2(maxcol2),:)=empbar;
                        else % the matching score is lower than propm*nnz, so do not pair
                             nsimiMM(maxind2(maxcol2),maxcol2)=empbar;
                        end
                    else % there is at least one empty barcode
                        nsimiMM(maxind2(maxcol2),maxcol2)=empbar; % update nsimiMM only
                    end
                else % the two barcodes do not belong to the same generation
                    nsimiMM(maxind2(maxcol2),maxcol2)=empbar;
                end
                [Mc2,maxind2] = max(nsimiMM); % max of each column of TsimiM
                [M2,maxcol2] = max(Mc2);
         end
            % update the new M                 
                
                M=max(M1,M2);
        
     end

     
    for i=1:lv % collect the remaining cells
        if (celltag(i,2)==maxgen) % singleton
            singleton = [singleton; i];% collect singleton barcodes, even if it is empty barcode

        elseif (celltag(i,2)>0) % unchanged
            tvec1 = trim(tbarC(:,i),2*ldtl+n);
            colv1 = tvec1(tvec1~=0);
            tvec2 = trim(tbarC2(:,i),2*ldtl+n);
            colv2 = tvec1(tvec1~=0);
            if (~isempty(colv1) || ~isempty(colv2))
                unchanged = [unchanged; i celltag(i,2)]; % here the index is new
                warning('you have unchanged barcodes due to different generation number');
            else
                error('Empty unchanged barcode');
            end            
        end
    end
    

    
    % now ncount is the length of paired
    ctunch = size(unchanged,1); % length of unchanged

    %check validity of the new generation
    TsimiM = simiM+simiMM; % total similarity of both matrices
    lv = length(singleton);
    if (lv>0)
        tpars = [ncelltag; unchanged; zeros(lv,2)];
        tpars(ncount+ctunch+1:end,2) = maxgen-1;
        while (~valcheck(tpars) && lv>1)% cell numbers in each generation do not pass validity check, make one pair

           

            nsimiM = 2*empbar*ones(lv,lv);
        
           for i = 1:lv
               for j = i+1:lv
                    if (maxgen == it)
                        nsimiM(i,j)=TsimiM(indcell(singleton(i)),indcell(singleton(j)));
                    else
                        nsimiM(i,j) = TsimiM(singleton(i),singleton(j));
                    end
                end
            end
            [Mc,maxind] = max(nsimiM); % max of each column
            [~,maxcol] = max(Mc);
            paired = [paired; singleton(maxcol) singleton(maxind(maxcol))];
            ncount = ncount + 1;
            ncelltag = [ncelltag; ncount maxgen-1];
            celltag(singleton(maxcol),2)=0; % remove cell i
            celltag(singleton(maxind(maxcol)),2)=0; % remove cell curpar
%             nsimiM(:,maxcol)=empbar;
%             nsimiM(maxind(maxcol),:)=empbar;

            remc = 0; % remaing cell count
            newsing = zeros(lv-2,1); % new singleton cell list

            for j=1:lv
                if ( celltag(singleton(j),2)>0)
                    remc = remc+1;
                    newsing(remc) = singleton(j);
                end
            end
            lv = lv-2;
            tpars = [ncelltag; unchanged; zeros(lv,2)];
            tpars(ncount+ctunch+1:end,2) = maxgen-1;      
            singleton = newsing; % update singleton list
            newsing = [];
        end
    end
  
   
    % passed validity check
    
    for i=1:ctunch% add unchanged cells, now ncelltag has paired and unchanged
        ncount = ncount + 1;
        ncelltag = [ncelltag; ncount unchanged(i,2)];
    end
    % now ncount is the sum of paired and unchanged
    

    % need to rebuild parent nodes and compare to root-----------------

    thisgen = zeros(2*ldtl+n,ncount+lv); % barcodes in this iteration
    % it consists of three parts: paired, unchanged yet, and singleton
    tsize = size(paired,1);
    curlin = zeros(ncount+lv,2); % current lineage

    for i = 1:tsize % paired--------------------part I
        % extract lineage information
        curlin(i,1)=paired(i,1);
        curlin(i,2)=paired(i,2);

        colv1 = trim(tbarC(:,paired(i,1)),2*ldtl+n);
        colv2 = trim(tbarC(:,paired(i,2)),2*ldtl+n);
        [~, outv]=vcomp5(colv1,colv2);
        tlength = size(outv,1);
        if (tlength>2*ldtl+n)
            error('Too long barcode, check vcomp5.');
        end
        parent = zeros(tlength,1);
            for j = 1:tlength % reconstruct temporary parent node
                if (outv(j,1)==outv(j,2))
                    parent(j) = outv(j,1);
                elseif (rand<0.5)
                    parent(j) = outv(j,1);
                else
                    parent(j) = outv(j,2);
                end
            end

            % try a different way of alignment-collapse and redo vcomp_cns
            colv1 = parent(parent~=0); 
            if (~isempty(colv1)) %  rebuild the parent node only if the child is not empty
                
%                 % now colv1 is the temporary parent
%                 [~, outv]=vcomp_cons(barM(ldtl+1:ldtl+n,1),colv1); % compare to root
%                 
%                 tlength = size(outv,1);
%                 parent = outv(:,2);
%                 for j = 1:tlength % reconstruct parent node
%                     if (outv(j,1)~=parent(j)) % with a certain probability to replace entry with the root entry
%                         if(rand<1/maxgen)
%                             parent(j)=outv(j,1);
%                         end
%                     end
%                 end
                
                %thisgen(:,i) =  [zeros(ldtl,1);parent;zeros(n+ldtl-tlength(1),1)];
                %parent = trim(parent,2*ldtl+n);
                parent = trim(colv1,2*ldtl+n);
                tlength = length(parent);
                thisgen(:,i) =  [parent;zeros(n+2*ldtl-tlength,1)]; % fill up to make equal length
            end
    end % end paired part I---------------------------------------------

    %======repeat for 2nd barcode---------------------------------
    thisgen2 = zeros(2*ldtl+n,ncount+lv); % barcodes in this iteration
    for i = 1:tsize % paired--------------------part I
        % extract lineage information

        colv1 = trim(tbarC2(:,paired(i,1)),2*ldtl+n);
        colv2 = trim(tbarC2(:,paired(i,2)),2*ldtl+n);
        [~, outv]=vcomp5(colv1,colv2);
        tlength = size(outv,1);
        if (tlength>2*ldtl+n)
            error('Too long barcode, check vcomp5.');
        end
        parent = zeros(tlength,1);
            for j = 1:tlength % reconstruct temporary parent node
                if (outv(j,1)==outv(j,2))
                    parent(j) = outv(j,1);
                elseif (rand<0.5)
                    parent(j) = outv(j,1);
                else
                    parent(j) = outv(j,2);
                end
            end

            % try a different way of alignment-collapse and redo vcomp_cns
            colv1 = parent(parent~=0); 
            if (~isempty(colv1))
               
                % now colv1 is the temporary parent
%                 [~, outv]=vcomp_cons(barM2(ldtl+1:ldtl+n,1),colv1); % compare to root
%                 
%                 tlength = size(outv,1);
%                 parent = outv(:,2);
%                 for j = 1:tlength % reconstruct parent node
%                     if (outv(j,1)~=parent(j)) % with a certain probability to replace entry with the root entry
%                         if(rand<1/maxgen)
%                             parent(j)=outv(j,1);
%                         end
%                     end
%                 end
                %thisgen(:,i) =  [zeros(ldtl,1);parent;zeros(n+ldtl-tlength(1),1)];
                %parent = trim(parent,2*ldtl+n);
                parent = trim(colv1,2*ldtl+n);
                tlength = length(parent);
                thisgen2(:,i) =  [parent;zeros(n+2*ldtl-tlength,1)]; % fill up to make equal length
            end
    end % end paired part I---------------------------------------------

    
    % unchanged barcodes part II--------------------------------------
    for i=tsize+1:ncount
        % extract information, this cell does not belong to this
        % generation, leave it empty
        curlin(i,2)=unchanged(i-tsize,1);
        
        parent = trim(tbarC(:,unchanged(i-tsize,1)),2*ldtl+n);
        tlength = size(parent,1);
        %thisgen(:,i) =  [zeros(ldtl,1);parent;zeros(n+ldtl-tlength(1),1)];
        thisgen(:,i) =  [parent;zeros(n+2*ldtl-tlength,1)];
    end
    
%------repeat for 2nd barcode---------------------------
    % unchanged barcodes part II--------------------------------------
    for i=tsize+1:ncount
        % extract information, this cell does not belong to this
        % generation, leave it empty
       
        parent = trim(tbarC2(:,unchanged(i-tsize,1)),2*ldtl+n);
        tlength = size(parent,1);
        %thisgen(:,i) =  [zeros(ldtl,1);parent;zeros(n+ldtl-tlength(1),1)];
        thisgen2(:,i) =  [parent;zeros(n+2*ldtl-tlength,1)];
    end


    % singleton barcodes part III-------------------------------------     
    if (lv>0) % there are remaining singleton cells
        for i = 1:lv
            ncount = ncount + 1;

            curlin(ncount,1)=singleton(i);% extract information
            ncelltag = [ncelltag; ncount maxgen-1];
%            colv1 = trim(tbarC(:,indcell(singleton(i))),2*ldtl+n);

            colv1 = trim(tbarC(:,singleton(i)),2*ldtl+n);            
            colv2 = colv1(colv1~=0);
%             [~, outv]=vcomp5(barM(ldtl+1:ldtl+n,1),colv1); % compare to root
            if (~isempty(colv2))
%                 [~, outv]=vcomp_cons(barM(ldtl+1:ldtl+n,1),colv2); % compare to root
%                 
%                 tlength = size(outv,1);
%                 parent = outv(:,2);
%                 for j = 1:tlength % reconstruct parent node
%                     if (outv(j,1)~=parent(j)) % with a certain probability to replace entry with the root entry
%                         if(rand<1/maxgen)
%                             parent(j)=outv(j,1);
%                         end
%                     end
%                 end
                %thisgen(:,ncount) =
                %[zeros(ldtl,1);parent;zeros(n+ldtl-tlength(1),1)];
                %parent=colv1;
                parent=colv2;
                tlength = length(parent);
                thisgen(:,ncount) = [parent;zeros(n+2*ldtl-tlength,1)];
            end

      % repeat for 2nd barcode-----------------------------------   
            colv1 = trim(tbarC2(:,singleton(i)),2*ldtl+n);            
            colv2 = colv1(colv1~=0);
            %             [~, outv]=vcomp5(barM(ldtl+1:ldtl+n,1),colv1); % compare to root
            if (~isempty(colv2))
%                 [~, outv]=vcomp_cons(barM2(ldtl+1:ldtl+n,1),colv2); % compare to root
%                 
%                 tlength = size(outv,1);
%                 parent = outv(:,2);
%                 for j = 1:tlength % reconstruct parent node
%                     if (outv(j,1)~=parent(j)) % with a certain probability to replace entry with the root entry
%                         if(rand<1/maxgen)
%                             parent(j)=outv(j,1);
%                         end
%                     end
%                 end
                %thisgen(:,ncount) = [zeros(ldtl,1);parent;zeros(n+ldtl-tlength(1),1)];
                %parent=colv1;
                parent=colv2;
                tlength = length(parent);
                thisgen2(:,ncount) = [parent;zeros(n+2*ldtl-tlength,1)];
            end
        end        
    end
      
%     %-----check rebuilding parent barcode
%     figure
%     heatmap(tbarC);
%     xlabel('1st');
% 
%     figure
%     heatmap(tbarC2);
%     xlabel('2nd');


    celltag = zeros(ncount,2); % update cell tags------------------
    celltag=ncelltag;
%     celltag(:,1) = indcell;
    tbarC = zeros(n+2*ldtl,ncount);
    tbarC = thisgen;

    tbarC2 = zeros(n+2*ldtl,ncount);
    tbarC2 = thisgen2;

    indcell = [1:ncount]';

    maxgen = maxgen - 1;
    %tbarC = thisgen;
    C{maxgen}=curlin;    
end

%C

% build new lineage tree

for i = it-1:-1:it-trbk % this is the dimension of C, C{i} consists of 2 columns, do backward
    clen = length(C{i}(:,1)); % length of column in C{i}
    if (i==it-1)% leaves
        tempvec=zeros(clen,2); % convert columns of C into vectors, need to find its original cell index
        for j = 1:clen
            % convert to its original cell index, which contains the
            % topology of the true tree
            if (C{i}(j,1)>0)
                tempvec(j,1)=indice{i}(C{i}(j,1)); 
            end
            if (C{i}(j,2)>0)
                tempvec(j,2)=indice{i}(C{i}(j,2)); 
            end
            ntree{i}(j,:)=sort(tempvec(j,:)); % newtree
        end
        
    else % internal
        newC=zeros(clen,2^(it-i)); % record all the decendants of each node in its original index, i.e., clade
        for j = 1:clen
            if (C{i}(j,1)>0 && C{i}(j,2)>0) % a pair, combine the decendants
                newC(j,:)=sort([tempvec(indice{i}(C{i}(j,1)),:) tempvec(indice{i}(C{i}(j,2)),:)]);               
            end
            if (C{i}(j,1)==0) %unchanged cells
                newC(j,:)=sort([tempvec(indice{i}(C{i}(j,2)),:) zeros(1,2^(it-i-1))]);
            end
            if (C{i}(j,2)==0) %singleton cells
                newC(j,:)=sort([tempvec(indice{i}(C{i}(j,1)),:) zeros(1,2^(it-i-1))]); % to make it a vector with full length
            end
            
        end
        tempvec = newC;
        ntree{i} = newC;
    end
end

% calculate similarity score between two trees----------------------------
comptree = 0; % number of correct nodes, i.e., correct parent with correct decendants/clade
splitm = 0; % splitting nodes match, 

nodctnew = 0; % new node counts, depending on trbk. nodct was the total number of nodes in the original tree.
spltct = 0; % splitting nodes count, the nodes that actually divide
div_gen_all = zeros(1,trbk); % from original tree
div_gen_mat = div_gen_all; % from rebuilt tree
inter_gen_all = div_gen_all; % from original tree
inter_gen_mat = div_gen_all; % from rebuilt tree
for i = it-1:-1:it-trbk
    ogen = otree{i}; % original generation
    splitind = splt{i}; % split node indicator
    div_gen_all(i)=sum(splitind);
    spltct = spltct + div_gen_all(i); 
    inter_gen_all(i) = size(ogen,1);
    nodctnew = nodctnew + inter_gen_all(i);
    ngen = ntree{i}; % new generation
    [m1, m2] = findm(ogen,ngen,splitind); %m1 is matched nodes, m2 is matched splitting nodes
    comptree = comptree + m1;
    splitm = splitm + m2;
    inter_gen_mat(i) = m1;
    div_gen_mat(i) = m2;
end 

% showstring = [num2str(comptree),'/',num2str(nodctnew), ',', num2str(splitm), '/', num2str(spltct)];
% disp(showstring)
% 
% disp('All internal nodes by generation')
% disp(inter_gen_all)
% disp('All matched internal nodes by generation')
% disp(inter_gen_mat)
% disp('All dividing nodes by generation')
% disp(div_gen_all)
% disp('All matched dividing nodes by generation')
% disp(div_gen_mat)


fy1=comptree/nodctnew*100;
if (spltct==0) % no splitting nodes    
    disp('No spliting nodes')
    fy2 = 0;
else
    fy2 = splitm/spltct*100;
end
