% This program simulates barcodes evolution with irregular cell division
% and cell death, and with varying crispr-cas9 speed.
function [fy1,fy2] = funbarnewRMPNF(n,it,propm,ss,mupb,ins_sub,lgdelprob,divp,clive, pulse, trbk)

% RMP method, nonfiltered

ldtl = ceil(1*n); % lead and tail space on the ends of barcode
% assign a large negative number to indicate empty barcodes
empbar = -4*(n+2*ldtl);

barM = zeros(n+2*ldtl,2^(it+1)-1); % barcode matrix
barM(ldtl+1:ldtl+n,1)=ceil(4*rand(n,1)); % initial barcode A,C,G,T = 1,2,3,4
sted = ones(2,2^(it+1)-1); % record start and end position of each barcode
                            % mostly due to insertion
sted(1,:)=(ldtl+1);
sted(2,:) = ldtl+n;

inmupb = mupb;

% start mutation, construct the evolution tree-----------------------------
for k=1:it    % cell division iteration = generation number
    simtree{k}=zeros(2^k,2); % first column is cell index in this generation, second column is cell generation number
    simtree{k}(:,1)=1:2^k;  % generation number=0 means dead cell

    % test pulse dox concentration
	if (pulse == 1)
     		if (mod(k,2)==0)
         		mupb = 0.005; % base cut rate, suppose no dox
     		else
         		mupb = inmupb; % with pulses of dox
     		end
	end
    
    if (k==1) % no need to check parent's life, assume the root is a live cell
            for j=1:2^(k-1) 
                % first child---------------- 
                    rliv = rand; %probability of a live child
                    if (rliv<clive) % live child, generate it
                        [barM(:,2^k+2*(j-1)), sted(1,2^k+2*(j-1)), sted(2,2^k+2*(j-1))]=genchild(n,ldtl,sted(1,2^(k-1)+j-1),sted(2,2^(k-1)+j-1),barM(:,2^(k-1)+j-1),mupb,lgdelprob,ins_sub);
                        if (nnz(barM(:,2^k+2*(j-1)))>0) % live cell
                            simtree{k}(2*j-1,2)=k;
                        end
                    end
                    % 2nd child-------------------------------------------    
                    rliv = rand; %probability of a live child
                    if (rliv<clive) % live child, generate it
                        [barM(:,2^k+2*j-1), sted(1,2^k+2*j-1), sted(2,2^k+2*j-1)]=genchild(n,ldtl,sted(1,2^(k-1)+j-1),sted(2,2^(k-1)+j-1),barM(:,2^(k-1)+j-1),mupb,lgdelprob,ins_sub);
                        if (nnz(barM(:,2^k+2*j-1))>0) % live cell
                            simtree{k}(2*j,2)=k;
                        end
                    end
            end % j for this generation
    else %k>1, need to check parent's life
        for j=1:2^(k-1) 
            if (simtree{k-1}(j,2)>0) % for each live cell in parent generation
                rdiv = rand; % probability of division
                if (rdiv<divp) % division occurs
                    % first child---------------- 
                    rliv = rand; %probability of a live child
                    if (rliv<clive) % live child, generate it
                        [barM(:,2^k+2*(j-1)), sted(1,2^k+2*(j-1)), sted(2,2^k+2*(j-1))]=genchild(n,ldtl,sted(1,2^(k-1)+j-1),sted(2,2^(k-1)+j-1),barM(:,2^(k-1)+j-1),mupb,lgdelprob,ins_sub);
                        if (nnz(barM(:,2^k+2*(j-1)))>0) % live cell
                            simtree{k}(2*j-1,2)=k;
                        end
                    end
                    % 2nd child-------------------------------------------    
                    rliv = rand; %probability of a live child
                    if (rliv<clive) % live child, generate it
                        [barM(:,2^k+2*j-1), sted(1,2^k+2*j-1), sted(2,2^k+2*j-1)]=genchild(n,ldtl,sted(1,2^(k-1)+j-1),sted(2,2^(k-1)+j-1),barM(:,2^(k-1)+j-1),mupb,lgdelprob,ins_sub);
                        if (nnz(barM(:,2^k+2*j-1))>0) % live cell
                            simtree{k}(2*j,2)=k;
                        end
                    end
                else % single child
                     % only child---------------- 
                    rliv = rand; %probability of a live child
                    if (rliv<clive) % live child, generate it
                        [barM(:,2^k+2*(j-1)), sted(1,2^k+2*(j-1)), sted(2,2^k+2*(j-1))]=genchild(n,ldtl,sted(1,2^(k-1)+j-1),sted(2,2^(k-1)+j-1),barM(:,2^(k-1)+j-1),mupb,lgdelprob,ins_sub);
                        if (nnz(barM(:,2^k+2*(j-1)))>0) % live cell
                            simtree{k}(2*j-1,2)=k;
                        end
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




%reverse alignment between root and leaves, consecutive matches preferred-----------------------
revtrac = zeros(n+2*ldtl, 2^it+1); % reverse trace
revtrac(:,1)=barM(:,1); % show the root/first cell barcode
simiM2 = zeros(2^it,1);
tvec1 = barM(1+ldtl:n+ldtl,1);
colv1 = tvec1(tvec1~=0);
for i=1:2^it
    if(simtree{it}(i,2)>0) % this is a live cell
        % collapse columns        
        tvec2 = barM(sted(1,2^it+i-1):sted(2,2^it+i-1),2^it+i-1);        
        colv2 = tvec2(tvec2~=0);
        
        if (isempty(colv2))% 
            %warning('empty barcode');
            %i
            simiM2(i) = empbar;
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
% % draw reconstructed leaves
% figure
% heatmap(revtrac);
% title('Reverse alignment consecutive matches preferred')


%clear tvec1 colv1;

% figure
% plot(simiM2,'*');
% title('Matching scores: Root vs leaves');

% compare similarity between cells in the last generation, leaves------------
simiM = zeros(2^it,2^it); 
for i=1:2^it
    if (simtree{it}(i,2)>0) % live cell
        for j=i+1:2^it
            if (simtree{it}(j,2)>0) % live cell
                tvec1 = revtrac(:,i+1);
                tvec2 = revtrac(:,j+1);
                % trim leading and tail 0s--------
                    colv1 = trim(tvec1,2*ldtl+n);
                        %----------------                
                    colv2 = trim(tvec2,2*ldtl+n);

                % compare
                simiM(i,j)=vcomp5(colv1,colv2);
                %clear tvec1 tvec2 colv1 colv2;
            end
        end
    end
end


% make a symmetric matrix
% simiM=simiM+simiM'-eye(size(simiM)).*diag(simiM);
simiM=simiM+simiM';
for i=1:2^it
    simiM(i,i)=empbar;
end
% figure
% heatmap(simiM);
% title('Similarity among leaves');


% reconstruct lineage tree-------------------------------
celltag = genN(simiM, simiM2, it,ss,simtree{it}); % need to determine each barcodes's level/generation........., genN sorted the cells according to similarity to root
indcell = celltag(:,1); % contains all live cells
celltag(:,2)=it;% this is for testing, assuming all leaves have the same level
%--------------------------------------------------------------------
valflag = valcheck(celltag);
if (valflag == 0) % not valid assignment
    error('Not valid assignment of generation numbers');
end

%---------------------build the original tree-----------
[otree, splt, nodct] = btree(celltag); % original tree and number of nodes
%------------------------------------------------------

maxgen = max(celltag(:,2)); % current max generation number
m = length(indcell);

% record the barcodes of leaves
tbarC = zeros(n+2*ldtl,m);
for i = 1:m
    tbarC(:,i) = revtrac(:,indcell(i)+1); % record the sampled (not necessarily all) barcodes of the 'last' generation, i.e., available information at the beginning
                            % reorder the cells in tbarC
end

% figure
% heatmap(tbarC);

%lineage = zeros(length(indcell),2*(maxgen-1)); % both columns are nonempty=paired, left column nonempty = singleton, right column nonempty = unchanged.
C={}; % use cell array to store lineage
indice={}; % record cell orders
while maxgen > 1 % start rebuilding the binomial tree---------------------
    
%     figure
%     heatmap(tbarC);

    indice{maxgen-1}=indcell;

    m=length(indcell);

    
    if (maxgen<it)
        nsimiM = zeros(m,m); % new similarity between cells
    end
    
    paired = [];
    lv = length(indcell); % length of vector
    ncelltag = []; % new cell tag
    ncount = 0; % new cell count
    singleton = [];
    unchanged = []; % unchanged barcodes due to lower generation
    for i = 1:lv % look for cells in the last generation
        bmatch = empbar; % current best match
        curpar = 0; % current pair index
        if (celltag(i,2)==maxgen)
            %tvec1 = trim(tbarC(:,indcell(i)),2*ldtl+n);
            tvec1 = trim(tbarC(:,i),2*ldtl+n);
            colv1 = tvec1(tvec1~=0);
            if (~isempty(colv1)) % found one valid cell, then look for its pair   
                %-----------search for its pair
                for j = i+1:lv
                    if (celltag(j,2)==maxgen) % same generation
                        %tvec2 = trim(tbarC(:,indcell(j)),2*ldtl+n);
                        tvec2 = trim(tbarC(:,j),2*ldtl+n);
                        colv2 = tvec2(tvec2~=0);
                        if (~isempty(colv2)) % nonempty
                            
                            if (maxgen==it)% simiM was calculated before, so make computing faster
                                tscore = simiM(indcell(i),indcell(j)); % but it uses original index
                            else
                                tscore = vcomp5(colv1,colv2);
                                nsimiM(i,j)=tscore;
                                nsimiM(j,i)=tscore;
                            end
                            
                            if tscore > bmatch
                                bmatch = tscore;
                                curpar = j; %current pair
                            end
                        end
                    end
                end
                %--------------------------
                if (bmatch>propm*nnz(colv1)) % valid pair, this threshold depends on the score matrix in vcomp and vcomp5
                    paired = [paired; i curpar];
                    ncount = ncount + 1;
                    ncelltag = [ncelltag; ncount maxgen-1];
                    celltag(i,2)=0; % remove cell i
                    celltag(curpar,2)=0; % remove cell curpar
                else % cannot find a pair for i
                     singleton = [singleton; i];
                end
                
            end
        elseif (celltag(i,2)>0) % generation < maxgen
%            unchanged = [unchanged; celltag(i,:)]; % here the index is original in tbarC
            unchanged = [unchanged; i celltag(i,2)]; % here the index is new
            warning('you have unchanged barcodes due to different generation number');
        else % already deleted cells
        end
    end
    
    % now ncount is the length of paired
    ctunch = size(unchanged,1); % length of unchanged

    %check validity of the new generation
    lv = length(singleton);
    if (lv>0)
        tpars = [ncelltag; unchanged; zeros(lv,2)];
        tpars(ncount+ctunch+1:end,2) = maxgen-1;
        while (~valcheck(tpars) && lv>1)% cell numbers in each generation do not pass validity check, make one pair
            % look for the best pair
%             curpar = 0;           
            bmatch = empbar;
            curpar1 = 0;
            curpar2 = 0;
            if (maxgen == it)
                for i=1:lv-1
                    for j = i+1:lv
                        if (simiM(indcell(singleton(i)),indcell(singleton(j)))>bmatch)
                            bmatch=simiM(indcell(singleton(i)),indcell(singleton(j)));
                            curpar1 = i;
                            curpar2 = j;
                        end
                    end
                end
            else
                for i=1:lv-1
                    for j = i+1:lv
                        if (nsimiM(singleton(i),singleton(j))>bmatch)
                            bmatch=nsimiM(singleton(i),singleton(j));
                            curpar1 = i;
                            curpar2 = j;
                        end
                    end
                end
            end
%             %tvec1 = trim(tbarC(:,indcell(singleton(1))),2*ldtl+n);
%             tvec1 = trim(tbarC(:,singleton(1)),2*ldtl+n);
%             colv1 = tvec1(tvec1~=0);
%             for j = 2:lv
%                 %tvec2 = trim(tbarC(:,indcell(singleton(j))),2*ldtl+n);
%                 tvec2 = trim(tbarC(:,singleton(j)),2*ldtl+n);
%                 colv2 = tvec2(tvec2~=0);
%                 tscore = vcomp5(colv1,colv2);
%                 if tscore > bmatch
%                    bmatch = tscore;
%                    curpar = j;
%                 end
%             end
            % now curpar1 and curpar2 are the pair
            ncount = ncount + 1;
            %paired = [paired; singleton(1) singleton(curpar)];
            paired = [paired; singleton(curpar1) singleton(curpar2)];
            ncelltag = [ncelltag; ncount maxgen-1];
            remc = 0; % remaing cell count
            newsing = zeros(lv-2,1); % new singleton cell list
%             for j=2:lv
%                 if (j~=curpar)
%                     remc = remc+1;
%                     newsing(remc) = singleton(j);
%                 end
%             end
            for j=1:lv
                if (j~=curpar1 && j~=curpar2)
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
    
    % need to rebuild parent nodes and compare to root
    thisgen = zeros(2*ldtl+n,ncount+lv); % barcodes in this iteration
    % it consists of three parts: paired, unchanged yet, and singleton
    tsize = size(paired,1);
    curlin = zeros(ncount+lv,2); % current lineage
    for i = 1:tsize % paired--------------------part I
        % extract lineage information
%         lineage(i,2*maxgen-3)=indcell(paired(i,1));
%         lineage(i,2*maxgen-2)=indcell(paired(i,2));
%         curlin(i,1)=indcell(paired(i,1));
%         curlin(i,2)=indcell(paired(i,2));

        curlin(i,1)=paired(i,1);
        curlin(i,2)=paired(i,2);

%         colv1 = trim(tbarC(:,indcell(paired(i,1))),2*ldtl+n);
%         colv2 = trim(tbarC(:,indcell(paired(i,2))),2*ldtl+n);
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
%             colv1 = trim(parent,2*ldtl+n); 
%             % now colv1 is the temporary parent
%             [~, outv]=vcomp5(barM(ldtl+1:ldtl+n,1),colv1); % compare to root

            % try a different way of alignment-collapse and redo vcomp_cns
            colv1 = parent(parent~=0); 
            % now colv1 is the temporary parent
            [~, outv]=vcomp_cons(barM(ldtl+1:ldtl+n,1),colv1); % compare to root
            
            tlength = size(outv,1);
            parent = outv(:,2);
%             for j = 1:tlength % reconstruct parent node
%                 if (outv(j,1)~=parent(j)) % with a certain probability to replace entry with the root entry
%                     if(rand<1/maxgen)
%                         parent(j)=outv(j,1);
%                     end
%                 end
%             end
            %thisgen(:,i) =  [zeros(ldtl,1);parent;zeros(n+ldtl-tlength(1),1)];
            parent = trim(parent,2*ldtl+n);
            tlength = length(parent);
            thisgen(:,i) =  [parent;zeros(n+2*ldtl-tlength,1)]; % fill up to make equal length
    end % end paired part I---------------------------------------------
    
    % unchanged barcodes part II--------------------------------------
    for i=tsize+1:ncount
        % extract information, this cell does not belong to this
        % generation, leave it empty
%         lineage(i,2*maxgen-3)=0;
%         lineage(i,2*maxgen-2)=unchanged(i-tsize(1),1);
        curlin(i,2)=unchanged(i-tsize,1);
        
        parent = trim(tbarC(:,unchanged(i-tsize,1)),2*ldtl+n);
        tlength = size(parent,1);
        %thisgen(:,i) =  [zeros(ldtl,1);parent;zeros(n+ldtl-tlength(1),1)];
        thisgen(:,i) =  [parent;zeros(n+2*ldtl-tlength,1)];
    end
    
    % singleton barcodes part III-------------------------------------     
    if (lv>0) % there are remaining singleton cells
        for i = 1:lv
            ncount = ncount + 1;
%             lineage(ncount,2*maxgen-3)=indcell(singleton(i));% extract information
%            curlin(ncount,1)=indcell(singleton(i));% extract information
            curlin(ncount,1)=singleton(i);% extract information
            ncelltag = [ncelltag; ncount maxgen-1];
%            colv1 = trim(tbarC(:,indcell(singleton(i))),2*ldtl+n);

            colv1 = trim(tbarC(:,singleton(i)),2*ldtl+n);            
            colv2 = colv1(colv1~=0);
%             [~, outv]=vcomp5(barM(ldtl+1:ldtl+n,1),colv1); % compare to root

            [~, outv]=vcomp_cons(barM(ldtl+1:ldtl+n,1),colv2); % compare to root
            
            tlength = size(outv,1);
            parent = outv(:,2);
%             for j = 1:tlength % reconstruct parent node
%                 if (outv(j,1)~=parent(j)) % with a certain probability to replace entry with the root entry
%                     if(rand<1/maxgen)
%                         parent(j)=outv(j,1);
%                     end
%                 end
%             end
            %thisgen(:,ncount) = [zeros(ldtl,1);parent;zeros(n+ldtl-tlength(1),1)];
            thisgen(:,ncount) = [parent;zeros(n+2*ldtl-tlength,1)];
        end        
    end
      

    simiM2 = zeros(ncount,1);
%     tvec1 = barM(1+ldtl:n+ldtl,1); % root
%     colv1 = tvec1(tvec1~=0);

     colv1 = barM(1+ldtl:n+ldtl,1); % root
    for i=1:ncount % compare to root
        tvec2 = thisgen(:,i);        
        colv2 = trim(tvec2,2*ldtl+n); % since we know some structure about this barcode, use vcomp5
        simiM2(i)=vcomp5(colv1,colv2);
    end
    [~, indcell] = sort(simiM2,'descend'); % sorted cell according to matching scores to the root
    celltag = zeros(ncount,2); % update cell tags------------------
    celltag(:,1) = indcell;
    tbarC = zeros(n+2*ldtl,ncount);
    for i = 1:ncount
        celltag(i,2)=ncelltag(indcell(i),2);
        tbarC(:,i) = thisgen(:,indcell(i));
    end
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
%             if (C{i}(j,1)>0 && C{i}(j,2)>0) % pairs
%                 if (mod(indice{i}(C{i}(j,1)),2)==1 && indice{i}(C{i}(j,2))==indice{i}(C{i}(j,1))+1)
%                         comptree = comptree + 1; % the original tree is a full binomial tree
%                 end
%                 if (mod(indice{i}(C{i}(j,1)),2)==0 && indice{i}(C{i}(j,2))==indice{i}(C{i}(j,1))-1)
%                         comptree = comptree + 1;
%                 end
%             end
%             
%             if (C{i}(j,1)==0) % unchanged cells----------------------
%                 
%             end
%             
%             if (C{i}(j,2)==0) % singleton cells
%                 comptree = comptree + 1; % suppose a correct singleton node is found
%                 for k=1:clen
%                     if (C{i}(k,2)==0 && k~=j)
%                         if (mod(indice{i}(C{i}(j,1)),2)==1 && indice{i}(C{i}(k,1))==indice{i}(C{i}(j,1))+1)
%                             comptree = comptree - 1; % unmatched pair that should have been paired
%                         end
%                         if (mod(indice{i}(C{i}(j,1)),2)==0 && indice{i}(C{i}(k,1))==indice{i}(C{i}(j,1))-1)
%                             comptree = comptree - 1; % unmatched pair that should have been paired
%                         end
%                     end
%                 end
%             end
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
            
%              if (mod(newC(j,1),2^(it-i))==1 && newC(j,2^(it-i)) == newC(j,1)+2^(it-i)-1) % if form a full segment
%                     comptree = comptree + 1;
%              end
                
        end
        tempvec = newC;
        ntree{i} = newC;
    end
end

% calculate similarity score between two trees----------------------------
comptree = 0; % number of correct nodes, i.e., correct parent with correct decendants/clade
splitm = 0; % splitting nodes match, 
% since otree is full, and ntree may not be full, may skip the following
% comparison

% if size(otree,2)~=size(ntree,2) 
%     error('These two trees do not have the same level.');
% end
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
% scale it
%comptree/(2^it-2)*100
%fy=comptree/nodct*100;


fy1=comptree/nodctnew*100;
if (spltct==0) % no splitting nodes    
    disp('No spliting nodes')
    fy2 = 0;
else
    fy2 = splitm/spltct*100;
end