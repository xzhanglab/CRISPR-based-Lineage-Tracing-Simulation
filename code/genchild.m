% This program generates a child barcode from a parent barcode
function [y, nstp, nendp] = genchild(n,ldtl,stp,endp,parbarc,mupb,lgdelprob,ins_sub)
% k is division round
% j is cell index in this generation
% stp =sted(1,2^(k-1)+j-1) (for first child) is the start position in the barcode
% endp = sted(2,2^(k-1)+j-1) (for first child)is the end position
% parbarc = barM(:,2^(k-1)+j-1) is the parent barcode
cutsites = []; % record all cut sites
                    for i=stp:endp % for each position of its parent barcode
                       if (parbarc(i)>0) % valid letter, i.e., undeleted position
                            mut = rand; % indicator of mutation/cut
                            if (mut<=mupb) % mutation/cut occurs
                                cutsites = [cutsites i]; %dynamically extend cutsites
                            end
                       end
                    end
y=parbarc; % full parent barcode


                    if (length(cutsites)>=2)% more than 2 cuts
                        if (rand<lgdelprob) % perform one large deletion
                            randi = randsample(length(cutsites),2); % randomly select two sites to have a large deletion
                            y(cutsites(min(randi)):cutsites(max(randi))) = 0; 
                            cutsites(min(randi):max(randi)) = 0; % clear all the involved cutsites
                        end
                    end
                    
                    newcut = cutsites(cutsites~=0);
                    for i=1:length(newcut) % process all other cut sites
                        mutp1 = newcut(i);
                        % an effective cut, determine if insert 1,2,3, or sub
                            mut = rand; % indicator of insertion amount
                            if (mut<ins_sub(1)) % perfect repair
                            elseif (mut<ins_sub(2))% insert 1, adjust/shift segment of vector
                                if (stp<n+2*ldtl-endp) % insert downstream
                                    y(mutp1+1:endp+1)=y(mutp1:endp);
                                    endp = endp+1;
                                    newcut(i:end)=newcut(i:end)+1; % update rest cut sites, only for downstream insertion
                                else % insert upstream
                                    y(stp-1:mutp1-1)=y(stp:mutp1);
                                    stp=stp-1;
                                end
                                % insert a random letter at this position
                                y(mutp1)=ceil(4*rand);
                            elseif (mut<ins_sub(3)) % insert 2
                                if (stp<n+2*ldtl-endp) % insert downstream
                                    y(mutp1+2:endp+2)=y(mutp1:endp);
                                    endp=endp+2;
                                    % insert 2 random letters downstream
                                    y(mutp1:mutp1+1)=ceil(4*rand(2,1));
                                    newcut(i:end)=newcut(i:end)+2; % update rest cut sites, only for downstream insertion
                                else % insert upstream
                                    y(stp-2:mutp1-2)=y(stp:mutp1);
                                    stp=stp-2;
                                    % insert 2 random letters uptream
                                    y(mutp1-1:mutp1)=ceil(4*rand(2,1));
                                end
                            elseif (mut<ins_sub(4)) % insert 3
                                if (stp<n+2*ldtl-endp) % insert downstream
                                    y(mutp1+3:endp+3)=y(mutp1:endp);
                                    endp=endp+3;
                                    % insert 3 random letters downstream
                                    y(mutp1:mutp1+2)=ceil(4*rand(3,1));
                                    newcut(i:end)=newcut(i:end)+3; % update rest cut sites, only for downstream insertion
                                else % insert upstream
                                    y(stp-3:mutp1-3)=y(stp:mutp1);
                                    stp=stp-3;
                                    % insert 2 random letters uptream
                                    y(mutp1-2:mutp1)=ceil(4*rand(3,1));
                                end
                            elseif (mut<ins_sub(5)) % substitution, i.e., single nucleotide mutation
                                mua = ceil(3*rand); % amount of mutation
                                y(mutp1)= y(mutp1)+mua;
                                if (y(mutp1)>4)
                                    y(mutp1)=y(mutp1)-4;
                                end
                            else % single nucleotide deletion
                                y(mutp1)=0;
                            end
                    end % end for all other cut sites
nstp = stp;
nendp = endp;