% Simulation of CRISPR-Cas9 Barcode editting
clear;
clc;
%================parameter setting========================================
n = 100; % initial barcode length
it = 5; % cell division rounds
propm = 0.7; % proportion of nonzero counts of a barcode to be considered as a matched pair
ss = 1; % sample size proportion, a number between 0 and 1.
mupb = 0.1; % probability of a cut at each position

ins_sub = [0.7 0.8 0.83 0.85 0.9]; % probabilities of perfect repair (0,0.7), inserting 1 (0.7,0.8), 2 (0.8,0.83), 3 (0.83,0.85) nucleotides, substitution (0.85,0.9), or single nucleotide deletion (0.9,1), respectively. The actual probability equals the length of corresponding interval.
lgdelprob = 0.15; % this is the probability to have a large deletion, given more than 2 cut sites.

pulse = 0; % if pulse = 1, pulse induction; if pulse = 0, constant dox level

divp = 1; % probability of division
        %if division occurs, get two children, then check probability of
        %survival, then mutate.
        %if not dividing, get one child, check probability of surviving,
        %then simulate mutation.
clive = 1; % probability of cell survival/live. Recommended value is 1. NOTE: program may generate error message if all leaf barcodes are empty (all dead cells).


trbk = it-1; % this is the number of generations to trace back, trbk<it, trbk is at most it-1.

%==================1 barcode==============================================
        %[accuracyRMP(1),accuracyRMP(2)]=funbarnew(n,it,propm,ss,mupb,ins_sub,lgdelprob,divp,clive, pulse, trbk); % RMP method
        % [accuracyRMPNF(1), accuracyRMPNF(2)]=funbarnewRMPNF(n,it,propm,ss,mupb,ins_sub,lgdelprob,divp,clive, pulse, trbk); % RMPNF method
        % [accuracyNBJ(1), accuracyNBJ(2)]=funbarNBJ(n,it,propm,ss,mupb,ins_sub,lgdelprob,divp,clive, pulse, trbk); % NBJ method
         [accuracyNBJNF(1), accuracyNBJNF(2)]=funbarNBJNF(n,it,propm,ss,mupb,ins_sub,lgdelprob,divp,clive, pulse, trbk); % NBJNF method
 
    disp('Accuracy of all internal nodes');
    %accuracyRMP(1)
    % accuracyRMPNF(1)
    % accuracyNBJ(1)
     accuracyNBJNF(1)

    disp('Accuracy of dividing nodes');
    %accuracyRMP(2)
    % accuracyRMPNF(2)
    % accuracyNBJ(2)
     accuracyNBJNF(2)

% =========================2 barcodes===================================
propmi = 0.4; % for individual barcode

 % [accuracyNBJ(1), accuracyNBJ(2)]=funbarNBJ2(n,it,propm,propmi,ss,mupb,ins_sub,lgdelprob,divp,clive, pulse, trbk); % NBJ method
 % [accuracyNBJNF(1), accuracyNBJNF(2)]=funbarNBJ2trbk(n,it,propm,propmi,ss,mupb,ins_sub,lgdelprob,divp,clive, pulse, trbk); % NBJNF method
 %  disp('Accuracy of all internal nodes');
 %  accuracyNBJ(1)
 %  accuracyNBJNF(1)
 %  disp('Accuracy of dividing nodes');
 %  accuracyNBJ(2)
 %  accuracyNBJNF(2)
