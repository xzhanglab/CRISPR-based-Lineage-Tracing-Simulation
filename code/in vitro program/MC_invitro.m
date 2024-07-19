% run read tree files with different methods, Monte Carlo
clear;
clc;

MCit = 50; % monte carlo run round

Taccu = zeros(76,9);

for fnumb = 1:76 
    fnumb
    Taccu(fnumb,1:2) = readtree(MCit, 0.7,fnumb);
    Taccu(fnumb,3:4) =readtreeRMP_NF(MCit, 0.7,fnumb);
    Taccu(fnumb,5:6) =readtreeNBJ(MCit, 0.3,fnumb);
    Taccu(fnumb,7:9) =readtreeNBJ_NF(MCit, 0.3,fnumb);
end
Taccu