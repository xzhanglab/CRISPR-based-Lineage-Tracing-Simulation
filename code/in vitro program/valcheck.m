% This function checks validity of the cell numbers in each generation
function y = valcheck(v)
maxgen = max(v(:,2)); % record the deepest level
celln = zeros(maxgen,1);

lv = length(v(:,1));
for i = 1:lv
    celln(v(i,2))=celln(v(i,2))+1;
end

tempsum = 0;
for i = 1:maxgen
    tempsum = tempsum + celln(i)*2^(maxgen-i);
end

if tempsum>2^maxgen
    y=0; % not valid assignments of the level numbers
else
    y=1; % valid
end