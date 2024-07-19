%This function trim the leading and tail zeros of a vector
function y = trim(v,totl) % input: vector, and total possible length
sti = 1; % start index

while ((v(sti)==0) && (sti<totl))
    sti = sti+1;
end

edi = length(v); % end index
while((v(edi)==0)&&(edi>1))
 edi = edi-1;
end

if (sti<=edi)
y = v(sti:edi);
else
%warning('Null cell');
y = 0;
end