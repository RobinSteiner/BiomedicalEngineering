% Function used for task 5

function [M] = dimreduction(A)

[m,n,l]=size(A); %m=128, n=2, l=14 -> 14 times 128x2 matrix ( in the order [x1, y1, x2, y2, ...] )
M = zeros(m*n,l); % dim(M)=256x14

for i=1:l
    for k = 1:m
        M(k*2-1,i) = A(k,1,i) ;
        M(k*2,i) = A(k,2,i);
    end
end

end