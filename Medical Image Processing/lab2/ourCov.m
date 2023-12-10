% Task 1.(a)

function[covarianceMatrix] = ourCov(data)

[m,n] = size(data); % m rows, n columns 
covarianceMatrix = zeros(n); % nxn Matrix

% covariance matrix calculation
for i = 1:n
    for j = 1:n
        covarianceMatrix(i,j) = covariance(data(:,i), data(:,j));
    end
end
end

% function calculating the covariance of 2 vectors A and B
function c = covariance(A, B)
l = length(A);
s = 0;
for k = 1:l
    s = s + (A(k) - mean(A)) * ((B(k)) - mean(B));
end

c = s/(l-1);

end