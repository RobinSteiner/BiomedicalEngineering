% Task 2.(a)

function[eValSorted, eVecSorted] = ourPca(data)
[eVec, eVal] = eig(ourCov(data));

% sorting the eigenvectors by the eigenvalues descending
[eValSorted, ind] = sort(diag(eVal), 'descend');
eVecSorted = eVec(:, ind);

end
