% Task 4.

function pca3D(D)

d = D';

% we first perform PCA :
[eigVal, eigVec] = ourPca(d);

% and then plot the data and the Eigenvectors - the projection
plot3DPCA(d, mean(d), eigVec, eigVal, 1, 0);

% and then do the reconstruction
plot3DPCA(d, mean(d), eigVec, eigVal, 1, 1);

end