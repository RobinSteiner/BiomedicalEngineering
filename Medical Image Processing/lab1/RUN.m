load daten.mat %contains data1, data2, data3, data4
load daten3d.mat %contains data (3D)
load shapes.mat %contains matrix aligned 

% Task 1 - covariance Matrix
% Computing the covariance matrices for datan.mat using ourCov.m
Covariancematrix_1 = ourCov(data1')
Covariancematrix_2 = ourCov(data2')
Covariancematrix_3 = ourCov(data3')
Covariancematrix_4 = ourCov(data4')

% Visulazing the data
dataPlots

% Task 2 - PCA
% Computing the PCA 
[eigenvalues_1,eigenvectors_1] = ourPCA(data1')
[eigenvalues_2,eigenvectors_2] = ourPCA(data2')
[eigenvalues_3,eigenvectors_3] = ourPCA(data3')
[eigenvalues_4,eigenvectors_4] = ourPCA(data4')

% plot the results 
pcaPlots

% Task 3 - Subspace projection
% Projecting data2 to main vector PC1 and to side vector PC2
% Reconstructing the projection and plotting results
subspaceprojection

% Task 4 - Investigation in 3D
% Perform the PCA & plot the data & eigenvectors of the data in daten3d.mat
% Projecting the data in the subspace constructed by the first 2
% Eigenvectors and reconstruct the points in the original space
pca3D(data)

% Task 5 - Shape modeling
% Perform PCA on shape.mat (matrix aligned with nPoints x nDimensions x nShapes)
[eigenvalues_aligned,eigenvectors_aligned] = ourPca(dimreduction(aligned));

% Implement function to compute and plot new shapes
shapeModeling


