% Task 2.(b) 
% Plotting for the 2D data set

data1Transposed = data1';
[eVal, eVec] = ourPca(data1Transposed);
plot2DPCA(data1Transposed, mean(data1Transposed), data1Transposed, eVec, eVal, 1, 0);

data2Transposed = data2';
[eVal, eVec] = ourPca(data2Transposed);
plot2DPCA(data2Transposed, mean(data2Transposed), data2Transposed, eVec, eVal, 1, 0);

data3Transposed = data3';
[eVal, eVec] = ourPca(data3Transposed);
plot2DPCA(data3Transposed, mean(data3Transposed), data3Transposed, eVec, eVal, 1, 0);

data4Transposed = data4';
[eVal, eVec] = ourPca(data4Transposed);
plot2DPCA(data4Transposed, mean(data4Transposed), data4Transposed, eVec, eVal, 1, 0);

