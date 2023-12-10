% Task 3.

plotProjectionAndCalculateError(data2', "Data 2");

%plotProjectionAndCalculateError(data1', "Data 1");
%plotProjectionAndCalculateError(data3', "Data 3");
%plotProjectionAndCalculateError(data4', "Data 4");

function plotProjectionAndCalculateError(data, title)
[eVal, eVec] = ourPca(data);

% transformation into PC base
dataInPCSpace = data * eVec;

% projecting onto the first and second PC by simply omitting the other axis
% and shifting to the mean along the omitted axis
meanPC = mean(dataInPCSpace, 1);
dataPC1Projected = [dataInPCSpace(:,1), zeros(50,1) + meanPC(2)];
dataPC2Projected = [zeros(50,1) + meanPC(1), dataInPCSpace(:,2)];

% calculating the average error of the projection onto PC1 and PC2
[n,~] = size(data); % m rows, n columns
errorPC1Projection = sum(abs(dataInPCSpace(:,2) - meanPC(2))) / n;
errorPC2Projection = sum(abs(dataInPCSpace(:,1) - meanPC(1))) / n;
fprintf('mean error of %s projection onto PC1: %f\n', title, errorPC1Projection)
fprintf('mean error of %s projection onto PC2: %f\n', title, errorPC2Projection)

% reconstructing into the inital data space
dataPC1Reconstructed = dataPC1Projected * eVec';
dataPC2Reconstructed = dataPC2Projected * eVec';

% plotting for both PC1 and PC2
plot2DPCA(data, mean(data), dataPC1Reconstructed, eVec, eVal, 1, 1);
plot2DPCA(data, mean(data), dataPC2Reconstructed, eVec, eVal, 1, 1)
end