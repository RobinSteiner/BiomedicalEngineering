% Task 5.(c)

load shapes.mat

dataMatrix = dimreduction(aligned); %converting 128x2x14 matrix into 256x14 in order to be able to perform pca

[eigVal, eigVec]= ourPca(dataMatrix');

m = mean(dataMatrix')';

[nEigenvectors, ~] = size(m);
meanShape = zeros(nEigenvectors/2,2);

meanShape(:,1) = m(1:2:end);
meanShape(:,2) = m(2:2:end);

% b: Plotting with different coefficients for one of the first 3 PCs 
stddeviations = sqrt(abs(eigVal))';

plotShape(meanShape, eigVec, [-3*stddeviations(1) 0 0]);
plotShape(meanShape, eigVec, [3*stddeviations(1) 0 0]);

plotShape(meanShape, eigVec, [0 -3*stddeviations(2) 0]);
plotShape(meanShape, eigVec, [0 3*stddeviations(2) 0]);

plotShape(meanShape, eigVec, [0 0 -3*stddeviations(3)]);
plotShape(meanShape, eigVec, [0 0 3*stddeviations(3)]);

% c: generating and plotting bones with different numbers of PCs using a
% random generated coefficient matrix b
b = randn(1,nEigenvectors).*stddeviations;
plotShape(meanShape, eigVec, b);

sumStd = sum(stddeviations);
percentages = cumsum(stddeviations./sumStd*100);


% Define the thresholds
thresholds = [80, 90, 95, 100];


% Find the first elements above the thresholds and plot
for threshold = thresholds
    indices = find(percentages >= threshold);
    if ~isempty(indices)
        fprintf("Using the first %d principal components captures %d%% of the total variance of the data\n", indices(1), threshold)
        plotShape(meanShape, eigVec, b(1:indices(1)));
    end
end



