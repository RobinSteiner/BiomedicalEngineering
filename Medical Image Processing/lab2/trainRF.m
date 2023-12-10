function [rf] = trainRF(images, masks)

nTrees = 100;
features = [];
labels = [];

for i = 1:numel(images)
    image = images{i};
    mask = masks{i};
    
    % Compute features for the image
    imageFeatures = cache(@computeFeatures, image);
    
    % Get foreground indices
    foregroundIndices = find(mask);
    
    % Sample a subset of background pixels matching the number of foreground pixels
    numForegroundPixels = numel(foregroundIndices);
    backgroundIndices = find(~mask);
    backgroundIndices = backgroundIndices(randperm(numel(backgroundIndices), numForegroundPixels));
    
    % Add features and labels to the training data
    features = [features, imageFeatures(:, [foregroundIndices; backgroundIndices])];
    labels = [labels; ones(numForegroundPixels, 1); zeros(numForegroundPixels, 1)];
end

% Transpose the feature matrix
features = features';

% Train a Random Forest classifier
rf = TreeBagger(nTrees, features, labels, 'OOBVarImp', 'on');

end