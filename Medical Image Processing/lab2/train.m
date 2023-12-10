% Exercise 2.4: Shape Particle Filters

% 2.4.a)

function [rf, meanShape, eVecSorted, eValSorted] = train(images, masks, aligned)
   
    % Training classifier based on random forest 
    rf = trainRF(images, masks);

    
    % This second part is used in order to generate a PCA shape model
    % based on the aligned dataset

    [PointsNb, ~, samples] = size(aligned(:,:,1:30));
    
    shapes2D = reshape(aligned(:,:,1:30), [2 * PointsNb, samples] );

    meanShape = sum(shapes2D,2) / samples;    

    [eValSorted, eVecSorted] = ourPca(shapes2D.');
    eValSorted(samples + 1: end,:) = [];
    eVecSorted(:,samples+1:end) = [];
end




function MaskPred = predictSegmentation(rf, testImage)
    features = cache(@computeFeatures, testImage);
    prediction = str2double(predict(rf, features.'));
    MaskPred = reshape(prediction, size(testImage));
end





