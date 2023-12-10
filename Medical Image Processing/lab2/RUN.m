% Task 1
% ---------------------------------------------------------------
shape = generateShape([1 1 1], 1, 0, 0, 0);
shape1 = generateShape([1 1 1], 2, 45, 0, 0);
shape2 = generateShape([1 1 1], 1, 90, -20, 0);
shape3 = generateShape([1 1 1], 1.5, -45, 30, 0);
shape4 = generateShape([1 1 1], 2, 30, 20, -300);
shape5 = generateShape([1 1 1], 0.5, 30, 20, -300);
figure("Name", "Shap Transformation",'NumberTitle','off')
plot(shape(:,1), shape(:,2),'o')
hold on
plot(shape1(:,1), shape1(:,2),'o')
hold on
plot(shape2(:,1), shape2(:,2),'o')
hold on
plot(shape3(:,1), shape3(:,2),'o')
hold on
plot(shape4(:,1), shape4(:,2),'o')
hold on
plot(shape5(:,1), shape5(:,2),'o')
axis equal
% ---------------------------------------------------------------


% Task 2
% ---------------------------------------------------------------
load handdata.mat images masks aligned
features = cache(@computeFeatures,images{1,2});

figure("Name", "Feature 1 - Grey values",'NumberTitle','off')
imagesc(squeeze(features(1,:,:)))

axis image

figure("Name", "Feature 2 - Gradient in x-direction",'NumberTitle','off')
imagesc(squeeze(features(2,:,:)))

axis image

figure("Name", "Feature 3 - Gradient in y-direction",'NumberTitle','off')
imagesc(squeeze(features(3,:,:)))

axis image

figure("Name", "Feature 4 - Magnitude of Gradient",'NumberTitle','off')
imagesc(squeeze(features(4,:,:)));

axis image

figure("Name", "Feature 5 - Haar-like for gray values",'NumberTitle','off')
imagesc(squeeze(features(5,:,:)));

axis image

figure("Name", "Feature 6 - Haar-like for gradient-magnitude",'NumberTitle','off')
imagesc(squeeze(features(6,:,:)));

axis image

figure("Name", "Feature 7 - Indices of grey-values in x direction",'NumberTitle','off')
imagesc(squeeze(features(7,:,:)));

axis image

figure("Name", "Feature 8 - Indices of grey values in y direction",'NumberTitle','off')
imagesc(squeeze(features(8,:,:)));

axis image
% ---------------------------------------------------------------



% Task 3
% ---------------------------------------------------------------
rf = trainRF(images, masks);
figure("Name", "Permuted Var Delta Error",'NumberTitle','off')
bar(rf.OOBPermutedVarDeltaError, 0.4);

figure("Name", "Out-of-bag error",'NumberTitle','off')
plot(oobError(rf));
% ---------------------------------------------------------------



% Task 4
% ---------------------------------------------------------------
[rf, meanShape, eVecSorted, eValSorted] = train(images, masks, aligned);

% ---------------------------------------------------------------



