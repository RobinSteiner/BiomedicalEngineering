% Task 5.(a)

function [newShape] = generateShape(b, scaling, rotation, x_translation, y_translation) %output is new shape in matrix in form 128x2
load shapes.mat aligned

dataMatrix = dimreduction(aligned); %converting 128x2x14 matrix into 256x14

[~, eigVec]= ourPca(dataMatrix');
     
newShape_onecolumn = mean(dataMatrix')' + eigVec(:,1:length(b))*b';

%convert back into 128x2 matrix with x & y coordinates in the columns
[rows, ~] = size(newShape_onecolumn);
newShape = zeros(rows/2,2);

newShape(:,1) = newShape_onecolumn(1:2:end);
newShape(:,2) = newShape_onecolumn(2:2:end);

R = [cosd(rotation) -sind(rotation); 
    sind(rotation) cosd(rotation)];
S = [scaling 0; 
    0 scaling];

newShape = newShape*R*S;

newShape(:,1) = newShape(:,1)+y_translation;
newShape(:,2) = newShape(:,2)+x_translation;
end
