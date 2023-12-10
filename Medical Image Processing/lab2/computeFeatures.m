function [M] = computeFeatures(image)

[pixel_x,pixel_y]=size(image);                  

feat = 8;                                       %number of features for the exercise

M = zeros(feat,pixel_x,pixel_y);                %initialisation of feature Matrix #feature x size of image

%Feature 1 - grey values
M(1,:,:) = image;           

%Gradient in x- and y- direction and magnitude
[gradient_x,gradient_y] = imgradientxy(image);
[gradient_mag,~] = imgradient(gradient_x,gradient_y); 

%Feature 2 - Gradient in x-direction
M(2,:,:) = gradient_x;

%Feature 3 - Gradient in y-direction
M(3,:,:) = gradient_y;

%Feature 4 - Magnitude of Gradient
M(4,:,:) = gradient_mag;

%Haar-like features of the image with given function
%Feature 5 - Haar-like for gray values
%Feature 6 - Haar-like for gradient-magnitude

HLF_GreyVal = computeHaarLike(image);

HLF_gradient_mag = computeHaarLike(gradient_mag); %gives Matrix 20x43758 where rows stand for elements of the filter

%only use the first row and reshape into matrix 306x143
%for gray values
HLF_GreyVal_first = HLF_GreyVal(1,:);
M(5,:,:) = reshape(HLF_GreyVal_first,pixel_x,[]);

%for gradient-magnitude
HLF_gradient_mag_first = HLF_gradient_mag(1,:);
M(6,:,:) = reshape(HLF_gradient_mag_first,pixel_x,[]);

%Feature 7 - Indices of grey-values in x direction
for i=1:pixel_x
    M(7,i,:)=i;
end

%Feature 8 - Indices of grey values in y direction
for j=1:pixel_y
    M(8,:,j)=j;
end
end