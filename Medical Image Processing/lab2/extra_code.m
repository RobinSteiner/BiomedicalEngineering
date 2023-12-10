% Image and corresponding segmentation

n = randi([1,50])

M = masks{1,n};
[j,i,v] = find(M); % i and j coordinates of non zero values v 

C = [i j];

imshow(images{1,n})
hold on
plot(i, j,'.')