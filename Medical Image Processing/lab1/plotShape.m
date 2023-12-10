% Task 5.(b) 

function plotShape(meanShape, eigVectors, b, x)
figure
shape = generateShape(b);
plot(shape(:,1), shape(:,2),'o')
hold on
plot(meanShape(:,1), meanShape(:,2),'o') %to add on the mean shape 
axis equal

end