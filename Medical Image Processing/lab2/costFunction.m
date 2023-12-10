% 2.4.b)
function [d] = costFunction(p, classificationResult)
    
    % We have to extract the two matrices that we are going to compare by
    % calculating the euclidean distance between corresponding coordinates
    % values
    
    % First matrix : Compute the shape from the parameter vector p
    b = p(1);
    scaling = p(2); 
    rotation = p(3); 
    x_translation = p(4); 
    y_translation = p(5);
    
    S = generateShape(b, scaling, rotation, x_translation, y_translation)' % dimension : 2*128
    
    % Second matrix : Resulting from classification output
    % We need to export the coordinates of the pixels of value 1
    
    [j,i,v] = find(classificationResult); % i and j coordinates of non zero values v 
    C = [i j]
    
    % Compute the cost based on the comparison
    D = pdist2(S,C,'euclidean'); % S and C must be of same dimension
    d = trace(D)
    
    
    % Example
    X = [40 88; 51 88; 35 78; 36 75; 39 72; 44 71; 48 71; 52 74; 55 77];
    Y = [36 43; 48 42; 31 26; 33 28; 37 30; 40 31; 45 30; 48 28; 51 24];

    D = pdist2(X,Y,'euclidean');
    dist = trace(D)

    % if X and Y are very similar :
    X = [40 88; 51 88; 35 78; 36 75; 39 72; 44 71; 48 71; 52 74; 55 77];
    Y = [43 88; 52 85; 35 78; 36 75; 39 76; 34 73; 48 71; 52 74; 55 77];

    D = pdist2(X,Y,'euclidean')
    dist = trace(D)

    % d diminishes 
    
end



% 2.4.c) Cost function optimization
function optimizedParams = optimization(d, mi, ma, drawFunction)
    
    optimizedParams = optimize(d, mi, ma, drawFunction); 
    
end
