% Assignment 04 template script
%% Do not forget to rename the script before submission *******************
%% Section S1
preamble
load struct04_Sec1.mat S % load structure array S from file struct04_Sec1.mat
% *************************************************************************
% !!!!! struct04_Sec1.mat must be in the same directory as the script !!!!!
% *************************************************************************
result=find_match(S,'016')
% correct answer: 1     5     9    15
result=find_match(S,033282)
% correct answer: warning message and empty result
result=find_match(S,'Moha')
% correct answer: 3     24
result=find_match(S,'moha')
% correct answer: 3     24
result=find_match(S,'0663')
% correct answer: 2     8    14    16    20    23    24    25    30
result=find_match(S(1:20),'0663')
% correct answer: 2     8    14    16    20
result=find_match(S,'0663','Scode')
% correct answer: warning message and empty result
result=find_match(S,'sal','name')
% correct answer: 16    19
result=find_match(S,'016','matnum')
% correct answer: 1     5     9    15
result=find_match(S,'w','gender')
% correct answer: 1     8    15    22    27
result=find_match(S,'a','gender')
% correct answer: empty result
result=find_match(S,'y')
% correct answer: 3    17    18    20    30
result=find_match('S','y')
% correct answer: warning message and empty result
%
%% Section S2
preamble
% Generate data for the plots
[x, y] = meshgrid(linspace(-8, 8, 100), linspace(-8, 8, 100));
z = Function04(x, y);
% Create Figure 1: Surface Plot
figure(1);
surf(x, y, z);
xlabel('x');
ylabel('y');
zlabel('z');
title('Surface Plot of Function04');
colorbar;
% Create Figure 2: Contour Plot
figure(2);
contour(x, y, z, 21);
xlabel('x');
ylabel('y');
title('Contour Plot of Function04');
colorbar;
axis equal;
% Find the minimal and maximal function values
min_value = min(z(:));
max_value = max(z(:));
% Define contour levels equally spaced between min and max values
contour_levels = linspace(min_value, max_value, 21);
% Set contour levels for the contour plot
contour_levels = unique(round(contour_levels, 2));
contour_levels = contour_levels(contour_levels >= min_value & contour_levels <= max_value);
% Adjust contour levels
contour_levels = sort(contour_levels, 'ascend');
% Display the contour levels
caxis([min_value, max_value]);
%
%%
% local function to be completed
function result=find_match(struct,key,fname)
% input arguments: ********************************************
% key: search string (char array or string)
% struct: structure array; the fields of struct contain 
%         strings or character arrays
% fname: name of the field of struct that is to be searched, 
%        (char array or string);
%        if empty or missing: all fields of struct are searched;
%        if field does not match a field name of struct, 
%        print out a warning and return an empty result;
%        the matching is case insensitive.
%
% *******************************************************
% if any input argument is not of the correct type,     *
% return an empty result and print a warning message,   *
% but do not terminate the execution!                   *
% *******************************************************
%
% output argument: **************************************
% result: sorted vector of indices of the elements of struct 
%         where a partial or complete match is found
% Example: result=find_match(struct,"01",'myfield') returns the 
%          indices of all elements of struct where the string "01" 
%          is found in field 'myfield'
% Example: result=find_match(struct,'Ab') returns the indices of
%          all elements of struct where any of the strings 'AB', 
%          'Ab', 'aB', 'ab' is found in any field
result=[];
% Check input arguments
if nargin < 2
    disp('Error: Insufficient input arguments.');
    return;
end
if ~isstruct(struct)
    disp('Warning: Input argument "struct" should be a structure array.');
    return;
end
if ~ischar(key) && ~isstring(key)
    disp('Warning: Input argument "key" should be a char array or a string.');
    return;
end
if nargin < 3
    fname = '';
end
% Iterate through each element of the structure array
for i = 1:length(struct)
    element = struct(i);
    if ~isempty(fname) && ~isfield(element, fname)
        disp(['Warning: Field "', fname, '" does not exist']);
        return;
    end
    % Get the field names of the current element
    field_names = fieldnames(element);
    % Iterate through each field of the current element
    for j = 1:length(field_names)
        if ~isempty(fname) && ~strcmp(field_names{j}, fname)
            continue; % Skip fields that don't match fname
        end
        % Check if the key is found in the field value
        if contains(lower(element.(field_names{j})), lower(key))
            result = [result, i];
            break;
        end
    end
end
% Sort the result vector
result = sort(result);
return
end
%
% local function for section S2
%
function z = Function04(x,y)
s=(x.^2+y.^2)/4000;
p=cos(x).*cos(y/sqrt(2));
z = s - p + 1;
z=z.*exp(-0.02*(x.^2+y.^2));
end
