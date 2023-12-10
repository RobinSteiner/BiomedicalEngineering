preamble
%% more data structures
%
%% n-dimensional numerical arrays, n>2
% Example with n=3
R=rand(2,3,4) % 4 2x3 matrices
whos R
S=sum(R) % column sums of the four matrices
whos S
T=squeeze(S) % remove singleton dimensions
size(T)
S=sum(R,2) % row sums of the four matrices
whos S
T=squeeze(S) % remove singleton dimensions
size(T)
S=sum(R,3) % sum of the four matrices
whos S
%% Cell arrays
% Data structures that can store objects of different classes or sizes
clear
% cell array of matrices
C{1,1}=magic(4);
C{1,2}=magic(5);
C{2,1}=magic(6);
% creating a subarray
C(1,:) % this is again a cell array!
% accessing the cell content
C{1,:} % this are the contents(!) of the cells in the first row of C
C{1,1}(2,2) % an element of the matrix in C{1,1}
% cell arrays can contain anything!
D{1}=pi; % a numerical value
D{2}=1:10; % a row vector
D{3}=rand(5,1); % a column vector
D{4}=C; % a cell array
D(4)
D{4}
D{4}{1,2} % C{1,2}
D{4}{1,2}(3,3) % element (3,3) of C{1,2}
% display contents
celldisp(D)
% convert matrix to cell array
A=rand(4)
num2cell(A)
% break up matrix into smaller cells
E=mat2cell(A,[2 2],[2 2])
% convert back to matrix
B=cell2mat(E)
%%