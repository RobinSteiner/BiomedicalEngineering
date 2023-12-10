%% MATLAB version and installed toolboxes
%%
% this is a comment
%
%% this starts a new section, comment is shown in bold
ver %version
%% Diary
%
%if exist('Diary1.txt')
%    delete('Diary1.txt') % delete file if it exists
%end
diary Diary1.txt % create text file with session diary
ls
%
%% Directories and files
%%
% show user path=your home directory
userpath
% change directory to userpath
cd(userpath)
% print current (working) directory
pwd
% list files in current directory
ls
dir
% list only MATLAB relevant files
what
% list files in preference directory
dir(prefdir)
% >>>>>>>>>>>>>>>>> MATLAB window: layout, home tab, editor tab, view tab
pwd % print working directory=show current position in directory tree
% create subdirectory
mkdir dummy
% go down in the tree
cd dummy
pwd
% go up in the tree
cd ..
pwd
% delete directory
rmdir dummy
%% Scripts
%%
% To execute a script, type its file name
preamble
% Look into code
edit preamble
% this file with extension .m is a script=sequence of MATLAB command
% get filename of currently executing script
mfilename
% Suppress output by terminating line with a semicolon
fn=mfilename;
% continue command on next line
fn=...
    mfilename;
%% help functions
%%
% open documentation center
% doclookfor
% search for help topics
lookfor('hyperbolic tangent')
lookfor('copy')
% get help for a specific function or command
help tanh
% command recall - up arrow
% clear command window
clc
% execution interrupt - ctrl c
% list all integers from 1 to 10 million
% [1:1e7]'
%% Matlab functions
%%
% print full path of function file
which tanh % built-in function, source code not available
which fminunc % toolbox function, source code available
% look at code
edit fminunc
% create file
edit newfile % private script or function
% type file
type fminunc 
more on
 type fminunc
more off
% check whether file exists
exist fminunc
exist('fminunc')
% delete file
delete newfile.m 
%% Search path
%%
% you usually have private scripts and private functions, maybe private versions of Matlab functions
%% scripts and functions are identified by their file name!
% the search path determines the order in which directories are searched for scripts or functions
%
% print search path
path
% add to search path
help addpath
% remove from search path
help rmpath
%% Work space
%%
% list objects in work space memory
who
% list objects with details
whos
% clear memory
clear
who
%% Variables
%%
% assign numerical values to variables
a=1
A=2
my_pi=pi % built-in function
% display numerical variables
format short % default
my_pi
format long
my_pi
format bank
my_pi
format short
my_pi
% show only numerical value
disp(my_pi)
% check whether variable exists
exist my_pi
exist('my_pi')
exist('My_pi') % names are case sensitive
% read numerical variables from keyboard
v=input('Enter a number: ')
%% Numerical data types
%%
i=int8(-55)
whos i
i=uint8(-55)
whos i
i=int16(-12345)
whos i
i=uint16(100000)
whos i
i=int32(-1234567890)
whos i
i=uint32(1234567890)
whos i
i=int64(-1e12)
whos i
i=uint64(1e20)
whos i
x=single(sqrt(2))
whos x
x*x-2
y=double(sqrt(2))
y*y-2
% check whether numerical variable
isnumeric(i)
isnumeric(x)
isnumeric(y)
% check whether floating point representation
isfloat(i)
isfloat(x)
isfloat(y)
%% Complex numbers
%%
x=2.33 % default is double!
whos x
y=-4.67
clear i
z=x+i*y
whos z
% check whether real
isreal(z)
isreal(x)
% real and imaginary part
real(z)
imag(z)
real(x)
imag(x)
%% Special real numbers
%%
% eps=smallest number such that 1+eps>1
eps
1+eps-1
1+eps/2-1
% smallest real number
realmin
% largest real number
realmax
%% Inf and NaN
%%
% infinity
w=realmax*1.1
isinf(w)
% negative infinity
w=-w
isinf(w)
% division by zero
10/0
-10/0
(-10)/(-0)
% not a number
0/0
isnan(inf-inf)
inf+inf
-inf-inf
% Stopped here October 7
%% Logical variables
%%
% value is false
f=(2==3)
f=logical(0)
f=false
% value is true
t=(2<3)
t=logical(pi) % any argument not zero gives true
t=true
% logical operations
% negation
~f
~t
% disjunction
t|f
% conjunction
t&f
% exclusive or
xor(f,f)
xor(t,t)
xor(t,f)
% equality and inequalities
2<3
3<=3
4>-17
4>=2
13==14
13~=14
x=[1 3 5] % numerical vector
y=[2 3 4] % numerical vector
l=(x==y) % logical vector
any(l)
all(l)
%% Basic data structures
%%
% numerical scalars
xx=3.49
whos xx
% vectors
xe=[] % empty vector
whos xe
size(xe)
% row vectors
xr=[1 2 3 4 6]
whos xr
size(xr)
% sort
xs=sort(xr)
xs=sort(xr,'descend')
% column vectors
xc=xr' % transposition
whos xc
size(xc)
% vectors with constant step size
x1=[2:3:17]' % first element 2, step size 3
x2=18:-2:0 % first element 18, step size -2
x3=linspace(0,1,21) % first element 0, last element 1, 21 elements
h=diff(x3) % step sizes
h=unique(h) % unique values
diff(h) % differences of unique values
% Matrices
M=[] % empty matrix
whos M
size(M)
% explicit definition
M=[1 2 3;4,5,6;7 8 9] % semicolon: next row - comma or blank: next column
whos M
size(M)
M1=zeros(3,5) % fill matrix with 0
size(M1)
M2=M1' % transpose
size(M2)
M3=ones(5) % fill matrix with 1
size(M3)
M4=diag(xr) % diagonal matrix
size(M4)
% identity matrix
M5=eye(6)
size(M5)
% check whether empty
isempty(M5)
% concatenate matrices and vectors
M6=[xr;xr] % vertical concatenation of row vectors
size(M6)
M7=[xc, xc] % horizontal concatenation of column vectors vectors
size(M7)
rng('default')
rng(0) % set random generator to default sequence
R1=rand(4,6) % matrix of uniform random numbers
R2=rand(4,2)
R3=[R1 R2] % horizontal concatenation
try
    R4=[R1;R2] % vertical concatenation
catch
    disp('error')
    exist('R4')
end
% flip
R2
fliplr(R2)
flipud(R2)
% replicate matrices
A=magic(3) % magic square
AA=repmat(A,3,2) % 3x2 copies of A
% block diagonal matrix
B=blkdiag(A,A,AA)
%% Addressing matrix and vector elements
%
%% addressing with indices
% The indices can be scalars or vectors
x1=[2:3:17]'
% last element
x1(end)
% last but one element
x1(end-1)
% every other element
x1(1:2:end)
% extend vector
x1(end+1:end+4)=[33 34 35 36]
% delete elements in a vector
x1([2 4 end-1])=[]
% random matrix
rng(1)
R=rand(4,6)
% single element
R(2,5)
% row
R(2,:)
% column
R(:,3)
% submatrix
R(2:3,3:5)
% collection of elements
R(2,[1 3 4 5])
R(3,1:2:end)
R
Rsave=R
% linear form of matrix
Rlin=R(:)
% linear indices
Rlin(1:8)
%% Logical addressing
x2=R(2,:) % extract values >= 0.5
l=x2>=0.5
x2(l)
x2(x2>=0.7)
% find indices
ind=find(x2>0.7)
x2(ind)=1
% logical addressing in a matrix
% display nonzero elements
R(R~=0)
% alternative
nonzeros(R)
% change some elements
R(R>0.8)=0
R(isnan(R))=0
% find indices in a matrix
help find
ind=find(R) % find linear indices of nonzero elements
% find row and column indices of nonzero elements
[ind,jnd]=find(R)
% find also the values
[ind,jnd,val]=find(R)
ind=find(R>0.5) % linear indices of elements > 0.5
% switch sign of those elements
R(ind)=-R(ind)
% find row and column indices of elements <0.5
[ind,jnd]=find(R<-0.5)
% set negative elements to 0
ind=find(R<0);
R(ind)=0
%% Character arrays and strings
%%
% character array
c='I am learning MATLAB'
ischar(c)
% extend character arrays
c1=[c,' today']
% subarray
c1(3:7)
% string
str(1)=string(c)
% extend string
str(1)=str(1)+" today"
% string arrays
str(2)="I already know MATLAB"
str(2,1)="I have forgotten MATLAB"
isstring(str)
ischar(str)
% convert to character array
c2=char(str(1,1))
ischar(c2)
% read characters from keyboard
c3=input('Enter characters: ','s')
%% Concatenation of strings and character arrays
%%
s11='abc'
s12='DEF'
S1=strcat(s11,s12) %character vector
S_1=[s11 s12] %character vector.. same result
%
s31=["abc","def"]
s32=["ghi","jkl"]
S3=strcat(s31,s32) %string array
join(S3) % note the blank!
S3(1)+S3(2) % no blank!
%% Find strings
%%
str1="the hardest part mastering strfind is to think of a proper character VECTOR"
strfind(str1,'v')
strfind(str1,'V')
strfind(str1,'c')
strfind(str1," ")
% findstr
% not recommended, searches the longer of the two input arguments for any
% occurrences of the shorter one k=findstr(str1, str2)
%% Related functions: startswith, endswith, contains; see documentation
%
%% String comparison
%%
str1='aaa'
str2='bbb'
str3='cccc'
str1==str2
try
    str1==str3 % error!
catch
    string(str1)==string(str3) % works
end
strcmp(str1,str3) % works for character arrays and strings
%% Display numerical result with text
disp(['The number pi is appoximately equal to ',num2str(pi,'%12.8f')]) % see num2str
%% Switch off diary
%%
diary off
