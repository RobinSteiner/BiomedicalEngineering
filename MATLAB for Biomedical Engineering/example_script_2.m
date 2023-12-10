%% Standard preamble
preamble
%% open diary
% diary Diary2.txt % create text file with session diary
%% Arithmetical operations
x=rand % uniform random number in [0,1]
y=rand % uniform random number in [0,1]
x+y
plus(x,y) % functional form
x-y
minus(x,y)
x*y
times(x,y)
x/y
rdivide(x,y)
x\y
ldivide(x,y) % y/x=reciprocal of x/y
x^y
power(x,y)
%% Primes and factoring
n=100
primes(n) % primes up to n
isprime(37) % is a prime
isprime(39) % is not a prime
n=1234567890
factor(n) % prime factors
isprime(ans)
%% Rounding etc
round(534.1256) % round to nearest integer
round(534.7256) % round to nearest integer
round(-534.7256) % round to nearest integer
round(534.1256,1) % round to 1 digit after decimal point
round(534.1256,2) % round to 2 digits after decimal point
round(534.1256,-1) % round to the left of the decimal point
round(534.1256,-2) % round to the left of the decimal point
round(534.1256,6,'significant') % round to 6 significant digits
% floor and ceiling
floor(534.7256) % round to nearest integer towards -inf
ceil(534.7256) % round to nearest integer towards inf
floor(-534.7256) % round to nearest integer towards -inf
ceil(-534.7256) % round to nearest integer towards inf
%% Sign and absolute value
sign(-534.7256)
abs(-534.7256)
%% Roots
% square root
sqrt(x)
ans^2-x
% n-th root
nthroot(x,4)
ans^4-x
%% Trigonometric functions
[sin(1) cos(1) tan(1) cot(1)] % angle in rad
[sind(50) cosd(50) tand(50) cotd(50)] % angle in deg
%% Arc functions
[asin(1) acos(1) atan(0.5) acot(0.5)] % return angle in rad
[asind(1) acosd(1) atand(0.5) acotd(0.5)] % return angle in deg
%% exponential function
E=exp(1) % Euler's number
%% logarithms
log10(100) % decadic logarithm
log(E^3) % natural logarithm
log2(16) % dyadic logarithm
%% hyperbolic functions
[sinh(2) cosh(2) tanh(2) coth(2)]
%% area functions
[asinh(0.5) acosh(1.5) atanh(0.5) acoth(1.5)]
% 
%% All elementary functions take also complex arguments and vector/matrix arguments!
%
% complex arguments
z=11.2-23.6i % complex number
round(z)
round(z,-1)
floor(z)
ceil(z)
z=1+i;
abs(z) % norm of z
norm(z) % norm of z
angle(z) % angle of z
z/abs(z) % normalized z
sign(z) % normalized z
s=power(z,1/3)
s^3
exp(2i)
log(2i)
asin(2)
asin(2i)
% On vector/matrix arguments they are evaluated elementwise!
x=rand(1,10)-0.5
sign(x)
abs(x)
sin(x)
exp(x)
log(x)
A=rand(5)
asind(A)
%% Vector operations
x=randn(1,10) % row vector of standard normal random numbers
y=randn(1,10)
x+y % vector sum
x-y % vector difference
2*x % scalar product
dot(x,y) % inner (dot) product
cross(x(1:3),y(1:3)) % outer (cross) product, only for 3-vectors!
abs(x) % absolute values
norm(x) % Euclidean norm, length
sqrt(dot(x,x)) % same
u=x/norm(x) % unit vector
norm(u)
norm(x,1) % 1-norm of x
sum(abs(x)) % 1-norm of x
norm(x,inf) % Chebyshev norm of x
max(abs(x)) % Chebyshev norm of x
%% Pointwise=elementwise operations
x.*y % elementwise product
sum(x.*y)-dot(x,y) % almost 0
x*y'-dot(x,y) % exactly zero
x./y % elementwise quotient
z=x.^y % elementwise exponent. y can also be a scalar, see automatic expansion!
isreal(z) % check whether result is real
z(imag(z)==0) % find real elements in z
max(x) % maximal element
min(x) % minimal element
sum(x) % sum of elements
cumsum(x) % cumulative sums
prod(x) % product of elements
cumprod(x) % cumulative products
%% Automatic expansion
x
x+1 % scalar 1 is expanded to a vector of 1's
x+[1;2;3] % add 3 scalars -> matrix with 3 rows = [x+1;x+2;x+3];
v=y'
x+v % add column vector -> matrix <<<<<<<<<<<
repmat(x,length(y),1)+repmat(v,1,length(x)) % same result
x.*[1;2;3] % matrix with 3 rows =[1*x;2*x;3*x]
x.^[1;2;3] % matrix with 3 rows = [x.^1;x.^2;x.^3]
%% Matrix operations
n=5
A=rand(n) % nxn random matrix
B=rand(n)
A+B % matrix sum
A-B % matrix difference
A*B % matrix product
A^3 % matrix power A*A*A, only for square matrices!
mpower(A,3) % functional form
power(A,3) % elementwise power!
%
invA=inv(A) % inverse matrix, only for square matrices!
A*invA % identity matrix up to rounding errors
round(A*invA-eye(n),13) % null matrix up to rounding errors
det(A) % determinant, only for square matrices!
trace(A) % trace, only for square matrices!
rank(A) % rank, for arbitrary matrices
%
A/B % right division, B must be square
A*inv(B) % same, but different algorithm
A/B-A*inv(B) % not exactly zero
A\B % left division, A must be square
inv(A)*B % same, but different algorithm
A\B-inv(A)*B % not exactly zero
sum(A) % column sum -> row vector
sum(A,2) % row sum -> column vector
sum(sum(A)) % sum of all elements
sum(A(:)) % sum of all elements
eig(A) % eigenvalues
%
C=rand(3,5)
rank(C)
norm(C) % 2-norm of C = largest singular value
svd(C) % singular values
sv=sqrt(eig(C*C')) % singular values
max(sv) % 2-norm of C = largest singular value
norm(C,1) % 1-norm of C
max(sum(abs(C))) % 1-norm of C
norm(C,inf) % inf-norm of C
max(sum(abs(C'))) % inf-norm of C
norm(C,'fro') % Frobenius norm of C
norm(C(:),2) % Frobenius norm of C
%%
%
D=pinv(C) % Penrose-Moore inverse, pseudoinverse, C*D*C=C,D*C*D=D
round(C*D*C-C,13) % null matrix up to rounding errors
round(D*C*D-D,13) % null matrix up to rounding errors
round(pinv(A)-invA,13)  % null matrix up to rounding errors
%
A.*B % elementwise product
A./B % elementwise quotient
A.^3 % elementwise exponentiation
(A.^3).^(1/3)-A % null matrix up to rounding errors
sqrt(A) % elementwise square root
round(sqrt(A).^2-A,13) % null matrix up to rounding errors
% matrix functions
S=sqrtm(A) % matrix square root
round(S*S-A,13) % null matrix up to rounding errors
S=mpower(A,0.5) % same, but complex elements
C=mpower(A,1/3) % cube root
round(mpower(C,3)-A,13)
E=exp(A) % elementwise exponentiation
L=log(E) % elementwise logarithm
round(L-A,13) % null matrix up to rounding errors
E=expm(A) % matrix exponential
L=logm(E) % matrix logarithm
round(L-A,13) % null matrix up to rounding errors
%% Convert a matrix to a table
A=rand(4,5);
T=array2table(A,'VariableNames',["V1","V2","V3","V4","V5"],'RowNames',["R1","R2","R3","R4"])
%% Save (part of) work space
whos
% save entire workspace
save myspace
clear
load myspace
whos
% save just A and b
save mysystem A b
clear
load mysystem
whos
clear
load("mysystem","A","b")
whos
return