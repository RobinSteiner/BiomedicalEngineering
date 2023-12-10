preamble
%% open diary
% open_diary  % create text file with session diary
%
%% Functions 
%
%% Function handles = Pointers to functions
%
%% handle for built-in function
mysin=@sin
whos mysin
mysin(pi/3)
%  Use fplot to plot the sine function in the interval [-4,4]
%  fplot expects a handle = pointer as the first input argument
help fplot
try
    fplot(sin,[-4,4]) % kernel function
catch ME
    disp(ME.message)
    figure(1) % open figure window
    clf % clear figure window
    fplot(@sin,[-4,4]) % sine function is plotted in the interval [-4,4]
    hold on % keep previous plot
    fplot(mysin,[-4,4]) % plot again
end
%% handle for anonymous function
myfun=@(Z) 5*sin(Z).*Z; 
% myfun is the name of the handle, the function is vectorized and can be called
% with vector or matrix input
whos myfun
display(myfun)
methods(myfun)
func2str(myfun)
% Z is just a formal argument and does not occupy space in memory
% Its name is irrelevant!
%
myfun2=@(X) 5*sin(X).*X; % implements the same anonymous function!
y=linspace(0,2*pi,21);
myfun(y) % y is an actual argument in memory!
myfun2(y) % y is an actual argument in memory!
%% handle for parametric anonymous function
mypoly=@(X,C) C(1)*X^3+C(2)*X^2+C(3)*X+C(4); % polynomial with general coefficients
% X and C are formal arguments, they do not occupy space in memory!
% Their name is irrelevant!
display(mypoly)
% set values of the coefficients
c1=[1 -4 -2 4]
% evaluate mypoly
x=1;
mypoly(x,c1) % x and c are actual arguments in memory!
x=-2:1:5;
try
    mypoly(x,c1)
catch ME
    disp(ME.message)
    disp('Function handle mypoly is not vectorized')
end
% Vectorize mypoly
mypoly=@(X,C) C(1)*X.^3+C(2)*X.^2+C(3)*X+C(4);
mypoly(x,c1)
% Plot the function
figure(1)
clf
h1=fplot(@(X) mypoly(X,c1),[-2,5]) % fplot is a kernel function
% the argument of fplot is an anonymous function without handle(!) 
% it is a function of X alone, as reqired by fplot!
% h is a graphics handle (pointer to the graphics object)
% it allows you to change the properties of the object, for instance
set(h1, 'color','r')
hold on
%
% alternatively use plot instead of fplot
%
xx=linspace(-2,5,1000); % regular dense grid in the interval [-2,5]
plot(xx,mypoly(xx,c1),'color','r')
% axis labels, I use the LaTeX interpreter!
xlabel('$x$') % LaTeX math font
xlabel('x') % LaTex roman font
xlabel('x','interpreter','none') % Helvetica font
ylabel('$y$') % LaTeX math font
ylabel('y') % LaTex roman font
ylabel('y','interpreter','none') % Helvetica font
% how to change the default interpreter
get(0,'default') % 0 or groot is the root of the graphics tree
% if you do not know LaTeX syntax, change 'latex' to 'none' in startup.m
edit startup
htit=title('{\bf Polynomial function \texttt{mypoly}}')
% Use rgb triplet for setting the color
set(htit,'fontsize',18,'color',[0.3333 0.1059 0.5490]) % American violet= #551b8c 
legstr1=func2str(h1.Function)
legstr1="f(X)="+legstr1(5:end)
% c can be modified
c2=[1 -3 -2 3];
h2=fplot(@(X) mypoly(X,c2),[-2,5],'g') % plot with green line
legstr2=func2str(h2.Function);
legstr2="f(X)="+legstr2(5:end);
hleg=legend([h1 h2],[legstr1, legstr2],'location','northwest')
set(hleg,'fontsize',16)
% define cell array of function handles
myfuns{1}=@(X) mypoly(X,c1);
myfuns{2}=@(X) cos(X).*X;
myfuns{1}(-2:1:5)
myfuns{2}(y)
% anonymous function without named handle
figure(3)
clf
fplot(@(X) cos(X)+exp(-X),[0,2*pi])
xlim([0 2*pi])
hold on
h=yline(0);
set(h,'linestyle',':','color','black')
feval(@(X) cos(X)+exp(-X),y) % evaluate anonymous function
% Find zeros of the anonymous function
clear z
z(1)=fzero(@(X) cos(X)+exp(-X),2) % start search at x=2
z(2)=fzero(@(X) cos(X)+exp(-X),5) % start search at x=5
hold on
plot(z,0*z,'o','markersize',10)
xlabel('$x$')
ylabel('$y$')
title("Anonymous function")
hleg=legend('$y=f(x)$')
return
