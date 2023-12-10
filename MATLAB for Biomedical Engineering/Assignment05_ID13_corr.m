%% Assignment 05
%% Section S1
preamble
X=20;
xx=linspace(0,X,1000);
% Compute and plot derivatives
%
syms x n
f(x)=besselj(n,x);
[f0,f1]=compderiv();
hfig=plotderiv(xx,f0,f1);
%
n=1;
f(x)=besselj(n,x);
%
[f0,f1]=compderiv(eval(f),1);
hfig=plotderiv(xx,f0,f1);
%
[f0,f1]=compderiv(f,1,1.5);
hfig=plotderiv(xx,f0,f1);
%
[f0,f1,f2]=compderiv(f,1,'2');
hfig=plotderiv(xx,f0,f1,f2);
%
[f0,f1,f2]=compderiv(f,1,{'2'});
hfig=plotderiv(xx,f0,f1,f2);
%
f0=compderiv(f);
hfig=plotderiv(xx,f0);
%
[f0,f1,f2]=compderiv(f,1,2);
hfig=plotderiv(xx,f0,f1,f2);
%
n=2;
f(x)=besselj(n,x);
[f0,f1,f2,f3]=compderiv(f,[1,2,3]);
hfig=plotderiv(xx,f0,f1,f2,f3);
%
n=3;
f(x)=besselj(n,x);
[f0,f1,f2,f3,f4,f5]=compderiv(f,[1,2],{3,4});
hfig=plotderiv(xx,f0,f1,f2,f3,f4,f5);
return %%%%% end of script
%
%% local functions
%
function [f0,varargout]=compderiv(f,varargin)
%
%% Compute derivatives of a symbolic function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input 
% f...symbolic function, must not contain more than one symbolic variable
% varargin...Orders of the derivatives, must be non-negative integers
%            If it is empty, set it to zero and return only f0
%% Output 
% f0...handle of anonymous function that implements f
% varargout...handles of anonymous functions that implement
%             the requested derivatives or empty matrices if necessary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Check that:
% - at least one input argument is present; if not, print a message and
%   return empty matrices
% - f is a symbolic function; if not, print a messsage and return empty matrices
% - f contains only a single symbolic variable; if not, print a message and
%   return empty matrices
% - varargin can be converted to a numerical vector; if not, print a message, 
%   return f0 and empty matrices  
% - varargin contains only non-negative integers; if not, print a message, 
%   return f0 and empty matrices  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f0 = [];
varargout=cell(1, nargout-1);
if nargin==0
    disp(">>>>> compderiv: Supply at least one input argument")
    return
end
if isequal(class(f),'sym')
    disp(">>>>> compderiv: First input argument must be a symbolic function")
    return
end
if length(symvar(f)) > 1
    disp('>>>>> compderiv: The function f must contain only a single symbolic variable');
    return;
end
f0=matlabFunction(f);
if nargin>1
    varinp=varargin;
    indcell=cellfun(@iscell,varinp);
    try
        varinp(indcell)=cellfun(@cell2mat,varinp(indcell),'UniformOutput',false);
    catch
        disp(">>>>> compderiv: Different data types in cell aray")
        return
    end
    try
        orders=cell2mat(varinp);
    catch
        disp(">>>>> compderiv: Cannot convert variable arguments to a vector")
        return
    end
    if any(orders < 0) || ~isequal(orders, floor(orders))
        disp('>>>>> compderiv: The orders must be non-negative integers');
        return;
    end
    varargout(1:length(orders)) = arrayfun(@(x) matlabFunction(diff(f, x)), orders, 'UniformOutput', false);
end

end %%%%%% end of local function compderiv
% 
function hfig=plotderiv(xx,varargin) % RF, 2023-11-10
hfig=[];
if nargin==1,return,end
nhan=nargin-1;
handles=varargin;
for ihan=1:nhan
    han=handles{ihan};
    if ~isequal(class(han),'function_handle')
        continue
    else
        s(ihan)=string(inputname(ihan+1));
        if isempty(hfig),hfig=figure;hold on,end
        plot(xx,han(xx),'displayname',s(ihan))
    end
end
if isempty(hfig)
    disp('>>>>> plotderiv: No figure created')
    return
end
box on 
grid on
xlabel('$x$')
ylabel('$y$')
title("Function and derivatives: "+join(s,','));
legend
return
end %%%%% end of local function plotderiv


