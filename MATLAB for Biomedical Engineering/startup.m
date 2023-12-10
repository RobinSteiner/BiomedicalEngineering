% mystartup   Startup file
%   Change the name of this file to STARTUP.M. The file 
%   is executed when MATLAB starts up, if it exists 
%   anywhere on the path.  In this example, the
%   MAT-file generated during quitting using FINISHSAV
%   is loaded into MATLAB during startup.

%   Copyright 1984-2000 The MathWorks, Inc. 
%   $Revision: 1.4 $  $Date: 2000/06/01 16:19:26 $

% set default figure and axis properties
set(0,'DefaultFigurePaperType','a4');
set(0,'DefaultFigureUnits','centimeters')
set(0,'DefaultFigurePaperUnits','centimeters')
set(0,'DefaultAxesFontSize',14)
set(0,'DefaultAxesFontUnits','points')
set(0,'DefaultAxesFontName','Helvetica'); % 
set(0,'DefaultTextFontName','Helvetica'); % 
set(0,'DefaultTextFontSize',14);
set(0,'defaultTextInterpreter','latex')
set(0,'defaultLegendInterpreter','latex')
set(0,'defaultFigureWindowStyle','docked')
set(0,'defaultAxesTickLabelInterpreter','latex');  
% Set working directory to my directory for the VU
mypath="C:\Users\Robin\Documents\MATLAB\matlab_bme" % please replace the string with your actual path to MatlabVU
% the path must not contain blank spaces!!
%
userpath(mypath)
cd(mypath)
addpath(mypath)


