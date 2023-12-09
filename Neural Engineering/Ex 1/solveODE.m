function [] = solveODE()
%% ------------------------------------------------
%  Implementation of different ways of solving the ODE y' = -20*y
%  The analytical solution of the ODE is y(t) = exp(-20*t)
%  The computation is done analytically or by a fixed step size forward or backward Euler method

%% ------------- Step 1: PARAMETERS -------------
y0 = 1; % Initial condition
tStop = 1; % Total duration of simulation
tDt = 0.025; % Time step


%% ------------- Step 2: COMPUTE ADDITIONAL PARAMETERS -------------
% Compute time step
timeSteps = round(tStop/tDt)+1;
timeStep = 0:tDt:tStop;


%% ------------- Step 3: ANALYTICAL -------------
tDtAn = 0.001; % Fine time step for reference solution
timeStepAn = 0:tDtAn:tStop;
yAnVec = exp(-20*timeStepAn);


%% ------------- Step 4: FORWARD EULER ----------
yFEVec = zeros(1,timeSteps); % Allocate memory
yFEVec(1) = y0; % Set initial condition
% Loop over time
for t=1:timeSteps-1
    yFEVec(t+1) = yFEVec(t) + (-20*yFEVec(t))*tDt;
end


%% ------------- Step 5: BACKWARD EULER ---------
yBEVec = zeros(1,timeSteps); % Allocate memory
yBEVec(1) = y0; % Set initial condition

% Implicit BE to explicit form
% y1=y0+(-20*y1)*h
% y1+20*y1*h=y0
% y1*(1+20h)=y0
% y1=y0/(1+20h)

% Loop over time
for t=1:timeSteps-1
    yBEVec(t+1) = yBEVec(t)/(1+20*tDt);
end
   
   
%% --------- Step 6: VISUALIZE RESULTS ---------
% y over time
figure('Color','w','DefaultAxesFontSize',7)
hold all
grid on
set(gca,'TickDir','out')
plot(timeStepAn,yAnVec)
plot(timeStep,yFEVec)
plot(timeStep,yBEVec)
legend('Analytical','FE','BE')
xlim([0,tStop])
ylim([min(-0.1,min(yFEVec)),max(1.1,max(yFEVec))])
xlabel('Time')
ylabel('Y')
    
end
