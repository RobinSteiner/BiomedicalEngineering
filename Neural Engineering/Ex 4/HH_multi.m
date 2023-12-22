function [timeStep,vMat,mMat,hMat,nMat,iNaMat,iKMat,iLMat] = HH_multi()
%% ------------------------------------------------
%  Hodgkin & Huxley multi-compartment (stick) model
%  Run: [timeStep,vMat,mMat,hMat,nMat,iNaMat,iKMat,iLMat] = HH_multi()
%  The computation is done by a fixed step size forward or backward Euler method
%  Two stimulation modes: (1) Current injection into a single compartment ('iClamp')
%                         (2) Extracellular stimulation ('extracellular')
% ------------------------------------------------
%             
%                       _
%                  (2) / \     ||   
%                      \_/    _||_ (1)
%                             \  /                    
%   -------------------------- \/ --------------------------
%   |    1     |    2     |    3     |    4     |    5     |
%   --------------------------------------------------------
%
%
%% ------------- Step 1: PARAMETERS -------------
% Solver    
solver = 'BE'; % Type of solver, 'FE'==Forward Euler or 'BE'==Backward Euler

% Stimulus parameters
mode = 'extracellular'; % Stimulus mode, 'iClamp' or 'extracellular'
elecType = 'point'; % Electrode type, 'point' or 'disk' or 'FEM', for 'extracellular' only
dDisk = 50; % Electrode diameter, in um, for 'disk' only
rhoExt = 300; % Extracellular resistivity, in Ohm*cm, for 'point and 'disk' only
FEMfile = 'export.csv'; % File with FEM Solution, for 'FEM' only
I = -1; % Stimulus amplitude, in pA ('iClamp') or uA ('extracellular')
elecX = 0; % X-shift of cell relative to the electrode (at 0/0)), in um
elecY = 25; % Y-shift of cell relative to the electrode (at 0/0)), in um

% Temporal parameters
tStop = 10; % Total duration of simulation, in ms
tDel = 1; % Time when stimulus starts, in ms
tDur = 0.5; % Time of stimulus ON, in ms
tDt = 0.0125; % Time step, in ms

% Stick parameters, total length = (nComp-1) x lComp
nComp = 201; % Number of compartments
lComp = 10; % Compartment length, in um
rComp = 1; % Compartment radius, in um

c = 1; % Membrane specific capacitance, in uF/cm2, original=1
rhoA = 100; % Axial/Intracellular conductivity, in Ohm*cm, , original=100

% Temperature
temp = 6.3; % Model temperature, in Celsius, original=6.3

% Hodgkin & Huxley parameters
vInit = -65; % Membrane voltage to initialize the simulation, in mV
gNa = 120; % Sodium channel maximum conductivity, in mS/cm2, original=120
gK = 36; % Potassium channel maximum conductivity, in mS/cm2, original=36
gL = 0.3; % Leak channel maximum conductivity, in mS/cm2, original=0.3
eNa = 50; % Sodium reversal/equilibrium potential, in mV, original=50
eK = -77; % Potassium reversal/equilibrium potential, in mV, original=-77
eL = -54.3; % Leak reversal/equilibrium potential, in mV, original=-54.3

% Plot flags
plotV = 1; % Membrane voltage over time


%% ------------- Step 2: COMPUTE ADDITIONAL PARAMETERS -------------
% Central compartment of stick
centerComp = (nComp-1)/2+1;

% Axial resistance of each compartment, in kOhm
R = ((2*rhoA*lComp*1e-4)/(2*(rComp*1e-4)^2*pi))*1e-3;
% Compartment surface, in cm2
A = (2*rComp*pi*lComp)*1e-8;
% Compartment membrance capacitance, in uF
C = c*A;

% Compute time step
timeSteps = floor(tStop/tDt)+1;
timeStep = 0:tDt:tStop;

% Temperature adjustment
kTemp = 3^((temp-6.3)/10);

% Other constants
vAdd = 0.001;

% Set up tridiagonal matrix for fast computation of the inverse axial matrix (BE)
matInvDiag = [1+(tDt/C)*(1/R);ones(nComp-2,1)+(tDt/C)*(2/R);1+(tDt/C)*(1/R)]; % in (ms/uF)*(1/kOhm)==1/V
matInvOffDiag = ones(nComp-1,1)*-(tDt/C)*(1/R);
matInv = diag(matInvDiag,0) + diag(matInvOffDiag,-1) + diag(matInvOffDiag,1);

% Set up tridiagonal matrix for fast iAxial (FE) and iStim (extracellular) computation
matAxDiag = [-1/R;ones(nComp-2,1)*-2/R;-1/R]; % in 1/kOhm
matAxOffDiag = ones(nComp-1,1)*1/R;
matAxial = diag(matAxDiag,0) + diag(matAxOffDiag,-1) + diag(matAxOffDiag,1);

% Stimulus currents depending on stimulus type
switch mode
    case 'iClamp' % Simple current conversion
        iStim = zeros(nComp,1);
        iStim(centerComp) = 1e-6*I/A; % in 1e-6*pA/cm2==uA/cm2
    case 'extracellular' % Compute effect of extracellular electric fields on stick
        % Compute potentials at compartment centers, electrode is always located at elecX/elecY
        x = (-(centerComp-1)*lComp:lComp:(nComp-centerComp)*lComp)'+elecX; % in um
        y = ones(nComp,1)*elecY; % in um
        centers = [x,y];
        switch elecType % Compute extracellular potentials depending on electrode type
            case 'point'
                % Euklidean distance for each compartment center
                compDist = 1e-4*sqrt(sum(centers.^2,2)); % in cm
                % Analytical potentials (point source) for given distance
                potentials = 1e-3*(rhoExt*I)./(4*pi*compDist); % in 1e-3*(Ohm*cm*uA)/cm==mV
            case 'disk'
                % Axial and radial distance for each compartment center
                r = 1e-4*x; % in cm
                z = 1e-4*y; % in cm
                % Analytical potentials (disk electrode) for given distance
                rDisk = 1e-4*dDisk/2; % in cm
                potentials = 1e-3*(2*rhoExt*I)/(4*rDisk*pi) * asin((2*rDisk)./(sqrt((r-rDisk).^2+z.^2)+sqrt((r+rDisk).^2+z.^2))); % in 1e-3*(Ohm*cm*uA)/cm==mV
            case 'FEM'
                % Numerically computed potentials (arbitrary shape)
                imp = readmatrix(FEMfile); % Load exported solution
                xx = imp(:,11)*1e6; % FEM x-coordinates, conversion to um
                Ve = imp(:,10)*1e3; % FEM potentials, conversion to mV
                potentials = I*interp1(xx,Ve,x); % Interpolate to compartment centers
        end
        % Matrix times vector for stimulus current
        iStim = matAxial*potentials/A; % in ((1/kOhm)*mV)/cm2==uA/cm2
    otherwise
        fprintf('--- Unknown stimulus mode!')
        return
end

% Compute initial values
v0 = vInit; 
m0 = alphaM(v0,kTemp)/(alphaM(v0,kTemp)+betaM(v0,kTemp));
h0 = alphaH(v0,kTemp)/(alphaH(v0,kTemp)+betaH(v0,kTemp));
n0 = alphaN(v0,kTemp)/(alphaN(v0,kTemp)+betaN(v0,kTemp));

% Allocate memory for v, m, h, n
vMat = zeros(nComp,timeSteps);
mMat = zeros(nComp,timeSteps);
hMat = zeros(nComp,timeSteps);
nMat = zeros(nComp,timeSteps);

% Set Initial values
vMat(:,1) = v0;
mMat(:,1) = m0;
hMat(:,1) = h0;
nMat(:,1) = n0;


%% --------- Step 3: SOLVE ODEs & POSTPROCESSING ---------
tic 
for t=1:timeSteps-1
    % Get current states
    vVecT = vMat(:,t);
    mVecT = mMat(:,t);
    hVecT = hMat(:,t);
    nVecT = nMat(:,t);

    % Stimulus current
    iStimVec = zeros(nComp,1);
    if t>tDel/tDt && t<=(tDel+tDur)/tDt
        iStimVec = iStim; % in uA/cm2
    end
    
    % Ionic currents 
    % Sodium
    iNaVec = gNa*mVecT.^3.*hVecT.*(vVecT-eNa); % in (mS/cm2)*mV==uA/cm2
    % Potassium
    iKVec = gK*nVecT.^4.*(vVecT-eK); % in (mS/cm2)*mV==uA/cm2
    % Leak
    iLVec = gL*(vVecT-eL); % in (mS/cm2)*mV==uA/cm2
    iIonVec = iNaVec+iKVec+iLVec; % in uA/cm2

    % Update v, m, h and n
    switch solver
        case 'FE'
            % Axial currents
            iAxialVec = matAxial*vVecT/A; % in uA/cm2
            
            % Compute change of v
            vMat(:,t+1) = vVecT+(-iIonVec+iAxialVec+iStimVec)*(tDt/c); % in mV

            % Compute change of alphas and betas
            mMat(:,t+1) = mVecT+(alphaM(vVecT,kTemp).*(1-mVecT)-betaM(vVecT,kTemp).*mVecT)*tDt;
            hMat(:,t+1) = hVecT+(alphaH(vVecT,kTemp).*(1-hVecT)-betaH(vVecT,kTemp).*hVecT)*tDt;
            nMat(:,t+1) = nVecT+(alphaN(vVecT,kTemp).*(1-nVecT)-betaN(vVecT,kTemp).*nVecT)*tDt;
        case 'BE'
            % Additonal ionic contribution needed for BE
            % Sodium
            iNaAuxVec = gNa*mVecT.^3.*hVecT.*(vVecT+vAdd-eNa); % in (mS/cm2)*mV==uA/cm2
            % Potassium
            iKAuxVec = gK*nVecT.^4.*(vVecT+vAdd-eK); % in (mS/cm2)*mV==uA/cm2
            % Leak
            iLAuxVec = gL*(vVecT+vAdd-eL); % in (mS/cm2)*mV==uA/cm2
            % Sum
            rhsdidvVec = (iNaAuxVec-iNaVec+iKAuxVec-iKVec+iLAuxVec-iLVec)/vAdd; % in (uA/cm2)/mV

            % Compute change of v
            % Right hand side of matrix equation to be solved
            RHS = vVecT+(-iIonVec+rhsdidvVec.*vVecT+iStimVec)*(tDt/c);
            % Add ionic current contribution to left hand side
            matInv(1:1+size(matInv,1):end) = matInvDiag+rhsdidvVec*(tDt/c);
            
            % Solve matrix equation
            % vMat(:,t+1) = matInv\RHS; % in mV
            vMat(:,t+1) = sparse(matInv)\RHS; % in mV
        
            % Compute change of alphas and betas
            mMat(:,t+1) = (mVecT+tDt*alphaM(vMat(:,t+1),kTemp))./(1+tDt*(alphaM(vMat(:,t+1),kTemp)+betaM(vMat(:,t+1),kTemp)));
            hMat(:,t+1) = (hVecT+tDt*alphaH(vMat(:,t+1),kTemp))./(1+tDt*(alphaH(vMat(:,t+1),kTemp)+betaH(vMat(:,t+1),kTemp)));
            nMat(:,t+1) = (nVecT+tDt*alphaN(vMat(:,t+1),kTemp))./(1+tDt*(alphaN(vMat(:,t+1),kTemp)+betaN(vMat(:,t+1),kTemp)));
        otherwise
            fprintf('--- Unknown solver type!')
            return
    end
end
fprintf('--- Solving time was %.3f seconds.\n',toc) 

% Compute ionic current densities by using the state variables v, m, h and n
iNaMat = gNa*mMat.^3.*hMat.*(vMat-eNa); % in (mS/cm2)*mV==uA/cm2
iKMat = gK*nMat.^4.*(vMat-eK); % in (mS/cm2)*mV==uA/cm2
iLMat = gL*(vMat-eL); % in (mS/cm2)*mV==uA/cm2


%% --------- Step 4: VISUALIZE RESULTS ---------
if plotV==1
    % Membrane voltage over time
    figure('Position',[600,300,1000,420])
    subplot(1,2,1)
    hold all
    box off
    grid on
    set(gca,'TickDir','out')
    plot(timeStep,vMat,'b')
    plot(timeStep,vMat(centerComp,:),'r')
    rectangle('Position',[tDel,-100,tDur,10],'Edgecolor','r','Facecolor','r')
    xlim([0 tStop])
    ylim([-100 60])
    xlabel('Time (ms)')
    ylabel('Membrane voltage (mV)')
    
    % Membrane voltage over time (spatial)
    subplot(1,2,2)
    hold all
    box off
    set(gca,'TickDir','out')
    offset = -vInit;
    for i=1:nComp
        if i==centerComp
            plot(timeStep,vMat(i,:)+offset,'r');
        else
            plot(timeStep,vMat(i,:)+offset,'b');
        end
        offset = offset-lComp;
    end
    rectangle('Position',[tDel,offset-150,tDur,50],'Edgecolor','r','Facecolor','r')
    xlim([0 tStop])
    ylim([offset-150 150])
    xlabel('Time (ms)')
    ylabel('Location along stick (um)')
end


    %% ------------- DEFINE FUNCTIONS FOR ALPHAS AND BETAS -------------
    function am = alphaM(v,kT)
        am = kT*(0.1*(v+40))./(1-exp(-(v+40)/10));
    end
    function bm = betaM(v,kT)
        bm = kT*4*exp(-(v+65)/18);
    end
    function ah = alphaH(v,kT)  
        ah = kT*0.07*exp(-(v+65)/20);
    end
    function bh = betaH(v,kT)
        bh = kT*1./(1+exp(-(v+35)/10));
    end
    function an = alphaN(v,kT)
        an = kT*(0.01*(v+55))./(1-exp(-(v+55)/10));
    end
    function bn = betaN(v,kT)
        bn = kT*0.125*exp(-(v+65)/80);
    end
end