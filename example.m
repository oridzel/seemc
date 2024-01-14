clear

inputpar.matName = { 'SiO2_Drude' };
inputpar.isMetal = { false };
% inputpar.thickness = 50;
inputpar.numTrajectories = 1000;
inputpar.energy = [20:10:100 200 500 700 1000];

s = SEEMC(inputpar);
s.onlyEscaped = true;
s.cbRef = false;
s.trackTrajectories = false;

%% Run simulation
tic
s.simulate; 
toc

%% Trajectories
%{
s.getTrajectories;
s.plotTrajectories(1);
%}

%% Yields
% {
s.calculateYields;
s.plotYields;  
%}