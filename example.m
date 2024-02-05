clear

inputpar.matName = { 'Cu_FPA' };
inputpar.isMetal = { true };
% inputpar.thickness = 50;
inputpar.numTrajectories = 10000;
inputpar.energy = [50:10:100 200:50:500 600:100:1000];

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