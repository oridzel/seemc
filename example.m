clear

inputpar = struct;
inputpar.name = 'Cu_FPA';
inputpar.theta_0 = 0;
% inputpar.thickness = 50;
inputpar.numTrajectories = 10000;
% inputpar.energy = [50:10:100 200:50:500 600:100:1000];
inputpar.onlyEscaped = true;
inputpar.cbRef = false;
inputpar.trackTrajectories = false;

s = SEEMC(inputpar);

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