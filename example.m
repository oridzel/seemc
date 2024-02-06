clear

inputpar.matName = { 'Au_FPA' };
% inputpar.thickness = 50;
inputpar.numTrajectories = 1000;
inputpar.energy = [50:10:100 200:50:500 600:100:1000];

s = SEEMC(inputpar);
s.theta_0 = 0;
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