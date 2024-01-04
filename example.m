clear

inputpar.matName = { 'SiO2' };
inputpar.isMetal = { false };
% inputpar.thickness = 50;
inputpar.numTrajectories = 100;
inputpar.energy = 100;

s = SEEMC(inputpar);
s.onlyEscaped = false;
s.cbRef = true;
s.trackTrajectories = true;

%% Run simulation
tic
s.simulate; 
toc

%% Trajectories
s.getTrajectories;
s.plotTrajectories(1);