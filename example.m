clear

inputpar.matName = { 'Si_DFT_b1l1' };
inputpar.isMetal = { false };
% inputpar.thickness = 50;
inputpar.numTrajectories = 100;
inputpar.energy = 100; %[20:10:100 200 500 700 1000];

s = SEEMC(inputpar);
s.onlyEscaped = false;
s.cbRef = false;
s.trackTrajectories = true;

%% Run simulation
tic
s.simulate; 
toc

%% Trajectories
% {
s.getTrajectories;
s.plotTrajectories(1);
%}