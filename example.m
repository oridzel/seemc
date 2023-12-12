clear

inputpar.matName = 'Au';
inputpar.isMetal = true;
inputpar.numTrajectories = 100;
inputpar.trackTrajectories = false;
inputpar.energy = [100 102 105 110];

tic
s = SEEMC(inputpar);
s.simulate;
toc