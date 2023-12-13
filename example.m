clear

inputpar.matName = 'Au';
inputpar.isMetal = true;
inputpar.numTrajectories = 10000;
inputpar.trackTrajectories = false;
inputpar.energy = 100;
inputpar.onlyEscaped = false;

tic
s = SEEMC(inputpar);
s.simulate;
toc