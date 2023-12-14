clear

inputpar.matName = 'PMMA';
inputpar.isMetal = false;
inputpar.numTrajectories = 10000;
inputpar.trackTrajectories = false;
inputpar.energy = 70;
inputpar.onlyEscaped = true;

tic
s = SEEMC(inputpar);
s.simulate;
toc