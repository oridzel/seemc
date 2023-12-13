clear

inputpar.matName = 'Au';
inputpar.isMetal = true;
inputpar.numTrajectories = 50;
inputpar.trackTrajectories = true;
inputpar.energy = [100 500];
inputpar.onlyEscaped = false;

tic
s = SEEMC(inputpar);
s.simulate;
toc