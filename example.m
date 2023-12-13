clear

inputpar.matName = 'Au';
inputpar.isMetal = true;
inputpar.numTrajectories = 1000;
inputpar.trackTrajectories = false;
inputpar.energy = [20 50 100 200 500 1000];
inputpar.onlyEscaped = true;

tic
s = SEEMC(inputpar);
s.simulate;
toc