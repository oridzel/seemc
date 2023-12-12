clear

inputpar.matName = 'Au';
inputpar.isMetal = true;
inputpar.numTrajectories = 100;
inputpar.trackTrajectories = false;
inputpar.energy = 100;

s = SEEMC(inputpar);
s.simulate;