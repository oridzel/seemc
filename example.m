clear

inputpar.matName = 'Si_DFT_b1l1';
inputpar.isMetal = false;
inputpar.numTrajectories = 1000;
inputpar.trackTrajectories = false;
inputpar.energy = 100;
inputpar.onlyEscaped = true;

tic
s = SEEMC(inputpar);
s.simulate;
toc