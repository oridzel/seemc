clear

inputpar.matName = 'PMMA_Drude';
inputpar.isMetal = false;
inputpar.numTrajectories = 1000;
inputpar.energy = [20 50 100 200 250 300 350 400 450 500 750 1000 1500];

s = SEEMC(inputpar);
s.onlyEscaped = true;
s.cbRef = true;

tic
s.simulate;
toc