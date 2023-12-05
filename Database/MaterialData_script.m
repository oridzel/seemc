clear

currentPath = pwd;
% load Database/MaterialData.mat
MaterialData = {};

MaterialData.PMMA = Make_PMMA;
MaterialData.Si = Make_Si;
MaterialData.Si_DL = Make_Si_DL;
MaterialData.SiO2 = Make_SiO2;
MaterialData.Si_DFT_b1l1 = Make_Si_DFT_b1l1;
MaterialData.Si_DFT_b0l0 = Make_Si_DFT_b0l0;
MaterialData.Si_DFT_b1l0 = Make_Si_DFT_b1l0;

% cd currentPath
save MaterialData.mat