clear

current_full_path = dbstack('-completenames');
current_file_name = dbstack;
ind = strfind(current_full_path(1).file,current_file_name(1).file);

% MaterialData = {};
load MaterialData.mat

% MaterialData.Au = Make_Au;
% MaterialData.Au_DL = Make_Au_DL;
% MaterialData.Cu = Make_Cu;
MaterialData.Cu_DFT_b0l0 = Make_Cu_DFT_b0l0;
MaterialData.Cu_DFT_b1l0 = Make_Cu_DFT_b1l0;
MaterialData.Cu_DFT_b1l1 = Make_Cu_DFT_b1l1;
% MaterialData.Ag = Make_Ag; 
% MaterialData.W = Make_W;

% MaterialData.PMMA_Drude = Make_PMMA_Drude;
% MaterialData.PMMA_MLL = Make_PMMA_MLL;
% MaterialData.PS_Drude = Make_PS_Drude;
MaterialData.PS_MLL = Make_PS_MLL;
% MaterialData.SiO2_Drude = Make_SiO2_Drude;
% MaterialData.Si = Make_Si;
% MaterialData.Si_DL = Make_Si_DL;
% MaterialData.Si_DFT_b0l0 = Make_Si_DFT_b0l0;
% MaterialData.Si_DFT_b1l0 = Make_Si_DFT_b1l0;
% MaterialData.Si_DFT_b1l1 = Make_Si_DFT_b1l1;
% MaterialData.H2O = Make_H2O;

filename = [current_full_path(1).file(1:ind-2) filesep 'MaterialData.mat'];
save(filename,"MaterialData")