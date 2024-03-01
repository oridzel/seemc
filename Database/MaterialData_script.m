clear

current_full_path = dbstack('-completenames');
current_file_name = dbstack;
ind = strfind(current_full_path(1).file,current_file_name(1).file);

% MaterialData = {};
load MaterialData.mat

%% Metals
% MaterialData.Au_FPA = Make_Au_FPA;
% MaterialData.Au_Mermin = Make_Au_Mermin;

% MaterialData.Cu_FPA = Make_Cu_FPA;
% MaterialData.Cu_Mermin = Make_Cu_Mermin;
% MaterialData.Cu_DFT_b0l0 = Make_Cu_DFT_b0l0;
% MaterialData.Cu_DFT_b1l0 = Make_Cu_DFT_b1l0;
% MaterialData.Cu_DFT_b1l1 = Make_Cu_DFT_b1l1;

% MaterialData.Ag_Mermin = Make_Ag_Mermin;

% MaterialData.W_Mermin = Make_W_Mermin;

%% Insulators
MaterialData.PMMA_FPA = Make_PMMA_FPA;
MaterialData.PMMA_Drude = Make_PMMA_Drude;
MaterialData.PMMA_MLL = Make_PMMA_MLL;

MaterialData.PS_FPA = Make_PS_FPA;
MaterialData.PS_Drude = Make_PS_Drude;
MaterialData.PS_MLL = Make_PS_MLL;

MaterialData.SiO2_FPA = Make_SiO2_FPA;
MaterialData.SiO2_Drude = Make_SiO2_Drude;
MaterialData.SiO2_MLL = Make_SiO2_MLL;

MaterialData.Si_FPA = Make_Si_FPA;
MaterialData.Si_Mermin = Make_Si_Mermin;
MaterialData.Si_DFT_b0l0 = Make_Si_DFT_b0l0;
MaterialData.Si_DFT_b1l0 = Make_Si_DFT_b1l0;
MaterialData.Si_DFT_b1l1 = Make_Si_DFT_b1l1;

% MaterialData.H2O = Make_H2O;

%% Save database
filename = [current_full_path(1).file(1:ind-2) filesep 'MaterialData.mat'];
save(filename,"MaterialData")