clear

current_full_path = dbstack('-completenames');
current_file_name = dbstack;
ind = strfind(current_full_path(1).file,current_file_name(1).file);

MaterialData = {};

MaterialData.Au = Make_Au;
MaterialData.Au_DL = Make_Au_DL;
MaterialData.Cu = Make_Cu;
MaterialData.Ag = Make_Ag;

MaterialData.PMMA = Make_PMMA;
MaterialData.SiO2 = Make_SiO2;
MaterialData.Si = Make_Si;
MaterialData.Si_DL = Make_Si_DL;
MaterialData.Si_DFT_b1l1 = Make_Si_DFT_b1l1;

filename = [current_full_path(1).file(1:ind-2) filesep 'MaterialData.mat'];
save(filename,"MaterialData")