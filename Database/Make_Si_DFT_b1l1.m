function Si_DFT_b1l1 = Make_Si_DFT_b1l1

E0 = [1:100 150:50:500 600:100:2500 2750:250:5000 5500:500:30000];
N = 5000;

%% Basic
Si_DFT_b1l1.Mat = 'Si_DFT_b1l1';
Si_DFT_b1l1.M = 28.0855;
Si_DFT_b1l1.Z = 14;
Si_DFT_b1l1.Density = 2.33; %g/cm^3
Si_DFT_b1l1.Density = Si_DFT_b1l1.Density*10^-24/Si_DFT_b1l1.M*6.022*10^23; %#/A^3
Si_DFT_b1l1.NvTPP = 4;
Si_DFT_b1l1.NvSGS = 4;
Si_DFT_b1l1.Eg = 1.1;
Si_DFT_b1l1.Ep = 16.59;
Si_DFT_b1l1.Ef = 12.5;
Si_DFT_b1l1.Evb = 12.44;
Si_DFT_b1l1.Affinity = 4.05;

%% Elastic properties
% {
Si_DFT_b1l1.Elastic.x = zeros(numel(E0),1);
Si_DFT_b1l1.Elastic.l_el = zeros(numel(E0),1);
Si_DFT_b1l1.Elastic.l_tr = zeros(numel(E0),1);
Si_DFT_b1l1.Elastic.x = E0;
Si_DFT_b1l1.DECS.E0 = E0;
Si_DFT_b1l1.Composition.Z = 14;
Si_DFT_b1l1.Composition.index = 1;

tic;
[data] = ElsepaRunner.RunElsepa(Si_DFT_b1l1.Composition,E0);
toc
Si_DFT_b1l1.DECS.x = data(1).x;
for i = 1:numel(E0)
    Si_DFT_b1l1.Elastic.l_el(i) = 1/data(i).sigma_el/Si_DFT_b1l1.Density;
    Si_DFT_b1l1.Elastic.l_tr(i) = 1/data(i).sigma_tr1/Si_DFT_b1l1.Density;
    Si_DFT_b1l1.DECS.y(:,i) = data(i).y/trapz(data(i).x,data(i).y);
end
%}

%% Inelastic properties
osc.model = 'Mermin';
osc.A = [0.05 0.1 0.4 0.045 0.04 0.29 0.011 0.006 0.002]';
osc.G = [5 3 2.4 7 2 2.5 40 80 150]; 
osc.Om =[20 14 17.2 25 13.2 15.7 130 180 280];
osc.alpha = 1; 
osc.beps = 1;
osc.Ef = Si_DFT_b1l1.Ef; 
osc.qtran = 0.01:0.01:20;
osc.eloss = eps:.1:110;
osc.egap = Si_DFT_b1l1.Eg;

current_full_path = dbstack('-completenames');
current_file_name = dbstack;
if ispc
    ind = strfind(current_full_path(1).file,['Database\' current_file_name(1).file]);
    dirData = [current_full_path(1).file(1:ind-2) filesep 'Data\'];
elseif ismac || isunix
    ind = strfind(current_full_path(1).file,['Database/' current_file_name(1).file]);
    dirData = [current_full_path(1).file(1:ind-2) filesep 'Data/'];
end
load([dirData 'elf_si_dft_ocean_0_01_b1l1.mat'])
Si_DFT_b1l1.ELF = elf;
Si_DFT_b1l1.eloss = omega;
Si_DFT_b1l1.q = q;

Si_DFT_b1l1.DIIMFP = zeros(N,2,numel(E0));
Si_DFT_b1l1.l_in = zeros(numel(E0),1);
for i = 1:length(E0)
    if E0(i) > 2*Si_DFT_b1l1.Eg + Si_DFT_b1l1.Evb
        energy = E0(i) - Si_DFT_b1l1.Eg - Si_DFT_b1l1.Evb;
        osc.eloss = Si_DFT_b1l1.Eg:(energy-Si_DFT_b1l1.Eg)/(N-1):energy;
        Si_DFT_b1l1.DIIMFP(:,1,i) = osc.eloss;
        [iimfp, Si_DFT_b1l1.DIIMFP(:,2,i)] = ndiimfp(osc,E0(i),elf,q,omega);
        Si_DFT_b1l1.l_in(i) = 1/trapz(osc.eloss/h2ev,iimfp)*a0;
    else
        Si_DFT_b1l1.l_in(i) = Inf;
    end
end

%% Ionisation shells
Si_DFT_b1l1.Shells = {'2S1/2';'2P1/2';'2P3/2';'3S1/2';'3P1/2';};
Si_DFT_b1l1.EB = [154;104;104;13;8;];

