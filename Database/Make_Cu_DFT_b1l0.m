function Cu_DFT_b1l0 = Make_Cu_DFT_b1l0

E0 = [1:100 110:10:200 220:20:300 350:50:500 600:100:2500 2750:250:5000 5500:500:30000];
N = 4000;

%% Basic
Cu_DFT_b1l0.Mat = 'Cu_DFT_b1l0';
Cu_DFT_b1l0.M = 63.546;
Cu_DFT_b1l0.Z = 29;
Cu_DFT_b1l0.Density = 8.96; %g/cm^3
Cu_DFT_b1l0.Density = Cu_DFT_b1l0.Density*10^-24/Cu_DFT_b1l0.M*6.022*10^23; %#/A^3
Cu_DFT_b1l0.NvTPP = 11;
Cu_DFT_b1l0.Ep = 35.87;
Cu_DFT_b1l0.Ef = 8.7;
Cu_DFT_b1l0.Wf = 4.65;
Cu_DFT_b1l0.isMetal = true;

%% Elastic properties
% {
Cu_DFT_b1l0.Elastic.x = zeros(numel(E0),1);
Cu_DFT_b1l0.Elastic.l_el = zeros(numel(E0),1);
Cu_DFT_b1l0.Elastic.l_tr = zeros(numel(E0),1);
Cu_DFT_b1l0.Elastic.x = E0;
Cu_DFT_b1l0.DECS.E0 = E0;
Cu_DFT_b1l0.Composition.Z = Cu_DFT_b1l0.Z;
Cu_DFT_b1l0.Composition.index = 1;

tic;
[data] = ElsepaRunner.RunElsepa(Cu_DFT_b1l0.Composition,E0);
toc
Cu_DFT_b1l0.DECS.x = data(1).x;
for i = 1:numel(E0)
    Cu_DFT_b1l0.Elastic.sigma_el(i) = data(i).sigma_el;
    Cu_DFT_b1l0.Elastic.l_el(i) = 1/(data(i).sigma_el*a0^2)/Cu_DFT_b1l0.Density;
    Cu_DFT_b1l0.Elastic.l_tr(i) = 1/(data(i).sigma_tr1*a0^2)/Cu_DFT_b1l0.Density;
    Cu_DFT_b1l0.DECS.y(:,i) = data(i).y;
end
%}

%% Inelastic properties
osc.model = 'Mermin';
osc.A = [0.01 0.02 0.06 0.11 0.38 0.15 0.07 0.05 0.08 0.1]';
osc.G = [0.69 0.8 2.36 6.59 10.82 4.89 40.77 19.15 14.12 161.76]; 
osc.Om =[3.46 4.2 7.54 28.34 19.39 10.48 66.56 47.52 36.13 126.04];
osc.alpha = 1; 
osc.beps = 1;
osc.Ef = Cu_DFT_b1l0.Ef; 
osc.qtran = 0.01:0.01:20;
osc.eloss = eps:.1:110;
osc.egap = 0;

current_full_path = dbstack('-completenames');
current_file_name = dbstack;
if ispc
    ind = strfind(current_full_path(1).file,['Database\' current_file_name(1).file]);
    dirData = [current_full_path(1).file(1:ind-2) filesep 'Data\'];
elseif ismac || isunix
    ind = strfind(current_full_path(1).file,['Database/' current_file_name(1).file]);
    dirData = [current_full_path(1).file(1:ind-2) filesep 'Data/'];
end
dft_data = load([dirData 'elf_cu_dft+bse_20kp_ocean_0_01_b1l0.mat']);
Cu_DFT_b1l0.ELF = dft_data.elf;
Cu_DFT_b1l0.eloss = dft_data.omega;
Cu_DFT_b1l0.q = dft_data.q;

Cu_DFT_b1l0.DIIMFP = zeros(N,2,numel(E0));
Cu_DFT_b1l0.l_in = zeros(numel(E0),1);
for i = 1:length(E0)
    if E0(i) > Cu_DFT_b1l0.Ef
        energy = E0(i) - Cu_DFT_b1l0.Ef;
        eloss = eps:(energy-eps)/(N-1):energy;
        if energy < 10
            osc.eloss = eps:0.01:energy;
        else
            osc.eloss = eps:0.2:energy;
        end
        Cu_DFT_b1l0.DIIMFP(:,1,i) = eloss;
        [iimfp, diimfp] = ndiimfp(osc,E0(i),Cu_DFT_b1l0.ELF,Cu_DFT_b1l0.q,Cu_DFT_b1l0.eloss);
        Cu_DFT_b1l0.DIIMFP(:,2,i) = interp1(osc.eloss,diimfp,eloss);
        Cu_DFT_b1l0.l_in(i) = 1/trapz(osc.eloss/h2ev,iimfp)*a0;
    else
        Cu_DFT_b1l0.l_in(i) = Inf;
    end
end

%% Ionisation shells
Cu_DFT_b1l0.Shells = {'2S1/2';'2P1/2';'2P3/2';'3S1/2';'3P1/2';'3P3/2';'3D3/2';'3D5/2';'4S1/2';};
Cu_DFT_b1l0.EB = [1103;958;938;127;82;80;11;10;8;];

