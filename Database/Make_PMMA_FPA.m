function PMMA_FPA = Make_PMMA_FPA

E0 = [1:100 110:10:200 220:20:300 350:50:500 600:100:2500 2750:250:5000 5500:500:30000];
N = 4000;

%% Basic
PMMA_FPA = struct;
PMMA_FPA.Mat = 'PMMA_FPA';
PMMA_FPA.M = 100.114;
PMMA_FPA.Z = 54;
PMMA_FPA.Density = 1.188; %g/cm^3
PMMA_FPA.Density = PMMA_FPA.Density*10^-24/PMMA_FPA.M*6.022*10^23; %#/A^3
PMMA_FPA.NvTPP = 40;
PMMA_FPA.Eg = 6.7;
PMMA_FPA.Ep = 19.85;
PMMA_FPA.Evb = 15.8;
PMMA_FPA.Affinity = 4.5;
PMMA_FPA.Phonon.eps_zero = 3.9;
PMMA_FPA.Phonon.eps_inf = 2.2;
PMMA_FPA.Phonon.eloss = 0.1;
PMMA_FPA.isMetal = false;

%% Elastic properties
PMMA_FPA.Elastic.x = zeros(numel(E0),1);
PMMA_FPA.Elastic.l_el = zeros(numel(E0),1);
PMMA_FPA.Elastic.l_tr = zeros(numel(E0),1);
PMMA_FPA.Elastic.x = E0;
PMMA_FPA.DECS.E0 = E0;
PMMA_FPA.Composition.Z = [6 8 1];
PMMA_FPA.Composition.index = [5 2 8];

tic;
[data] = ElsepaRunner.RunElsepa(PMMA_FPA.Composition,E0,0);
toc
PMMA_FPA.DECS.x = data(1).x;
for i = 1:numel(E0)
    PMMA_FPA.Elastic.sigma_el(i) = data(i).sigma_el;
    PMMA_FPA.Elastic.l_el(i) = 1/(data(i).sigma_el*a0^2*PMMA_FPA.Density);
    PMMA_FPA.Elastic.l_tr(i) = 1/(data(i).sigma_tr1*a0^2*PMMA_FPA.Density);
    PMMA_FPA.DECS.y(:,i) = data(i).y;
end

%% Inelastic properties

osc.model = 'MerminLL';
osc.A = [0.06 0.04 0.06 0.16 0.05 0.07 0.08 0.06 0.01 0.02];
osc.G = [3.61  3.87  8.88  5.81  5.57  6.31 10.03  5.51 29.96 29.82];
osc.Om = [19.88 12.40 13.24 17.34 27.12 23.24 27.66 22.36 62.55 46.31];
osc.u = 6.2212;
osc.alpha = 1; 
osc.beps = 1;
osc.Ef = PMMA_FPA.Evb; 
osc.qtran = 0.01:0.01:20;
osc.eloss = eps:.1:110;
osc.egap = PMMA_FPA.Eg;

current_full_path = dbstack('-completenames');
current_file_name = dbstack;
if ispc
    ind = strfind(current_full_path(1).file,['Database\' current_file_name(1).file]);
    dirData = [current_full_path(1).file(1:ind-2) filesep 'Data\'];
elseif ismac || isunix
    ind = strfind(current_full_path(1).file,['Database/' current_file_name(1).file]);
    dirData = [current_full_path(1).file(1:ind-2) filesep 'Data/'];
end

% opt_elf = load([dirData '/opt/pmma.diel']);
% q = eps:0.1:10;
% eloss = eps:.5:500;
% tic
% [PMMA_FPA.ELF,~,~] = fpa_vector(q*a0,eloss/h2ev,opt_elf(:,1),opt_elf(:,4));
% toc
fpa_data = load([dirData 'elf_pmma_fpa.mat']);
PMMA_FPA.eloss = fpa_data.omega;
PMMA_FPA.q = fpa_data.q;
PMMA_FPA.ELF = fpa_data.elf;
PMMA_FPA.DIIMFP = zeros(N,2,numel(E0));
PMMA_FPA.l_in = zeros(numel(E0),1);
for i = 1:length(E0)
    if E0(i) > 2*PMMA_FPA.Eg + PMMA_FPA.Evb
        energy = E0(i) - PMMA_FPA.Eg - PMMA_FPA.Evb;
        eloss = PMMA_FPA.Eg:(energy-PMMA_FPA.Eg)/(N-1):energy;
        if energy < 10
            omega = PMMA_FPA.Eg:0.01:energy;
        else
            omega = PMMA_FPA.Eg:0.2:energy;
        end
        PMMA_FPA.DIIMFP(:,1,i) = eloss;
        [iimfp, diimfp_] = diimfp_interp_fpa(E0(i),omega,PMMA_FPA.ELF,PMMA_FPA.q,PMMA_FPA.eloss,osc);
        PMMA_FPA.DIIMFP(:,2,i) = interp1(omega,diimfp_,eloss);
        ind = isnan(PMMA_FPA.DIIMFP(:,2,i));
        PMMA_FPA.DIIMFP(ind,2,i) = diimfp_(end);
        PMMA_FPA.l_in(i) = 1/trapz(omega/h2ev,iimfp)*a0;
    else
        PMMA_FPA.l_in(i) = Inf;
    end
end

