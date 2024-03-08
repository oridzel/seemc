function SiO2_FPA = Make_SiO2_FPA

E0 = [1:100 110:10:200 220:20:300 350:50:500 600:100:2500 2750:250:5000 5500:500:30000];
N = 4000;

%% Basic
SiO2_FPA = struct;
SiO2_FPA.Mat = 'SiO2_FPA';
SiO2_FPA.M = 60.008;
SiO2_FPA.Z = 30;
SiO2_FPA.Density = 2.19; %g/cm^3
SiO2_FPA.Density = SiO2_FPA.Density*10^-24/SiO2_FPA.M*6.022*10^23; %#/A^3
SiO2_FPA.NvTPP = 16;
SiO2_FPA.Eg = 9.1;
SiO2_FPA.Ep = 22.02;
SiO2_FPA.Evb = 10;
SiO2_FPA.Affinity = 1.1;
SiO2_FPA.Phonon.eps_zero = 3.84; % static dielectric constant
SiO2_FPA.Phonon.eps_inf = 2.25; % high-frequency dielectric constant (the square of the static refractive index)
SiO2_FPA.Phonon.eloss = [63 153]*1e-3;
SiO2_FPA.isMetal = false;

%% Elastic properties
SiO2_FPA.Elastic.x = zeros(numel(E0),1);
SiO2_FPA.Elastic.l_el = zeros(numel(E0),1);
SiO2_FPA.Elastic.l_tr = zeros(numel(E0),1);
SiO2_FPA.Elastic.x = E0;
SiO2_FPA.DECS.E0 = E0;
SiO2_FPA.Composition.Z = [14 8];
SiO2_FPA.Composition.index = [1 2];

tic;
[data] = ElsepaRunner.RunElsepa(SiO2_FPA.Composition,E0,0);
toc
SiO2_FPA.DECS.x = data(1).x;
for i = 1:numel(E0)
    SiO2_FPA.Elastic.sigma_el(i) = data(i).sigma_el;
    SiO2_FPA.Elastic.l_el(i) = 1/(data(i).sigma_el*a0^2*SiO2_FPA.Density);
    SiO2_FPA.Elastic.l_tr(i) = 1/(data(i).sigma_tr1*a0^2*SiO2_FPA.Density);
    SiO2_FPA.DECS.y(:,i) = data(i).y;
end

%% Inelastic properties

osc.model = 'MerminLL';
osc.A = [0.06 0.03 0.03 0.04 0.18 0.11 0.09 0.01 0.07 0.01 0.01 0.03 0.01];
osc.G = [2.00  2.00  2.00  3.98  5.27  5.89  7.62 13.62 12.52 16.41  5.26 11.13 14.06];
osc.Om = [10.00 12.00 15.00 14.27 19.23 22.72 26.24 60.45 32.31 74.19 47.94 40.89 52.92];
osc.u = 9.0747;
osc.alpha = 1; 
osc.beps = 1;
osc.Ef = SiO2_FPA.Evb; 
osc.qtran = 0.01:0.01:20;
osc.eloss = eps:.1:110;
osc.egap = SiO2_FPA.Eg;

current_full_path = dbstack('-completenames');
current_file_name = dbstack;
if ispc
    ind = strfind(current_full_path(1).file,['Database\' current_file_name(1).file]);
    dirData = [current_full_path(1).file(1:ind-2) filesep 'Data\'];
elseif ismac || isunix
    ind = strfind(current_full_path(1).file,['Database/' current_file_name(1).file]);
    dirData = [current_full_path(1).file(1:ind-2) filesep 'Data/'];
end

% opt_elf = load([dirData '/opt/sio2.diel']);
% q = eps:0.1:10;
% eloss = eps:.5:500;
% tic
% [PMMA_FPA.ELF,~,~] = fpa_vector(q*a0,eloss/h2ev,opt_elf(:,1),opt_elf(:,4));
% toc
fpa_data = load([dirData 'elf_sio2_fpa.mat']);
SiO2_FPA.eloss = fpa_data.omega;
SiO2_FPA.q = fpa_data.q;
SiO2_FPA.ELF = fpa_data.elf;
SiO2_FPA.DIIMFP = zeros(N,2,numel(E0));
SiO2_FPA.l_in = zeros(numel(E0),1);
for i = 1:length(E0)
    if E0(i) > 2*SiO2_FPA.Eg + SiO2_FPA.Evb
        energy = E0(i) - SiO2_FPA.Eg - SiO2_FPA.Evb;
        eloss = SiO2_FPA.Eg:(energy-SiO2_FPA.Eg)/(N-1):energy;
        if energy < 10
            omega = SiO2_FPA.Eg:0.01:energy;
        else
            omega = SiO2_FPA.Eg:0.2:energy;
        end
        SiO2_FPA.DIIMFP(:,1,i) = eloss;
        [iimfp, diimfp_] = diimfp_interp_fpa(E0(i),omega,SiO2_FPA.ELF,SiO2_FPA.q,SiO2_FPA.eloss,osc);
        SiO2_FPA.DIIMFP(:,2,i) = interp1(omega,diimfp_,eloss);
        ind = isnan(SiO2_FPA.DIIMFP(:,2,i));
        SiO2_FPA.DIIMFP(ind,2,i) = diimfp_(end);
        SiO2_FPA.l_in(i) = 1/trapz(omega/h2ev,iimfp)*a0;
    else
        SiO2_FPA.l_in(i) = Inf;
    end
end
