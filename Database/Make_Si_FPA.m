function Si_FPA = Make_Si_FPA

E0 = [1:100 110:10:200 220:20:300 350:50:500 600:100:2500 2750:250:5000 5500:500:30000];
N = 4000;

%% Basic
Si_FPA.Mat = 'Si_FPA';
Si_FPA.M = 28.0855;
Si_FPA.Z = 14;
Si_FPA.Density = 2.33; %g/cm^3
Si_FPA.Density = Si_FPA.Density*10^-24/Si_FPA.M*6.022*10^23; %#/A^3
Si_FPA.NvTPP = 4;
Si_FPA.NvSGS = 4;
Si_FPA.Eg = 1.1;
Si_FPA.Ep = 16.59;
Si_FPA.Ef = 12.5;
Si_FPA.Evb = 12.44;
Si_FPA.Affinity = 4.05;
Si_FPA.isMetal = false;

%% Elastic properties
Si_FPA.Elastic.x = zeros(numel(E0),1);
Si_FPA.Elastic.l_el = zeros(numel(E0),1);
Si_FPA.Elastic.l_tr = zeros(numel(E0),1);
Si_FPA.Elastic.x = E0;
Si_FPA.DECS.E0 = E0;
Si_FPA.Composition.Z = Si_FPA.Z;
Si_FPA.Composition.index = 1;

tic;
[data] = ElsepaRunner.RunElsepa(Si_FPA.Composition,E0);
toc
Si_FPA.DECS.x = data(1).x;
for i = 1:numel(E0)
    Si_FPA.Elastic.sigma_el(i) = data(i).sigma_el;
    Si_FPA.Elastic.l_el(i) = 1/(data(i).sigma_el*a0^2*Si_FPA.Density);
    Si_FPA.Elastic.l_tr(i) = 1/(data(i).sigma_tr1*a0^2*Si_FPA.Density);
    Si_FPA.DECS.y(:,i) = data(i).y;
end

%% Inelastic properties
osc.model = 'Mermin';
osc.A = [0.05 0.1 0.4 0.045 0.04 0.29 0.011 0.006 0.002]';
osc.G = [5 3 2.4 7 2 2.5 40 80 150]; 
osc.Om =[20 14 17.2 25 13.2 15.7 130 180 280];
osc.alpha = 1; 
osc.beps = 1;
osc.Ef = Si_FPA.Ef; 
osc.qtran = 0.01:0.01:20;
osc.eloss = eps:.1:110;
osc.egap = Si_FPA.Eg;

current_full_path = dbstack('-completenames');
current_file_name = dbstack;
if ispc
    ind = strfind(current_full_path(1).file,['Database\' current_file_name(1).file]);
    dirData = [current_full_path(1).file(1:ind-2) filesep 'Data\'];
elseif ismac || isunix
    ind = strfind(current_full_path(1).file,['Database/' current_file_name(1).file]);
    dirData = [current_full_path(1).file(1:ind-2) filesep 'Data/'];
end

% opt_elf = load([dirData '/opt/si.diel']);
% q = eps:0.1:10;
% eloss = eps:.5:500;
% tic
% [Au_FPA.ELF,~,~] = fpa_vector(q*a0,eloss/h2ev,opt_elf(:,1),opt_elf(:,4));
% toc
fpa_data = load([dirData 'elf_si_fpa.mat']);
Si_FPA.eloss = fpa_data.omega;
Si_FPA.q = fpa_data.q;
Si_FPA.ELF = fpa_data.elf;
Si_FPA.DIIMFP = zeros(N,2,numel(E0));
Si_FPA.l_in = zeros(numel(E0),1);
for i = 1:length(E0)
    if E0(i) > 2*Si_FPA.Eg + Si_FPA.Evb
        energy = E0(i) - Si_FPA.Eg - Si_FPA.Evb;
        eloss = Si_FPA.Eg:(energy-Si_FPA.Eg)/(N-1):energy;
        if energy < 10
            omega = Si_FPA.Eg:0.01:energy;
        else
            omega = Si_FPA.Eg:0.2:energy;
        end
        Si_FPA.DIIMFP(:,1,i) = eloss;
        [iimfp, diimfp_] = diimfp_interp_fpa(E0(i),omega,Si_FPA.ELF,Si_FPA.q,Si_FPA.eloss,osc);
        Si_FPA.DIIMFP(:,2,i) = interp1(omega,diimfp_,eloss);
        ind = isnan(Si_FPA.DIIMFP(:,2,i));
        Si_FPA.DIIMFP(ind,2,i) = diimfp_(end);
        Si_FPA.l_in(i) = 1/trapz(omega/h2ev,iimfp)*a0;
    else
        Si_FPA.l_in(i) = Inf;
    end
end

%% Ionisation shells
Si_FPA.Shells = {'2S1/2';'2P1/2';'2P3/2';'3S1/2';'3P1/2';};
Si_FPA.EB = [154;104;104;13;8;];

