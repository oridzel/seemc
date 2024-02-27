function PS_FPA = Make_PS_FPA

E0 = [1:100 110:10:200 220:20:300 350:50:500 600:100:2500 2750:250:5000 5500:500:30000];
N = 4000;

%% Basic
PS_FPA = struct;
PS_FPA.Mat = 'PS_FPA';
PS_FPA.M = 104.1440;
PS_FPA.Z = 56;
PS_FPA.Density = 1.05; %g/cm^3
PS_FPA.Density = PS_FPA.Density*10^-24/PS_FPA.M*6.022*10^23; %#/A^3
PS_FPA.NvTPP = 40;
PS_FPA.Eg = 4.5;
PS_FPA.Ep = 18.3;
PS_FPA.Evb = 17;
PS_FPA.Affinity = 3.5;
PS_FPA.Phonon.eps_zero = 2.5; % static dielectric constant
PS_FPA.Phonon.eps_inf = 1.01; % high-frequency dielectric constant (the square of the static refractive index)
PS_FPA.Phonon.eloss = 0.3;
PS_FPA.isMetal = false;

%% Elastic properties
PS_FPA.Elastic.x = zeros(numel(E0),1);
PS_FPA.Elastic.l_el = zeros(numel(E0),1);
PS_FPA.Elastic.l_tr = zeros(numel(E0),1);
PS_FPA.Elastic.x = E0;
PS_FPA.DECS.E0 = E0;
PS_FPA.Composition.Z = [6 1];
PS_FPA.Composition.index = [8 8];

tic;
[data] = ElsepaRunner.RunElsepa(PS_FPA.Composition,E0,0);
toc
PS_FPA.DECS.x = data(1).x;
for i = 1:numel(E0)
    PS_FPA.Elastic.sigma_el(i) = data(i).sigma_el;
    PS_FPA.Elastic.l_el(i) = 1/(data(i).sigma_el*a0^2*PS_FPA.Density);
    PS_FPA.Elastic.l_tr(i) = 1/(data(i).sigma_tr1*a0^2*PS_FPA.Density);
    PS_FPA.DECS.y(:,i) = data(i).y;
end

%% Inelastic properties

osc.model = 'MerminLL';
osc.A = [0.045 0.013 0.011 0.053 0.168 0.093 0.061 0.018 0.066 0.045 0.021 0.01];
osc.G = [0.791 1.931 1.55 3.763 5.415 4.305 4.611 7.793 5.175 5.926 8.526 13.125];
osc.Om = [6.749 11.774 9.924 14.751 21.154 17.898 23.453 37.127 25.588 28.927 32.841 46.14];
osc.u = 1.6415029200190399;
osc.alpha = 1; 
osc.beps = 1;
osc.Ef = PS_FPA.Evb; 
osc.qtran = 0.01:0.01:20;
osc.eloss = eps:.1:110;
osc.egap = PS_FPA.Eg;

current_full_path = dbstack('-completenames');
current_file_name = dbstack;
if ispc
    ind = strfind(current_full_path(1).file,['Database\' current_file_name(1).file]);
    dirData = [current_full_path(1).file(1:ind-2) filesep 'Data\'];
elseif ismac || isunix
    ind = strfind(current_full_path(1).file,['Database/' current_file_name(1).file]);
    dirData = [current_full_path(1).file(1:ind-2) filesep 'Data/'];
end

% opt_elf = load([dirData '/opt/polystyrene.diel']);
% q = eps:0.1:10;
% eloss = eps:.5:500;
% tic
% [PMMA_FPA.ELF,~,~] = fpa_vector(q*a0,eloss/h2ev,opt_elf(:,1),opt_elf(:,4));
% toc
fpa_data = load([dirData 'elf_ps_fpa.mat']);
PS_FPA.eloss = fpa_data.omega;
PS_FPA.q = fpa_data.q;
PS_FPA.ELF = fpa_data.elf;
PS_FPA.DIIMFP = zeros(N,2,numel(E0));
PS_FPA.l_in = zeros(numel(E0),1);
for i = 1:length(E0)
    if E0(i) > 2*PS_FPA.Eg + PS_FPA.Evb
        energy = E0(i) - PS_FPA.Eg - PS_FPA.Evb;
        eloss = PS_FPA.Eg:(energy-PS_FPA.Eg)/(N-1):energy;
        if energy < 10
            omega = PS_FPA.Eg:0.01:energy;
        else
            omega = PS_FPA.Eg:0.2:energy;
        end
        PS_FPA.DIIMFP(:,1,i) = eloss;
        [iimfp, diimfp_] = diimfp_interp_fpa(E0(i),omega,PS_FPA.ELF,PS_FPA.q,PS_FPA.eloss,osc);
        PS_FPA.DIIMFP(:,2,i) = interp1(omega,diimfp_,eloss);
        ind = isnan(PS_FPA.DIIMFP(:,2,i));
        PS_FPA.DIIMFP(ind,2,i) = diimfp_(end);
        PS_FPA.l_in(i) = 1/trapz(omega/h2ev,iimfp)*a0;
    else
        PS_FPA.l_in(i) = Inf;
    end
end
