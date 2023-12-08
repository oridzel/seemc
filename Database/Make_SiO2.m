function SiO2 = Make_SiO2

E0 = [1:100 150:50:500 600:100:2500 2750:250:5000 5500:500:30000];
N = 4000;

%% Basic
SiO2 = struct;
SiO2.Mat = 'SiO2';
SiO2.M = 60.008;
SiO2.Z = 30;
SiO2.Density = 2.19; %g/cm^3
SiO2.Density = SiO2.Density*10^-24/SiO2.M*6.022*10^23; %#/A^3
SiO2.NvTPP = 16;
SiO2.Eg = 9.1;
SiO2.Ep = 22.02;
SiO2.Evb = 10;
SiO2.Affinity = 1.1;
SiO2.Phonon.eps_zero = 3.84;
SiO2.Phonon.eps_inf = 2.25;
SiO2.Phonon.eloss = 0.1;

%% Elastic properties
% {
SiO2.Elastic.x = zeros(numel(E0),1);
SiO2.Elastic.l_el = zeros(numel(E0),1);
SiO2.Elastic.l_tr = zeros(numel(E0),1);
SiO2.Elastic.x = E0;
SiO2.DECS.E0 = E0;
SiO2.Composition.Z = [14 8];
SiO2.Composition.index = [1 2];

tic;
[data] = ElsepaRunner.RunElsepa(SiO2.Composition,E0);
toc
SiO2.DECS.x = data(1).x;
for i = 1:numel(E0)
    SiO2.Elastic.l_el(i) = 1/data(i).sigma_el/SiO2.Density;
    SiO2.Elastic.l_tr(i) = 1/data(i).sigma_tr1/SiO2.Density;
    SiO2.DECS.y(:,i) = data(i).y/trapz(data(i).x,data(i).y);
end
%}

%% Inelastic properties
osc.model = 'Drude';
osc.A = [18.16 25.94 39.79 42.6 91.29 134.48 186.26]';
osc.G = [1.5 2.78 3.51 3.92 8.53 30.18 79.87]; 
osc.Om =[10.58 12.13 13.72 16.9 20.61 35.17 69.21];
osc.alpha = 1; 
osc.beps = 1;
osc.Ef = SiO2.Evb; 
osc.qtran = 0.01:0.01:20;
osc.eloss = eps:.1:110;
osc.egap = SiO2.Eg;

SiO2.ELF = eps_sum_allwq(osc,'bulk');
SiO2.eloss = osc.eloss;
SiO2.q = osc.qtran;
SiO2.DIIMFP = zeros(N,2,numel(E0));
SiO2.l_in = zeros(numel(E0),1);
for i = 1:length(E0)
    if E0(i) > 2*SiO2.Eg + SiO2.Evb
        energy = E0(i) - SiO2.Eg - SiO2.Evb;
        osc.eloss = SiO2.Eg:(energy-SiO2.Eg)/(N-1):energy;
        SiO2.DIIMFP(:,1,i) = osc.eloss;
        [iimfp, SiO2.DIIMFP(:,2,i)] = ndiimfp(osc,E0(i));
        SiO2.l_in(i) = 1/trapz(osc.eloss/h2ev,iimfp)*a0;
    else
        SiO2.l_in(i) = Inf;
    end
end

