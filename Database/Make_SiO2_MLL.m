function SiO2_MLL = Make_SiO2_MLL

E0 = [1:100 150:50:500 600:100:2500 2750:250:5000 5500:500:30000];
N = 4000;

%% Basic
SiO2_MLL = struct;
SiO2_MLL.Mat = 'SiO2_MLL';
SiO2_MLL.M = 60.008;
SiO2_MLL.Z = 30;
SiO2_MLL.Density = 2.19; %g/cm^3
SiO2_MLL.Density = SiO2_MLL.Density*10^-24/SiO2_MLL.M*6.022*10^23; %#/A^3
SiO2_MLL.NvTPP = 16;
SiO2_MLL.Eg = 9.1;
SiO2_MLL.Ep = 22.02;
SiO2_MLL.Evb = 10;
SiO2_MLL.Affinity = 1.1;
SiO2_MLL.Phonon.eps_zero = 3.84; % static dielectric constant
SiO2_MLL.Phonon.eps_inf = 2.25; % high-frequency dielectric constant (the square of the static refractive index)
SiO2_MLL.Phonon.eloss = 0.1;

%% Elastic properties
% {
SiO2_MLL.Elastic.x = zeros(numel(E0),1);
SiO2_MLL.Elastic.l_el = zeros(numel(E0),1);
SiO2_MLL.Elastic.l_tr = zeros(numel(E0),1);
SiO2_MLL.Elastic.x = E0;
SiO2_MLL.DECS.E0 = E0;
SiO2_MLL.Composition.Z = [14 8];
SiO2_MLL.Composition.index = [1 2];

tic;
[data] = ElsepaRunner.RunElsepa(SiO2_MLL.Composition,E0,0);
toc
SiO2_MLL.DECS.x = data(1).x;
for i = 1:numel(E0)
    SiO2_MLL.Elastic.l_el(i) = 1/data(i).sigma_el/SiO2_MLL.Density;
    SiO2_MLL.Elastic.l_tr(i) = 1/data(i).sigma_tr1/SiO2_MLL.Density;
    SiO2_MLL.DECS.y(:,i) = data(i).y/trapz(data(i).x,data(i).y);
end
%}

%% Inelastic properties
osc.model = 'MerminLL';
osc.A = [0.06 0.03 0.03 0.04 0.18 0.11 0.09 0.01 0.07 0.01 0.01 0.03 0.01];
osc.G = [2.00  2.00  2.00  3.98  5.27  5.89  7.62 13.62 12.52 16.41  5.26 11.13 14.06];
osc.Om = [10.00 12.00 15.00 14.27 19.23 22.72 26.24 60.45 32.31 74.19 47.94 40.89 52.92];
osc.u = 9.0747;
osc.alpha = 0; 
osc.beps = 1;
osc.Ef = SiO2_MLL.Evb; 
osc.qtran = 0.01:0.01:20;
osc.eloss = eps:.1:110;
osc.egap = SiO2_MLL.Eg;

SiO2_MLL.ELF = eps_sum_allwq(osc,'bulk');
SiO2_MLL.eloss = osc.eloss;
SiO2_MLL.q = osc.qtran;
SiO2_MLL.DIIMFP = zeros(N,2,numel(E0));
SiO2_MLL.l_in = zeros(numel(E0),1);
for i = 1:length(E0)
    if E0(i) > 2*SiO2_MLL.Eg + SiO2_MLL.Evb
        energy = E0(i) - SiO2_MLL.Eg - SiO2_MLL.Evb;
        osc.eloss = SiO2_MLL.Eg:(energy-SiO2_MLL.Eg)/(N-1):energy;
        SiO2_MLL.DIIMFP(:,1,i) = osc.eloss;
        [iimfp, diimfp] = ndiimfp(osc,E0(i));
        SiO2_MLL.DIIMFP(:,2,i) = diimfp./trapz(osc.eloss,diimfp);
        SiO2_MLL.l_in(i) = 1/trapz(osc.eloss/h2ev,iimfp)*a0;
    else
        SiO2_MLL.l_in(i) = Inf;
    end
end

