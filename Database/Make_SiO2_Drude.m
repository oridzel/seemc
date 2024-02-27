function SiO2_Drude = Make_SiO2_Drude

E0 = [1:100 110:10:200 220:20:300 350:50:500 600:100:2500 2750:250:5000 5500:500:30000];
N = 4000;

%% Basic
SiO2_Drude = struct;
SiO2_Drude.Mat = 'SiO2_Drude';
SiO2_Drude.M = 60.008;
SiO2_Drude.Z = 30;
SiO2_Drude.Density = 2.19; %g/cm^3
SiO2_Drude.Density = SiO2_Drude.Density*10^-24/SiO2_Drude.M*6.022*10^23; %#/A^3
SiO2_Drude.NvTPP = 16;
SiO2_Drude.Eg = 9.1;
SiO2_Drude.Ep = 22.02;
SiO2_Drude.Evb = 10;
SiO2_Drude.Affinity = 1.1;
SiO2_Drude.Phonon.eps_zero = 3.84; % static dielectric constant
SiO2_Drude.Phonon.eps_inf = 2.25; % high-frequency dielectric constant (the square of the static refractive index)
SiO2_Drude.Phonon.eloss = 0.1;
SiO2_Drude.isMetal = false;

%% Elastic properties
% {
SiO2_Drude.Elastic.x = zeros(numel(E0),1);
SiO2_Drude.Elastic.l_el = zeros(numel(E0),1);
SiO2_Drude.Elastic.l_tr = zeros(numel(E0),1);
SiO2_Drude.Elastic.x = E0;
SiO2_Drude.DECS.E0 = E0;
SiO2_Drude.Composition.Z = [14 8];
SiO2_Drude.Composition.index = [1 2];

tic;
[data] = ElsepaRunner.RunElsepa(SiO2_Drude.Composition,E0,0);
toc
SiO2_Drude.DECS.x = data(1).x;
for i = 1:numel(E0)
    SiO2_Drude.Elastic.sigma_el(i) = data(i).sigma_el;
    SiO2_Drude.Elastic.l_el(i) = 1/(data(i).sigma_el*a0^2*SiO2_Drude.Density);
    SiO2_Drude.Elastic.l_tr(i) = 1/(data(i).sigma_tr1*a0^2*SiO2_Drude.Density);
    SiO2_Drude.DECS.y(:,i) = data(i).y;
end
%}

%% Inelastic properties
osc.model = 'Drude';
osc.A = [18.16 25.94 39.79 42.6 91.29 134.48 186.26]';
osc.G = [1.5 2.78 3.51 3.92 8.53 30.18 79.87]; 
osc.Om =[10.58 12.13 13.72 16.9 20.61 35.17 69.21];
osc.alpha = 1; 
osc.beps = 1;
osc.Ef = SiO2_Drude.Evb; 
osc.qtran = 0.01:0.01:20;
osc.eloss = eps:.1:110;
osc.egap = SiO2_Drude.Eg;

SiO2_Drude.ELF = eps_sum_allwq(osc,'bulk');
SiO2_Drude.eloss = osc.eloss;
SiO2_Drude.q = osc.qtran;
SiO2_Drude.DIIMFP = zeros(N,2,numel(E0));
SiO2_Drude.l_in = zeros(numel(E0),1);
for i = 1:length(E0)
    if E0(i) > 2*SiO2_Drude.Eg + SiO2_Drude.Evb
        energy = E0(i) - SiO2_Drude.Eg - SiO2_Drude.Evb;
        eloss = SiO2_Drude.Eg:(energy-SiO2_Drude.Eg)/(N-1):energy;
        if energy < 10
            osc.eloss = SiO2_Drude.Eg:0.01:energy;
        else
            osc.eloss = SiO2_Drude.Eg:0.2:energy;
        end
        SiO2_Drude.DIIMFP(:,1,i) = eloss;
        [iimfp, diimfp] = ndiimfp(osc,E0(i));
        SiO2_Drude.DIIMFP(:,2,i) = interp1(osc.eloss,diimfp,eloss);
        SiO2_Drude.l_in(i) = 1/trapz(osc.eloss/h2ev,iimfp)*a0;
    else
        SiO2_Drude.l_in(i) = Inf;
    end
end

