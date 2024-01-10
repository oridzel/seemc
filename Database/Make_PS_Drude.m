function PS_Drude = Make_PS_Drude

E0 = [1:100 150:50:500 600:100:2500 2750:250:5000 5500:500:30000];
N = 4000;

%% Basic
PS_Drude = struct;
PS_Drude.Mat = 'PS_Drude';
PS_Drude.M = 104.1440;
PS_Drude.Z = 56;
PS_Drude.Density = 1.05; %g/cm^3
PS_Drude.Density = PS_Drude.Density*10^-24/PS_Drude.M*6.022*10^23; %#/A^3
PS_Drude.NvTPP = 40;
PS_Drude.Eg = 4.5;
PS_Drude.Ep = 18.3;
PS_Drude.Evb = 17;
PS_Drude.Affinity = 3.5;
PS_Drude.Phonon.eps_zero = 2.5; % static dielectric constant
PS_Drude.Phonon.eps_inf = 1.01; % high-frequency dielectric constant (the square of the static refractive index)
PS_Drude.Phonon.eloss = 0.3;

%% Elastic properties
PS_Drude.Elastic.x = zeros(numel(E0),1);
PS_Drude.Elastic.l_el = zeros(numel(E0),1);
PS_Drude.Elastic.l_tr = zeros(numel(E0),1);
PS_Drude.Elastic.x = E0;
PS_Drude.DECS.E0 = E0;
PS_Drude.Composition.Z = [6 1];
PS_Drude.Composition.index = [8 8];

tic;
[data] = ElsepaRunner.RunElsepa(PS_Drude.Composition,E0,0);
toc
PS_Drude.DECS.x = data(1).x;
for i = 1:numel(E0)
    PS_Drude.Elastic.l_el(i) = 1/data(i).sigma_el/PS_Drude.Density;
    PS_Drude.Elastic.l_tr(i) = 1/data(i).sigma_tr1/PS_Drude.Density;
    PS_Drude.DECS.y(:,i) = data(i).y/trapz(data(i).x,data(i).y);
end

%% Inelastic properties
osc.model = 'Drude';
osc.A = [46.95 14.28 34.72 24.28 4.78 31.07 47.63 34.49 31.32 25.56 18.04 13.16 7.39]';
osc.G = [0.53 1.17 3.06 4.72 2.04 3.33 4.59 5 6.15 6.41 5.8 5.45 4.81]; 
osc.Om =[5.7 9.07 9.87 10.52 10.71 11.89 14.27 17.06 20.33 23.87 27.57 31.41 35.34];
osc.alpha = 1;
osc.beps = 1;
osc.Ef = PS_Drude.Evb; 
osc.qtran = 0.01:0.01:20;
osc.eloss = eps:.1:110;
osc.egap = PS_Drude.Eg;

PS_Drude.ELF = eps_sum_allwq(osc,'bulk');
PS_Drude.eloss = osc.eloss;
PS_Drude.q = osc.qtran;
PS_Drude.DIIMFP = zeros(N,2,numel(E0));
PS_Drude.l_in = zeros(numel(E0),1);
for i = 1:length(E0)
    if E0(i) > 2*PS_Drude.Eg + PS_Drude.Evb
        energy = E0(i) - PS_Drude.Eg - PS_Drude.Evb;
        osc.eloss = PS_Drude.Eg:(energy-PS_Drude.Eg)/(N-1):energy;
        PS_Drude.DIIMFP(:,1,i) = osc.eloss;
        [iimfp, PS_Drude.DIIMFP(:,2,i)] = ndiimfp(osc,E0(i));
        PS_Drude.l_in(i) = 1/trapz(osc.eloss/h2ev,iimfp)*a0;
    else
        PS_Drude.l_in(i) = Inf;
    end
end

