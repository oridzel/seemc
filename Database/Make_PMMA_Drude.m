function PMMA_Drude = Make_PMMA_Drude

E0 = [1:100 150:50:500 600:100:2500 2750:250:5000 5500:500:30000];
N = 4000;

%% Basic
PMMA_Drude = struct;
PMMA_Drude.Mat = 'PMMA_Drude';
PMMA_Drude.M = 100.114;
PMMA_Drude.Z = 54;
PMMA_Drude.Density = 1.188; %g/cm^3
PMMA_Drude.Density = PMMA_Drude.Density*10^-24/PMMA_Drude.M*6.022*10^23; %#/A^3
PMMA_Drude.NvTPP = 40;
PMMA_Drude.Eg = 6.7;
PMMA_Drude.Ep = 19.85;
PMMA_Drude.Evb = 15.8;
PMMA_Drude.Affinity = 4.5;
PMMA_Drude.Phonon.eps_zero = 3.9;
PMMA_Drude.Phonon.eps_inf = 2.2;
PMMA_Drude.Phonon.eloss = 0.1;

%% Elastic properties
PMMA_Drude.Elastic.x = zeros(numel(E0),1);
PMMA_Drude.Elastic.l_el = zeros(numel(E0),1);
PMMA_Drude.Elastic.l_tr = zeros(numel(E0),1);
PMMA_Drude.Elastic.x = E0;
PMMA_Drude.DECS.E0 = E0;
PMMA_Drude.Composition.Z = [6 8 1];
PMMA_Drude.Composition.index = [5 2 8];

tic;
[data] = ElsepaRunner.RunElsepa(PMMA_Drude.Composition,E0,0);
toc
PMMA_Drude.DECS.x = data(1).x;
for i = 1:numel(E0)
    PMMA_Drude.Elastic.l_el(i) = 1/data(i).sigma_el/PMMA_Drude.Density;
    PMMA_Drude.Elastic.l_tr(i) = 1/data(i).sigma_tr1/PMMA_Drude.Density;
    PMMA_Drude.DECS.y(:,i) = data(i).y/trapz(data(i).x,data(i).y);
end

%% Inelastic properties
osc.model = 'Drude';
osc.A = [11.96 25.16 22.138 26.7 27.7 29.5 27.6 23.5 39.2 28.4 28.5 22.4 16.76 19.9 41.95]';
osc.G = [2.1 2.2 2.1 2.6 3.1 4 5.1 6.1 9.3 11.6 15.8 16.6 21.2 22.3 41.6]; 
osc.Om =[9.3 10.5 11.9 13.2 14.6 16.5 18.7 21.2 24.2 28.5 32 36.3 44.6 48.4 62];
osc.alpha = 1; 
osc.beps = 1;
osc.Ef = PMMA_Drude.Evb; 
osc.qtran = 0.01:0.01:20;
osc.eloss = eps:.1:110;
osc.egap = PMMA_Drude.Eg;

PMMA_Drude.ELF = eps_sum_allwq(osc,'bulk');
PMMA_Drude.eloss = osc.eloss;
PMMA_Drude.q = osc.qtran;
PMMA_Drude.DIIMFP = zeros(N,2,numel(E0));
PMMA_Drude.l_in = zeros(numel(E0),1);
for i = 1:length(E0)
    if E0(i) > 2*PMMA_Drude.Eg + PMMA_Drude.Evb
        energy = E0(i) - PMMA_Drude.Eg - PMMA_Drude.Evb;
        osc.eloss = PMMA_Drude.Eg:(energy-PMMA_Drude.Eg)/(N-1):energy;
        PMMA_Drude.DIIMFP(:,1,i) = osc.eloss;
        [iimfp, diimfp] = ndiimfp(osc,E0(i));
        PMMA_Drude.DIIMFP(:,2,i) = diimfp./trapz(osc.eloss,diimfp);
        PMMA_Drude.l_in(i) = 1/trapz(osc.eloss/h2ev,iimfp)*a0;
    else
        PMMA_Drude.l_in(i) = Inf;
    end
end

