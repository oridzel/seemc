function PMMA = Make_PMMA

E0 = [1:100 150:50:500 600:100:2500 2750:250:5000 5500:500:30000];
N = 5000;

%% Basic
PMMA = struct;
PMMA.Mat = 'PMMA';
PMMA.M = 100.114;
PMMA.Z = 54;
PMMA.Density = 1.188; %g/cm^3
PMMA.Density = PMMA.Density*10^-24/PMMA.M*6.022*10^23; %#/A^3
PMMA.NvTPP = 40;
PMMA.Eg = 6.7;
PMMA.Ep = 19.85;
PMMA.Evb = 15.8;
PMMA.Affinity = 4.3;

%% Elastic properties
PMMA.Elastic.x = zeros(numel(E0),1);
PMMA.Elastic.l_el = zeros(numel(E0),1);
PMMA.Elastic.l_tr = zeros(numel(E0),1);
PMMA.Elastic.x = E0;
PMMA.DECS.E0 = E0;
PMMA.Composition.Z = [6 8 1];
PMMA.Composition.index = [5 2 8];

tic;
[data] = ElsepaRunner.RunElsepa(PMMA.Composition,E0);
toc
PMMA.DECS.x = data(1).x;
for i = 1:numel(E0)
    PMMA.Elastic.l_el(i) = 1/data(i).sigma_el/PMMA.Density;
    PMMA.Elastic.l_tr(i) = 1/data(i).sigma_tr1/PMMA.Density;
    PMMA.DECS.y(:,i) = data(i).y/trapz(data(i).x,data(i).y);
end

%% Inelastic properties
osc.model = 'Drude';
osc.A = [11.96 25.16 22.138 26.7 27.7 29.5 27.6 23.5 39.2 28.4 28.5 22.4 16.76 19.9 41.95]';
osc.G = [2.1 2.2 2.1 2.6 3.1 4 5.1 6.1 9.3 11.6 15.8 16.6 21.2 22.3 41.6]; 
osc.Om =[9.3 10.5 11.9 13.2 14.6 16.5 18.7 21.2 24.2 28.5 32 36.3 44.6 48.4 62];
osc.alpha = 1; 
osc.beps = 1;
osc.Ef = PMMA.Evb; 
osc.qtran = 0.01:0.01:20;
osc.eloss = eps:.1:110;
osc.egap = PMMA.Eg;

PMMA.ELF = eps_sum_allwq(osc,'bulk');
PMMA.eloss = osc.eloss;
PMMA.q = osc.qtran;
PMMA.DIIMFP = zeros(N,2,numel(E0));
PMMA.l_in = zeros(numel(E0),1);
for i = 1:length(E0)
    energy = E0(i) + PMMA.Eg + PMMA.Evb;
    if energy > 2*PMMA.Eg + PMMA.Evb
        osc.eloss = eps:(energy-eps)/(N-1):energy;
        PMMA.DIIMFP(:,1,i) = osc.eloss;
        [iimfp, PMMA.DIIMFP(:,2,i)] = ndiimfp(osc,energy);
        eloss_interp = PMMA.Eg:(energy-2*PMMA.Eg-PMMA.Evb)/N:energy-PMMA.Eg-PMMA.Evb;
        iimfp_interp = interp1(osc.eloss,iimfp,eloss_interp);
        PMMA.l_in(i) = 1/trapz(eloss_interp/h2ev,iimfp_interp)*a0;
    else
        PMMA.l_in(i) = Inf;
    end
end

