function PMMA_MLL = Make_PMMA_MLL

E0 = [1:100 110:10:200 220:20:300 350:50:500 600:100:2500 2750:250:5000 5500:500:30000];
N = 4000;

%% Basic
PMMA_MLL = struct;
PMMA_MLL.Mat = 'PMMA_MLL';
PMMA_MLL.M = 100.114;
PMMA_MLL.Z = 54;
PMMA_MLL.Density = 1.188; %g/cm^3
PMMA_MLL.Density = PMMA_MLL.Density*10^-24/PMMA_MLL.M*6.022*10^23; %#/A^3
PMMA_MLL.NvTPP = 40;
PMMA_MLL.Eg = 6.7;
PMMA_MLL.Ep = 19.85;
PMMA_MLL.Evb = 15.8;
PMMA_MLL.Affinity = 4.5;
PMMA_MLL.Phonon.eps_zero = 3.9;
PMMA_MLL.Phonon.eps_inf = 2.2;
PMMA_MLL.Phonon.eloss = 0.1;
PMMA_MLL.isMetal = false;

%% Elastic properties
PMMA_MLL.Elastic.x = zeros(numel(E0),1);
PMMA_MLL.Elastic.l_el = zeros(numel(E0),1);
PMMA_MLL.Elastic.l_tr = zeros(numel(E0),1);
PMMA_MLL.Elastic.x = E0;
PMMA_MLL.DECS.E0 = E0;
PMMA_MLL.Composition.Z = [6 8 1];
PMMA_MLL.Composition.index = [5 2 8];

tic;
[data] = ElsepaRunner.RunElsepa(PMMA_MLL.Composition,E0,0);
toc
PMMA_MLL.DECS.x = data(1).x;
for i = 1:numel(E0)
    PMMA_MLL.Elastic.sigma_el(i) = data(i).sigma_el;
    PMMA_MLL.Elastic.l_el(i) = 1/(data(i).sigma_el*a0^2*PMMA_MLL.Density);
    PMMA_MLL.Elastic.l_tr(i) = 1/(data(i).sigma_tr1*a0^2*PMMA_MLL.Density);
    PMMA_MLL.DECS.y(:,i) = data(i).y;
end

%% Inelastic properties

osc.model = 'MerminLL';
osc.A = [0.06 0.04 0.06 0.16 0.05 0.07 0.08 0.06 0.01 0.02];
osc.G = [3.61  3.87  8.88  5.81  5.57  6.31 10.03  5.51 29.96 29.82];
osc.Om = [19.88 12.40 13.24 17.34 27.12 23.24 27.66 22.36 62.55 46.31];
osc.u = 6.2212;
osc.alpha = 1; 
osc.beps = 1;
osc.Ef = PMMA_MLL.Evb; 
osc.qtran = 0.01:0.01:20;
osc.eloss = eps:.1:110;
osc.egap = PMMA_MLL.Eg;

PMMA_MLL.ELF = eps_sum_allwq(osc,'bulk');
PMMA_MLL.eloss = osc.eloss;
PMMA_MLL.q = osc.qtran;
PMMA_MLL.DIIMFP = zeros(N,2,numel(E0));
PMMA_MLL.l_in = zeros(numel(E0),1);
for i = 1:length(E0)
    if E0(i) > 2*PMMA_MLL.Eg + PMMA_MLL.Evb
        energy = E0(i) - PMMA_MLL.Eg - PMMA_MLL.Evb;
        eloss = PMMA_MLL.Eg:(energy-PMMA_MLL.Eg)/(N-1):energy;
        if energy < 10
            osc.eloss = PMMA_MLL.Eg:0.01:energy;
        else
            osc.eloss = PMMA_MLL.Eg:0.2:energy;
        end
        PMMA_MLL.DIIMFP(:,1,i) = eloss;
        [iimfp, diimfp] = ndiimfp(osc,E0(i));
        PMMA_MLL.DIIMFP(:,2,i) = interp1(osc.eloss,diimfp,eloss);
        PMMA_MLL.l_in(i) = 1/trapz(osc.eloss/h2ev,iimfp)*a0;
    else
        PMMA_MLL.l_in(i) = Inf;
    end
end

