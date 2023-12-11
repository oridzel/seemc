function Au_DL = Make_Au_DL

E0 = [1:100 150:50:500 600:100:2500 2750:250:5000 5500:500:30000];
N = 5000;

%% Basic
Au_DL.Mat = 'Au_DL';
Au_DL.M = 196.967;
Au_DL.Z = 79;
Au_DL.Density = 19.32; %g/cm^3
Au_DL.Density = Au_DL.Density*10^-24/Au_DL.M*6.022*10^23; %#/A^3
Au_DL.NvTPP = 11;
Au_DL.Ep = 29.92;
Au_DL.Ef = 9;
Au_DL.Wf = 5.2;

%% Elastic properties
% {
Au_DL.Elastic.x = zeros(numel(E0),1);
Au_DL.Elastic.l_el = zeros(numel(E0),1);
Au_DL.Elastic.l_tr = zeros(numel(E0),1);
Au_DL.Elastic.x = E0;
Au_DL.DECS.E0 = E0;
Au_DL.Composition.Z = Au_DL.Z;
Au_DL.Composition.index = 1;

tic;
[data] = ElsepaRunner.RunElsepa(Au_DL.Composition,E0);
toc
Au_DL.DECS.x = data(1).x;
for i = 1:numel(E0)
    Au_DL.Elastic.l_el(i) = 1/data(i).sigma_el/Au_DL.Density;
    Au_DL.Elastic.l_tr(i) = 1/data(i).sigma_tr1/Au_DL.Density;
    Au_DL.DECS.y(:,i) = data(i).y/trapz(data(i).x,data(i).y);
end
%}

%% Inelastic properties
osc.model = 'DrudeLindhard';
osc.A = [0.01 0.02 0.07 0.1 0.07 0.007 0.16 0.15 0.13 0.08 0.02 0.09 0.16 0.02 0.003 0.005 0.008]';
osc.G = [0.29 0.81 3.05 5.62 5.04 2.05 8.04 8.56 10.86 10.79 11.13 5.39 29.85 34.51 38.38 62.79 376.54]; 
osc.Om =[2.62 3.34 6.31 10.58 17.08 25.75 25.39 33.65 39.17 45.72 52.12 14.57 64.21 96.19 278.92 210.69 470.93];
osc.alpha = 1; 
osc.beps = 1;
osc.Ef = Au_DL.Ef; 
osc.qtran = 0.01:0.01:20;
osc.eloss = eps:.1:110;
osc.egap = eps;

Au_DL.ELF = eps_sum_allwq(osc,'bulk');
Au_DL.eloss = osc.eloss;
Au_DL.q = osc.qtran;
Au_DL.DIIMFP = zeros(N,2,numel(E0));
Au_DL.l_in = zeros(numel(E0),1);
for i = 1:length(E0)
    if E0(i) > Au_DL.Ef
        energy = E0(i) - Au_DL.Ef;
        osc.eloss = eps:(energy-eps)/(N-1):energy;
        Au_DL.DIIMFP(:,1,i) = osc.eloss;
        [iimfp, Au_DL.DIIMFP(:,2,i)] = ndiimfp(osc,E0(i));
        Au_DL.l_in(i) = 1/trapz(osc.eloss/h2ev,iimfp)*a0;
    else
        Au_DL.l_in(i) = Inf;
    end
end

%% Ionisation shells
Au_DL.Shells = {'4S1/2';'4P1/2';'4P3/2';'4D3/2';'4D5/2';'4F5/2';'4F7/2';'5S1/2';'5P1/2';'5P3/2';'5D3/2';'5D5/2';'6S1/2';};
Au_DL.EB = [763;646;549;356;339;91;87;113;76;60;10;10;9;];

