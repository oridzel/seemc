function Au = Make_Au

E0 = [1:100 150:50:500 600:100:2500 2750:250:5000 5500:500:30000];
N = 5000;

%% Basic
Au.Mat = 'Au';
Au.M = 196.967;
Au.Z = 79;
Au.Density = 19.32; %g/cm^3
Au.Density = Au.Density*10^-24/Au.M*6.022*10^23; %#/A^3
Au.NvTPP = 11;
Au.Ep = 29.92;
Au.Ef = 9;
Au.Wf = 5.2;

%% Elastic properties
% {
Au.Elastic.x = zeros(numel(E0),1);
Au.Elastic.l_el = zeros(numel(E0),1);
Au.Elastic.l_tr = zeros(numel(E0),1);
Au.Elastic.x = E0;
Au.DECS.E0 = E0;
Au.Composition.Z = Au.Z;
Au.Composition.index = 1;

tic;
[data] = ElsepaRunner.RunElsepa(Au.Composition,E0);
toc
Au.DECS.x = data(1).x;
for i = 1:numel(E0)
    Au.Elastic.l_el(i) = 1/data(i).sigma_el/Au.Density;
    Au.Elastic.l_tr(i) = 1/data(i).sigma_tr1/Au.Density;
    Au.DECS.y(:,i) = data(i).y/trapz(data(i).x,data(i).y);
end
%}

%% Inelastic properties
osc.model = 'Mermin';
% osc.A = [0.01 0.02 0.07 0.1 0.07 0.007 0.16 0.15 0.13 0.08 0.02 0.09 0.16 0.02 0.003 0.005 0.008]';
% osc.G = [0.29 0.81 3.05 5.62 5.04 2.05 8.04 8.56 10.86 10.79 11.13 5.39 29.85 34.51 38.38 62.79 376.54]; 
% osc.Om =[2.62 3.34 6.31 10.58 17.08 25.75 25.39 33.65 39.17 45.72 52.12 14.57 64.21 96.19 278.92 210.69 470.93];

osc.A = [0.0312 0.0429 0.1253 0.0469 0.2512 0.0792 0.1284 0.0517 0.0459 0.0515 0.1158 0.0300];
osc.G = [0.9698 3.2772 6.4314 3.5576 8.4554 6.9578 12.2837 12.6974 14.2658 4.4995 24.8269 25.0628];
osc.Om = [3.2380 6.2638 11.6974 17.5575 26.2625 33.3481 38.4167 45.6709 53.1109 15.2746 64.7666 84.5982];
osc.alpha = 1; 
osc.beps = 1;
osc.Ef = Au.Ef; 
osc.qtran = 0.01:0.01:20;
osc.eloss = eps:.1:110;
osc.egap = eps;

Au.ELF = eps_sum_allwq(osc,'bulk');
Au.eloss = osc.eloss;
Au.q = osc.qtran;
Au.DIIMFP = zeros(N,2,numel(E0));
Au.l_in = zeros(numel(E0),1);
for i = 1:length(E0)
    if E0(i) > Au.Ef
        energy = E0(i) - Au.Ef;
        osc.eloss = eps:(energy-eps)/(N-1):energy;
        Au.DIIMFP(:,1,i) = osc.eloss;
        [iimfp, diimfp] = ndiimfp(osc,E0(i));
        Au.DIIMFP(:,2,i) = diimfp./trapz(osc.eloss,diimfp);
        Au.l_in(i) = 1/trapz(osc.eloss/h2ev,iimfp)*a0;
    else
        Au.l_in(i) = Inf;
    end
end

%% Ionisation shells
Au.Shells = {'4S1/2';'4P1/2';'4P3/2';'4D3/2';'4D5/2';'4F5/2';'4F7/2';'5S1/2';'5P1/2';'5P3/2';'5D3/2';'5D5/2';'6S1/2';};
Au.EB = [763;646;549;356;339;91;87;113;76;60;10;10;9;];

