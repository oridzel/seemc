function W = Make_W

E0 = [1:100 150:50:500 600:100:2500 2750:250:5000 5500:500:30000];
N = 4000;

%% Basic
W.Mat = 'W';
W.M = 183.85;
W.Z = 74;
W.Density = 19.3; %g/cm^3
W.Density = W.Density*10^-24/W.M*6.022*10^23; %#/A^3
W.NvTPP = 6;
W.Ep = 22.86;
W.Ef = 10.1;
W.Wf = 4.55;

%% Elastic properties
% {
W.Elastic.x = zeros(numel(E0),1);
W.Elastic.l_el = zeros(numel(E0),1);
W.Elastic.l_tr = zeros(numel(E0),1);
W.Elastic.x = E0;
W.DECS.E0 = E0;
W.Composition.Z = W.Z;
W.Composition.index = 1;

tic;
[data] = ElsepaRunner.RunElsepa(W.Composition,E0);
toc
W.DECS.x = data(1).x;
for i = 1:numel(E0)
    W.Elastic.l_el(i) = 1/data(i).sigma_el/W.Density;
    W.Elastic.l_tr(i) = 1/data(i).sigma_tr1/W.Density;
    W.DECS.y(:,i) = data(i).y/trapz(data(i).x,data(i).y);
end
%}

%% Inelastic properties
osc.model = 'Mermin';
osc.A = [0.01 0.006 0.006 0.003 0.007 0.04 0.06 0.003 0.09 0.45 0.04 0.11 0.25 0.04]';
osc.G = [0.53 0.76 0.81 0.57 1.03 2.33 4.04 0.14 8.93 6.58 6.55 8.14 25.15 194.48]; 
osc.Om =[1.37 2.14 2.86 3.92 8.66 10.07 14.95 0.87 31.52 25.48 37.63 42.68 56.47 188.48];
osc.alpha = 1; 
osc.beps = 1;
osc.Ef = W.Ef; 
osc.qtran = 0.01:0.01:20;
osc.eloss = eps:.1:110;
osc.egap = eps;

W.ELF = eps_sum_allwq(osc,'bulk');
W.eloss = osc.eloss;
W.q = osc.qtran;
W.DIIMFP = zeros(N,2,numel(E0));
W.l_in = zeros(numel(E0),1);
for i = 1:length(E0)
    if E0(i) > W.Ef
        energy = E0(i) - W.Ef;
        osc.eloss = eps:(energy-eps)/(N-1):energy;
        W.DIIMFP(:,1,i) = osc.eloss;
        [iimfp, diimfp] = ndiimfp(osc,E0(i));
        W.DIIMFP(:,2,i) = diimfp./trapz(osc.eloss,diimfp);
        W.l_in(i) = 1/trapz(osc.eloss/h2ev,iimfp)*a0;
    else
        W.l_in(i) = Inf;
    end
end

%% Ionisation shells
W.Shells = {'3D3/2','3D5/2', '4S1/2', '4P1/2', '4P3/2', '4D3/2', '4D5/2', '4F5/2', '4F7/2', '5S1/2', '5P1/2', '5P3/2', '5D3/2', '6S1/2'};
W.EB = [1874; 1812; 599;495;428;261;248;39;36;80;51;41;9;8];

