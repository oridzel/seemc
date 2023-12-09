function Cu = Make_Cu

E0 = [1:100 150:50:500 600:100:2500 2750:250:5000 5500:500:30000];
N = 4000;

%% Basic
Cu.Mat = 'Cu';
Cu.M = 63.546;
Cu.Z = 29;
Cu.Density = 8.96; %g/cm^3
Cu.Density = Cu.Density*10^-24/Cu.M*6.022*10^23; %#/A^3
Cu.NvTPP = 11;
Cu.Ep = 35.87;
Cu.Ef = 8.7;
Cu.Wf = 4.65;

%% Elastic properties
% {
Cu.Elastic.x = zeros(numel(E0),1);
Cu.Elastic.l_el = zeros(numel(E0),1);
Cu.Elastic.l_tr = zeros(numel(E0),1);
Cu.Elastic.x = E0;
Cu.DECS.E0 = E0;
Cu.Composition.Z = Cu.Z;
Cu.Composition.index = 1;

tic;
[data] = ElsepaRunner.RunElsepa(Cu.Composition,E0);
toc
Cu.DECS.x = data(1).x;
for i = 1:numel(E0)
    Cu.Elastic.l_el(i) = 1/data(i).sigma_el/Cu.Density;
    Cu.Elastic.l_tr(i) = 1/data(i).sigma_tr1/Cu.Density;
    Cu.DECS.y(:,i) = data(i).y/trapz(data(i).x,data(i).y);
end
%}

%% Inelastic properties
osc.model = 'Mermin';
osc.A = [0.01 0.02 0.06 0.11 0.38 0.15 0.07 0.05 0.08 0.1]';
osc.G = [0.69 0.8 2.36 6.59 10.82 4.89 40.77 19.15 14.12 161.76]; 
osc.Om =[3.46 4.2 7.54 28.34 19.39 10.48 66.56 47.52 36.13 126.04];
osc.alpha = 1; 
osc.beps = 1;
osc.Ef = Cu.Ef; 
osc.qtran = 0.01:0.01:20;
osc.eloss = eps:.1:110;
osc.egap = eps;

Cu.ELF = eps_sum_allwq(osc,'bulk');
Cu.eloss = osc.eloss;
Cu.q = osc.qtran;
Cu.DIIMFP = zeros(N,2,numel(E0));
Cu.l_in = zeros(numel(E0),1);
for i = 1:length(E0)
    if E0(i) > Cu.Ef
        energy = E0(i) - Cu.Ef;
        osc.eloss = eps:(energy-eps)/(N-1):energy;
        Cu.DIIMFP(:,1,i) = osc.eloss;
        [iimfp, Cu.DIIMFP(:,2,i)] = ndiimfp(osc,E0(i));
        Cu.l_in(i) = 1/trapz(osc.eloss/h2ev,iimfp)*a0;
    else
        Cu.l_in(i) = Inf;
    end
end

%% Ionisation shells
Cu.Shells = {'2S1/2';'2P1/2';'2P3/2';'3S1/2';'3P1/2';'3P3/2';'3D3/2';'3D5/2';'4S1/2';};
Cu.EB = [1103;958;938;127;82;80;11;10;8;];

