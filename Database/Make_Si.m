function Si = Make_Si

E0 = [1:100 150:50:500 600:100:2500 2750:250:5000 5500:500:30000];
N = 5000;

%% Basic
Si.Mat = 'Si';
Si.M = 28.0855;
Si.Z = 14;
Si.Density = 2.33; %g/cm^3
Si.Density = Si.Density*10^-24/Si.M*6.022*10^23; %#/A^3
Si.NvTPP = 4;
Si.NvSGS = 4;
Si.Eg = 1.1;
Si.Ep = 16.59;
Si.Ef = 12.5;
Si.Evb = 12.44;
Si.Affinity = 4.05;

%% Elastic properties
Si.Elastic.x = zeros(numel(E0),1);
Si.Elastic.l_el = zeros(numel(E0),1);
Si.Elastic.l_tr = zeros(numel(E0),1);
Si.Elastic.x = E0;
Si.DECS.E0 = E0;
Si.Composition.Z = 14;
Si.Composition.index = 1;

tic;
[data] = ElsepaRunner.RunElsepa(Si.Composition,E0);
toc
Si.DECS.x = data(1).x;
for i = 1:numel(E0)
    Si.Elastic.l_el(i) = 1/data(i).sigma_el/Si.Density;
    Si.Elastic.l_tr(i) = 1/data(i).sigma_tr1/Si.Density;
    Si.DECS.y(:,i) = data(i).y/trapz(data(i).x,data(i).y);
end

%% Inelastic properties
osc.model = 'Mermin';
osc.A = [0.05 0.1 0.4 0.045 0.04 0.29 0.011 0.006 0.002]';
osc.G = [5 3 2.4 7 2 2.5 40 80 150]; 
osc.Om =[20 14 17.2 25 13.2 15.7 130 180 280];
osc.alpha = 1; 
osc.beps = 1;
osc.Ef = Si.Ef; 
osc.qtran = 0.01:0.01:20;
osc.eloss = eps:.1:110;
osc.egap = Si.Eg;

Si.ELF = eps_sum_allwq(osc,'bulk');
Si.eloss = osc.eloss;
Si.q = osc.qtran;
Si.DIIMFP = zeros(N,2,numel(E0));
Si.l_in = zeros(numel(E0),1);
for i = 1:length(E0)
    if E0(i) > 2*Si.Eg + Si.Evb
        energy = E0(i) - Si.Eg - Si.Evb;
        osc.eloss = Si.Eg:(energy-Si.Eg)/(N-1):energy;
        Si.DIIMFP(:,1,i) = osc.eloss;
        [iimfp, diimfp] = ndiimfp(osc,E0(i));
        Si.DIIMFP(:,2,i) = diimfp./trapz(osc.eloss,diimfp);
        Si.l_in(i) = 1/trapz(osc.eloss/h2ev,iimfp)*a0;
    else
        Si.l_in(i) = Inf;
    end
end

%% Ionisation shells
Si.Shells = {'2S1/2';'2P1/2';'2P3/2';'3S1/2';'3P1/2';};
Si.EB = [154;104;104;13;8;];

