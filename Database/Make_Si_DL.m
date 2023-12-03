function Si_DL = Make_Si_DL

E0 = [1:100 150:50:500 600:100:2500 2750:250:5000 5500:500:30000];
N = 2000;

%% Basic
Si_DL.Mat = 'Si';
Si_DL.M = 28.0855;
Si_DL.Z = 14;
Si_DL.Density = 2.33; %g/cm^3
Si_DL.Density = Si_DL.Density*10^-24/Si_DL.M*6.022*10^23; %#/A^3
Si_DL.NvTPP = 4;
Si_DL.NvSGS = 4;
Si_DL.Eg = 1.1;
Si_DL.Ep = 16.59;
Si_DL.Ef = 12.5;
Si_DL.Evb = 12;
Si_DL.Affinity = 4.05;

%% Elastic properties
Si_DL.Elastic.x = zeros(numel(E0),1);
Si_DL.Elastic.l_el = zeros(numel(E0),1);
Si_DL.Elastic.l_tr = zeros(numel(E0),1);
Si_DL.Elastic.x = E0;
Si_DL.DECS.E0 = E0;
Si_DL.Composition.Z = 14;
Si_DL.Composition.index = 1;

tic;
[data] = ElsepaRunner.RunElsepa(Si_DL.Composition,E0);
toc
Si_DL.DECS.x = data(1).x;
for i = 1:numel(E0)
    Si_DL.Elastic.l_el(i) = 1/data(i).sigma_el/Si_DL.Density;
    Si_DL.Elastic.l_tr(i) = 1/data(i).sigma_tr1/Si_DL.Density;
    Si_DL.DECS.y(:,i) = data(i).y/trapz(data(i).x,data(i).y);
end

%% Inelastic properties
osc.model = 'DrudeLindhard';
osc.A = [0.05 0.1 0.4 0.045 0.04 0.29 0.011 0.006 0.002]';
osc.G = [5 3 2.4 7 2 2.5 40 80 150]; 
osc.Om =[20 14 17.2 25 13.2 15.7 130 180 280];
osc.alpha = 1; 
osc.beps = 1;
osc.Ef = Si_DL.Ef; 
osc.qtran = 0.01:0.01:20;
osc.eloss = eps:.1:110;
osc.egap = Si_DL.Eg;

Si_DL.ELF = eps_sum_allwq(osc,'bulk');
Si_DL.eloss = osc.eloss;
Si_DL.q = osc.qtran;
Si_DL.DIIMFP = zeros(N,2,numel(E0));
Si_DL.l_in = zeros(numel(E0),1);
for i = 1:length(E0)
    energy = E0(i) + Si_DL.Eg + Si_DL.Evb;
    if energy > 2*Si_DL.Eg + Si_DL.Evb
        osc.eloss = eps:(energy-eps)/(N-1):energy;
        Si_DL.DIIMFP(:,1,i) = osc.eloss;
        [iimfp, Si_DL.DIIMFP(:,2,i)] = ndiimfp(osc,energy);
        eloss_interp = Si_DL.Eg:(energy-2*Si_DL.Eg-Si_DL.Evb)/N:energy-Si_DL.Eg-Si_DL.Evb;
        iimfp_interp = interp1(osc.eloss,iimfp,eloss_interp);
        Si_DL.l_in(i) = 1/trapz(eloss_interp/h2ev,iimfp_interp)*a0;
    else
        Si_DL.l_in(i) = Inf;
    end
end

%% Ionisation shells
Si_DL.Shells = {'2S1/2';'2P1/2';'2P3/2';'3S1/2';'3P1/2';};
Si_DL.EB = [154;104;104;13;8;];

