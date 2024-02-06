function Si_Mermin = Make_Si_Mermin

E0 = [1:100 110:10:200 220:20:300 350:50:500 600:100:2500 2750:250:5000 5500:500:30000];
N = 4000;

%% Basic
Si_Mermin.Mat = 'Si_Mermin';
Si_Mermin.M = 28.0855;
Si_Mermin.Z = 14;
Si_Mermin.Density = 2.33; %g/cm^3
Si_Mermin.Density = Si_Mermin.Density*10^-24/Si_Mermin.M*6.022*10^23; %#/A^3
Si_Mermin.NvTPP = 4;
Si_Mermin.NvSGS = 4;
Si_Mermin.Eg = 1.1;
Si_Mermin.Ep = 16.59;
Si_Mermin.Ef = 12.5;
Si_Mermin.Evb = 12.44;
Si_Mermin.Affinity = 4.05;
Si_Mermin.isMetal = false;

%% Elastic properties
Si_Mermin.Elastic.x = zeros(numel(E0),1);
Si_Mermin.Elastic.l_el = zeros(numel(E0),1);
Si_Mermin.Elastic.l_tr = zeros(numel(E0),1);
Si_Mermin.Elastic.x = E0;
Si_Mermin.DECS.E0 = E0;
Si_Mermin.Composition.Z = Si_Mermin.Z;
Si_Mermin.Composition.index = 1;

tic;
[data] = ElsepaRunner.RunElsepa(Si_Mermin.Composition,E0);
toc
Si_Mermin.DECS.x = data(1).x;
for i = 1:numel(E0)
    Si_Mermin.Elastic.sigma_el(i) = data(i).sigma_el;
    Si_Mermin.Elastic.l_el(i) = 1/data(i).sigma_el/Si_Mermin.Density;
    Si_Mermin.Elastic.l_tr(i) = 1/data(i).sigma_tr1/Si_Mermin.Density;
    Si_Mermin.DECS.y(:,i) = data(i).y;
end

%% Inelastic properties
osc.model = 'Mermin';
osc.A = [0.05 0.1 0.4 0.045 0.04 0.29 0.011 0.006 0.002]';
osc.G = [5 3 2.4 7 2 2.5 40 80 150]; 
osc.Om =[20 14 17.2 25 13.2 15.7 130 180 280];
osc.alpha = 1; 
osc.beps = 1;
osc.Ef = Si_Mermin.Ef; 
osc.qtran = 0.01:0.01:20;
osc.eloss = eps:.1:110;
osc.egap = Si_Mermin.Eg;

Si_Mermin.ELF = eps_sum_allwq(osc,'bulk');
Si_Mermin.eloss = osc.eloss;
Si_Mermin.q = osc.qtran;
Si_Mermin.DIIMFP = zeros(N,2,numel(E0));
Si_Mermin.l_in = zeros(numel(E0),1);
for i = 1:length(E0)
    if E0(i) > 2*Si_Mermin.Eg + Si_Mermin.Evb
        energy = E0(i) - Si_Mermin.Eg - Si_Mermin.Evb;
        eloss = Si_Mermin.Eg:(energy-Si_Mermin.Eg)/(N-1):energy;
        if energy < 10
            osc.eloss = Si_Mermin.Eg:0.01:energy;
        else
            osc.eloss = Si_Mermin.Eg:0.2:energy;
        end
        Si_Mermin.DIIMFP(:,1,i) = eloss;
        [iimfp, diimfp] = ndiimfp(osc,E0(i));
        Si_Mermin.DIIMFP(:,2,i) = interp1(osc.eloss,diimfp,eloss);
        Si_Mermin.l_in(i) = 1/trapz(osc.eloss/h2ev,iimfp)*a0;
    else
        Si_Mermin.l_in(i) = Inf;
    end
end

%% Ionisation shells
Si_Mermin.Shells = {'2S1/2';'2P1/2';'2P3/2';'3S1/2';'3P1/2';};
Si_Mermin.EB = [154;104;104;13;8;];

