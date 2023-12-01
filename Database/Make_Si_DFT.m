function Si_DFT = Make_Si_DFT

E0 = [1:100 150:50:500 600:100:2500 2750:250:5000 5500:500:30000];
N = 2000;

%% Basic
Si_DFT.Mat = 'Si';
Si_DFT.M = 28.0855;
Si_DFT.Z = 14;
Si_DFT.Density = 2.33; %g/cm^3
Si_DFT.Density = Si_DFT.Density*10^-24/Si_DFT.M*6.022*10^23; %#/A^3
Si_DFT.NvTPP = 4;
Si_DFT.NvSGS = 4;
Si_DFT.Eg = 1.1;
Si_DFT.Ep = 16.59;
Si_DFT.Ef = 12.5;
Si_DFT.Evb = 12;
Si_DFT.Affinity = 4.05;

%% Elastic properties
Si_DFT.Elastic.x = zeros(numel(E0),1);
Si_DFT.Elastic.l_el = zeros(numel(E0),1);
Si_DFT.Elastic.l_tr = zeros(numel(E0),1);
Si_DFT.Elastic.x = E0;
Si_DFT.DECS.E0 = E0;
Si_DFT.Composition.Z = 14;
Si_DFT.Composition.index = 1;

tic;
[data] = ElsepaRunner.RunElsepa(Si_DFT.Composition,E0);
toc
Si_DFT.DECS.x = data(1).x;
for i = 1:numel(E0)
    Si_DFT.Elastic.l_el(i) = 1/data(i).sigma_el/Si_DFT.Density;
    Si_DFT.Elastic.l_tr(i) = 1/data(i).sigma_tr1/Si_DFT.Density;
    Si_DFT.DECS.y(:,i) = data(i).y/trapz(data(i).x,data(i).y);
end

%% Inelastic properties
osc.model = 'Mermin';
osc.A = [0.05 0.1 0.4 0.045 0.04 0.29 0.011 0.006 0.002]';
osc.G = [5 3 2.4 7 2 2.5 40 80 150]; 
osc.Om =[20 14 17.2 25 13.2 15.7 130 180 280];
osc.alpha = 1; 
osc.beps = 1;
osc.Ef = Si_DFT.Ef; 
osc.qtran = 0.01:0.01:20;
osc.eloss = eps:.1:110;
osc.egap = Si_DFT.Eg;

load 'C:/Users/onr5/OneDrive - NIST/dev/m-scripts/elf_si_dft_ocean_0_01_b1l1.mat';
Si_DFT.ELF = elf;
Si_DFT.eloss = omega;
Si_DFT.q = q;

Si_DFT.DIIMFP = zeros(N,2,numel(E0));
Si_DFT.l_in = zeros(numel(E0),1);
for i = 1:length(E0)
    energy = E0(i) + Si_DFT.Eg + Si_DFT.Evb;
    if energy > 2*Si_DFT.Eg + Si_DFT.Evb
        osc.eloss = eps:(energy-eps)/(N-1):energy;
        Si_DFT.DIIMFP(:,1,i) = osc.eloss;
        if energy < 110
            [iimfp, Si_DFT.DIIMFP(:,2,i)] = ndiimfp(osc,energy,elf,q,omega);
        else
            [iimfp, Si_DFT.DIIMFP(:,2,i)] = ndiimfp(osc,energy);
        end
        eloss_interp = Si_DFT.Eg:(energy-2*Si_DFT.Eg-Si_DFT.Evb)/N:energy-Si_DFT.Eg-Si_DFT.Evb;
        iimfp_interp = interp1(osc.eloss,iimfp,eloss_interp);
        Si_DFT.l_in(i) = 1/trapz(eloss_interp/h2ev,iimfp_interp)*a0;
    else
        Si_DFT.l_in(i) = Inf;
    end
end

%% Ionisation shells
Si_DFT.Shells = {'2S1/2';'2P1/2';'2P3/2';'3S1/2';'3P1/2';};
Si_DFT.EB = [154;104;104;13;8;];

