function Ag_Mermin = Make_Ag_Mermin

E0 = [1:100 110:10:200 220:20:300 350:50:500 600:100:2500 2750:250:5000 5500:500:30000];
N = 4000;

%% Basic
Ag_Mermin.Mat = 'Ag';
Ag_Mermin.M = 107.8682;
Ag_Mermin.Z = 47;
Ag_Mermin.Density = 10.5; %g/cm^3
Ag_Mermin.Density = Ag_Mermin.Density*10^-24/Ag_Mermin.M*6.022*10^23; %#/A^3
Ag_Mermin.NvTPP = 11;
Ag_Mermin.Ep = 29.8;
Ag_Mermin.Ef = 7.2;
Ag_Mermin.Wf = 4.26;
Ag_Mermin.isMetal = true;

%% Elastic properties
% {
Ag_Mermin.Elastic.x = zeros(numel(E0),1);
Ag_Mermin.Elastic.l_el = zeros(numel(E0),1);
Ag_Mermin.Elastic.l_tr = zeros(numel(E0),1);
Ag_Mermin.Elastic.x = E0;
Ag_Mermin.DECS.E0 = E0;
Ag_Mermin.Composition.Z = Ag_Mermin.Z;
Ag_Mermin.Composition.index = 1;

tic;
[data] = ElsepaRunner.RunElsepa(Ag_Mermin.Composition,E0);
toc
Ag_Mermin.DECS.x = data(1).x;
for i = 1:numel(E0)
    Ag_Mermin.Elastic.sigma_el(i) = data(i).sigma_el;
    Ag_Mermin.Elastic.l_el(i) = 1/data(i).sigma_el/Ag_Mermin.Density;
    Ag_Mermin.Elastic.l_tr(i) = 1/data(i).sigma_tr1/Ag_Mermin.Density;
    Ag_Mermin.DECS.y(:,i) = data(i).y;
end
%}

%% Inelastic properties
osc.model = 'Mermin';
osc.A = [0.05 0.19 0.05 0.09 0.11 0.04 0.07 0.22 0.07 0.14 0.04 0.09 0.01]';
osc.G = [0.16 2.83 4.1 13.99 6.45 1.31 11.06 9.79 9.91 10.06 15.08 16.95 13.19]; 
osc.Om =[3.82 7.91 11.36 56.91 17.32 4.96 48.34 26.26 41.68 35 81.47 70 91.37];
osc.alpha = 1; 
osc.beps = 1;
osc.Ef = Ag_Mermin.Ef; 
osc.qtran = 0.01:0.01:20;
osc.eloss = eps:.1:110;
osc.egap = eps;

Ag_Mermin.ELF = eps_sum_allwq(osc,'bulk');
Ag_Mermin.eloss = osc.eloss;
Ag_Mermin.q = osc.qtran;
Ag_Mermin.DIIMFP = zeros(N,2,numel(E0));
Ag_Mermin.l_in = zeros(numel(E0),1);
for i = 1:length(E0)
    if E0(i) > Ag_Mermin.Ef
        energy = E0(i) - Ag_Mermin.Ef;
        eloss = eps:(energy-eps)/(N-1):energy;
        if energy < 10
            osc.eloss = eps:0.01:energy;
        else
            osc.eloss = eps:0.2:energy;
        end
        Ag_Mermin.DIIMFP(:,1,i) = eloss;
        [iimfp, diimfp] = ndiimfp(osc,E0(i));
        Ag_Mermin.DIIMFP(:,2,i) = interp1(osc.eloss,diimfp,eloss);
        Ag_Mermin.l_in(i) = 1/trapz(osc.eloss/h2ev,iimfp)*a0;
    else
        Ag_Mermin.l_in(i) = Inf;
    end
end

%% Ionisation shells
Ag_Mermin.Shells = {'3S1/2';'3P1/2';'3P3/2';'3D3/2';'3D5/2';'4S1/2';'4P1/2';'4P3/2';'4D3/2';'4D5/2';'5S1/2';};
Ag_Mermin.EB = [724;608;577;379;373;101;69;63;11;10;8;];

