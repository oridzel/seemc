function Cu_Mermin = Make_Cu_Mermin

E0 = [1:100 110:10:200 220:20:300 350:50:500 600:100:2500 2750:250:5000 5500:500:30000];
N = 4000;

%% Basic
Cu_Mermin.Mat = 'Cu_Mermin';
Cu_Mermin.M = 63.546;
Cu_Mermin.Z = 29;
Cu_Mermin.Density = 8.96; %g/cm^3
Cu_Mermin.Density = Cu_Mermin.Density*10^-24/Cu_Mermin.M*6.022*10^23; %#/A^3
Cu_Mermin.NvTPP = 11;
Cu_Mermin.Ep = 35.87;
Cu_Mermin.Ef = 8.7;
Cu_Mermin.Wf = 4.65;
Cu_Mermin.isMetal = true;

%% Elastic properties
% {
Cu_Mermin.Elastic.x = zeros(numel(E0),1);
Cu_Mermin.Elastic.l_el = zeros(numel(E0),1);
Cu_Mermin.Elastic.l_tr = zeros(numel(E0),1);
Cu_Mermin.Elastic.x = E0;
Cu_Mermin.DECS.E0 = E0;
Cu_Mermin.Composition.Z = Cu_Mermin.Z;
Cu_Mermin.Composition.index = 1;

tic;
[data] = ElsepaRunner.RunElsepa(Cu_Mermin.Composition,E0);
toc
Cu_Mermin.DECS.x = data(1).x;
for i = 1:numel(E0)
    Cu_Mermin.Elastic.sigma_el(i) = data(i).sigma_el;
    Cu_Mermin.Elastic.l_el(i) = 1/data(i).sigma_el/Cu_Mermin.Density;
    Cu_Mermin.Elastic.l_tr(i) = 1/data(i).sigma_tr1/Cu_Mermin.Density;
    Cu_Mermin.DECS.y(:,i) = data(i).y;
end
%}

%% Inelastic properties
osc.model = 'Mermin';
osc.A = [0.01 0.02 0.06 0.11 0.38 0.15 0.07 0.05 0.08 0.1]';
osc.G = [0.69 0.8 2.36 6.59 10.82 4.89 40.77 19.15 14.12 161.76]; 
osc.Om =[3.46 4.2 7.54 28.34 19.39 10.48 66.56 47.52 36.13 126.04];
osc.alpha = 1; 
osc.beps = 1;
osc.Ef = Cu_Mermin.Ef; 
osc.qtran = 0.01:0.01:20;
osc.eloss = eps:.1:110;
osc.egap = eps;

Cu_Mermin.ELF = eps_sum_allwq(osc,'bulk');
Cu_Mermin.eloss = osc.eloss;
Cu_Mermin.q = osc.qtran;
Cu_Mermin.DIIMFP = zeros(N,2,numel(E0));
Cu_Mermin.l_in = zeros(numel(E0),1);
for i = 1:length(E0)
    if E0(i) > Cu_Mermin.Ef
        energy = E0(i) - Cu_Mermin.Ef;
        eloss = eps:(energy-eps)/(N-1):energy;
        if energy < 10
            osc.eloss = eps:0.01:energy;
        else
            osc.eloss = eps:0.2:energy;
        end
        Cu_Mermin.DIIMFP(:,1,i) = eloss;
        [iimfp, diimfp] = ndiimfp(osc,E0(i));
        Cu_Mermin.DIIMFP(:,2,i) = interp1(osc.eloss,diimfp,eloss);
        Cu_Mermin.l_in(i) = 1/trapz(osc.eloss/h2ev,iimfp)*a0;
    else
        Cu_Mermin.l_in(i) = Inf;
    end
end

%% Ionisation shells
Cu_Mermin.Shells = {'2S1/2';'2P1/2';'2P3/2';'3S1/2';'3P1/2';'3P3/2';'3D3/2';'3D5/2';'4S1/2';};
Cu_Mermin.EB = [1103;958;938;127;82;80;11;10;8;];

