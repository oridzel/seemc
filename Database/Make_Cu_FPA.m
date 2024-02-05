function Cu_FPA = Make_Cu_FPA

E0 = [1:100 110:10:200 220:20:300 350:50:500 600:100:2500 2750:250:5000 5500:500:30000];
N = 4000;

%% Basic
Cu_FPA.Mat = 'Cu_FPA';
Cu_FPA.M = 63.546;
Cu_FPA.Z = 29;
Cu_FPA.Density = 8.96; %g/cm^3
Cu_FPA.Density = Cu_FPA.Density*10^-24/Cu_FPA.M*6.022*10^23; %#/A^3
Cu_FPA.NvTPP = 11;
Cu_FPA.Ep = 35.87;
Cu_FPA.Ef = 8.7;
Cu_FPA.Wf = 4.65;
Cu_FPA.isMetal = true;

%% Elastic properties
% {
Cu_FPA.Elastic.x = zeros(numel(E0),1);
Cu_FPA.Elastic.l_el = zeros(numel(E0),1);
Cu_FPA.Elastic.l_tr = zeros(numel(E0),1);
Cu_FPA.Elastic.x = E0;
Cu_FPA.DECS.E0 = E0;
Cu_FPA.Composition.Z = Cu_FPA.Z;
Cu_FPA.Composition.index = 1;

tic;
[data] = ElsepaRunner.RunElsepa(Cu_FPA.Composition,E0);
toc
Cu_FPA.DECS.x = data(1).x;
for i = 1:numel(E0)
    Cu_FPA.Elastic.sigma_el(i) = data(i).sigma_el;
    Cu_FPA.Elastic.l_el(i) = 1/(data(i).sigma_el*a0^2)/Cu_FPA.Density;
    Cu_FPA.Elastic.l_tr(i) = 1/(data(i).sigma_tr1*a0^2)/Cu_FPA.Density;
    Cu_FPA.DECS.y(:,i) = data(i).y;
end
%}

%% Inelastic properties
osc.model = 'Mermin';
osc.A = [0.01 0.02 0.06 0.11 0.38 0.15 0.07 0.05 0.08 0.1]';
osc.G = [0.69 0.8 2.36 6.59 10.82 4.89 40.77 19.15 14.12 161.76]; 
osc.Om =[3.46 4.2 7.54 28.34 19.39 10.48 66.56 47.52 36.13 126.04];
osc.alpha = 1; 
osc.beps = 1;
osc.Ef = Cu_FPA.Ef; 
osc.qtran = 0.01:0.01:20;
osc.eloss = eps:.5:110;
osc.egap = eps;

opt_elf = load("C:/Users/onr5/OneDrive - NIST/dev/seemc/Data/opt/cu.diel");
% q = 0:0.1:10;
% eloss = 0:.5:200;
tic
load cu_fpa_elf.mat
% [Cu_FPA.ELF,~,~] = fpa_vector(q*a0,eloss/h2ev,opt_elf(:,1),opt_elf(:,4));
toc
Cu_FPA.eloss = omega;
Cu_FPA.q = q;
Cu_FPA.ELF = elf;
Cu_FPA.DIIMFP = zeros(N,2,numel(E0));
Cu_FPA.l_in = zeros(numel(E0),1);
for i = 1:length(E0)
    if E0(i) > Cu_FPA.Ef
        energy = E0(i) - Cu_FPA.Ef;
        eloss = eps:(energy-eps)/(N-1):energy;
        omega = eps:0.2:energy;
        Cu_FPA.DIIMFP(:,1,i) = eloss;
        [iimfp, diimfp_] = diimfp_interp_fpa(E0(i),omega,Cu_FPA.ELF,Cu_FPA.q,Cu_FPA.eloss,osc);
        Cu_FPA.DIIMFP(:,2,i) = interp1(omega,diimfp_,eloss);
        % Cu_FPA.DIIMFP(:,2,i) = diimfp_interp./trapz(eloss,diimfp_interp);
        Cu_FPA.l_in(i) = 1/trapz(omega/h2ev,iimfp)*a0;
    else
        Cu_FPA.l_in(i) = Inf;
    end
end

%% Ionisation shells
Cu_FPA.Shells = {'2S1/2';'2P1/2';'2P3/2';'3S1/2';'3P1/2';'3P3/2';'3D3/2';'3D5/2';'4S1/2';};
Cu_FPA.EB = [1103;958;938;127;82;80;11;10;8;];

