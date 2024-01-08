function H2O = Make_H2O

E0 = [1:100 150:50:500 600:100:2500 2750:250:5000 5500:500:30000];
N = 4000;

%% Basic
H2O = struct;
H2O.Mat = 'H2O';
H2O.Z = 10;
H2O.M = 18.0152;
H2O.Density = 0.999973; %g/cm^3
H2O.Density = H2O.Density*10^-24/H2O.M*6.022*10^23; %#/A^3
H2O.NvTPP = 8;
H2O.Eg = 7.9;
H2O.Ep = 19.2;
H2O.Evb = 10;
H2O.Affinity = 1.3;
H2O.Phonon.eps_zero = 8.1; % static dielectric constant
H2O.Phonon.eps_inf = 1.7876; % high-frequency dielectric constant (the square of the static refractive index)
H2O.Phonon.eloss = 0.1;

%% Elastic properties
% {
H2O.Elastic.x = zeros(numel(E0),1);
H2O.Elastic.l_el = zeros(numel(E0),1);
H2O.Elastic.l_tr = zeros(numel(E0),1);
H2O.Elastic.x = E0;
H2O.DECS.E0 = E0;
H2O.Composition.Z = [1 8];
H2O.Composition.index = [2 1];

tic;
[data] = ElsepaRunner.RunElsepa(H2O.Composition,E0);
toc
H2O.DECS.x = data(1).x;
for i = 1:numel(E0)
    H2O.Elastic.l_el(i) = 1/data(i).sigma_el/H2O.Density;
    H2O.Elastic.l_tr(i) = 1/data(i).sigma_tr1/H2O.Density;
    H2O.DECS.y(:,i) = data(i).y/trapz(data(i).x,data(i).y);
end
%}

%% Inelastic properties
osc.model = 'Mermin';
osc.A = [0.2939 0.1608 0.0060 0.0014 0.0899 0.1499 0.0634 0.1459 0.0888];
osc.G = [5.4811 5.5443 0.3747 0.2509 0.2501 6.2774 1.8659 4.8673 7.7579];
osc.Om = [13.1653 16.0342 5.4358 45.1881 35.2055 18.4556 28.8435 9.9877 22.1646];
osc.alpha = 1; 
osc.beps = 1;
osc.Ef = H2O.Evb; 
osc.qtran = 0.01:0.01:20;
osc.eloss = eps:.1:110;
osc.egap = H2O.Eg;

H2O.ELF = eps_sum_allwq(osc,'bulk');
H2O.eloss = osc.eloss;
H2O.q = osc.qtran;
H2O.DIIMFP = zeros(N,2,numel(E0));
H2O.l_in = zeros(numel(E0),1);
for i = 1:length(E0)
    if E0(i) > 2*H2O.Eg + H2O.Evb
        energy = E0(i) - H2O.Eg - H2O.Evb;
        osc.eloss = H2O.Eg:(energy-H2O.Eg)/(N-1):energy;
        H2O.DIIMFP(:,1,i) = osc.eloss;
        [iimfp, H2O.DIIMFP(:,2,i)] = ndiimfp(osc,E0(i));
        H2O.l_in(i) = 1/trapz(osc.eloss/h2ev,iimfp)*a0;
    else
        H2O.l_in(i) = Inf;
    end
end

end