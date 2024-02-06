function Au_Mermin = Make_Au_Mermin

E0 = [1:100 110:10:200 220:20:300 350:50:500 600:100:2500 2750:250:5000 5500:500:30000];
N = 4000;

%% Basic
Au_Mermin.Mat = 'Au_Mermin';
Au_Mermin.M = 196.967;
Au_Mermin.Z = 79;
Au_Mermin.Density = 19.32; %g/cm^3
Au_Mermin.Density = Au_Mermin.Density*10^-24/Au_Mermin.M*6.022*10^23; %#/A^3
Au_Mermin.NvTPP = 11;
Au_Mermin.Ep = 29.92;
Au_Mermin.Ef = 9;
Au_Mermin.Wf = 5.2;
Au_Mermin.isMetal = true;

%% Elastic properties
% {
Au_Mermin.Elastic.x = zeros(numel(E0),1);
Au_Mermin.Elastic.l_el = zeros(numel(E0),1);
Au_Mermin.Elastic.l_tr = zeros(numel(E0),1);
Au_Mermin.Elastic.x = E0;
Au_Mermin.DECS.E0 = E0;
Au_Mermin.Composition.Z = Au_Mermin.Z;
Au_Mermin.Composition.index = 1;

tic;
[data] = ElsepaRunner.RunElsepa(Au_Mermin.Composition,E0);
toc
Au_Mermin.DECS.x = data(1).x;
for i = 1:numel(E0)
    Au_Mermin.Elastic.sigma_el(i) = data(i).sigma_el;
    Au_Mermin.Elastic.l_el(i) = 1/data(i).sigma_el/Au_Mermin.Density;
    Au_Mermin.Elastic.l_tr(i) = 1/data(i).sigma_tr1/Au_Mermin.Density;
    Au_Mermin.DECS.y(:,i) = data(i).y;
end
%}

%% Inelastic properties
osc.model = 'Mermin';
osc.A = [0.0312 0.0429 0.1253 0.0469 0.2512 0.0792 0.1284 0.0517 0.0459 0.0515 0.1158 0.0300];
osc.G = [0.9698 3.2772 6.4314 3.5576 8.4554 6.9578 12.2837 12.6974 14.2658 4.4995 24.8269 25.0628];
osc.Om = [3.2380 6.2638 11.6974 17.5575 26.2625 33.3481 38.4167 45.6709 53.1109 15.2746 64.7666 84.5982];
osc.alpha = 1; 
osc.beps = 1;
osc.Ef = Au_Mermin.Ef; 
osc.qtran = 0.01:0.01:20;
osc.eloss = eps:.1:110;
osc.egap = eps;

Au_Mermin.ELF = eps_sum_allwq(osc,'bulk');
Au_Mermin.eloss = osc.eloss;
Au_Mermin.q = osc.qtran;
Au_Mermin.DIIMFP = zeros(N,2,numel(E0));
Au_Mermin.l_in = zeros(numel(E0),1);
for i = 1:length(E0)
    if E0(i) > Au_Mermin.Ef
        energy = E0(i) - Au_Mermin.Ef;
        eloss = eps:(energy-eps)/(N-1):energy;
        if energy < 10
            osc.eloss = eps:0.01:energy;
        else
            osc.eloss = eps:0.2:energy;
        end
        Au_Mermin.DIIMFP(:,1,i) = eloss;
        [iimfp, diimfp] = ndiimfp(osc,E0(i));
        Au_Mermin.DIIMFP(:,2,i) = interp1(osc.eloss,diimfp,eloss);
        Au_Mermin.l_in(i) = 1/trapz(osc.eloss/h2ev,iimfp)*a0;
    else
        Au_Mermin.l_in(i) = Inf;
    end
end

%% Ionisation shells
Au_Mermin.Shells = {'4S1/2';'4P1/2';'4P3/2';'4D3/2';'4D5/2';'4F5/2';'4F7/2';'5S1/2';'5P1/2';'5P3/2';'5D3/2';'5D5/2';'6S1/2';};
Au_Mermin.EB = [763;646;549;356;339;91;87;113;76;60;10;10;9;];

