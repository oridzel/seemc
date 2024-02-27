function PS_MLL = Make_PS_MLL

E0 = [1:100 110:10:200 220:20:300 350:50:500 600:100:2500 2750:250:5000 5500:500:30000];
N = 4000;

%% Basic
PS_MLL = struct;
PS_MLL.Mat = 'PS_MLL';
PS_MLL.M = 104.1440;
PS_MLL.Z = 56;
PS_MLL.Density = 1.05; %g/cm^3
PS_MLL.Density = PS_MLL.Density*10^-24/PS_MLL.M*6.022*10^23; %#/A^3
PS_MLL.NvTPP = 40;
PS_MLL.Eg = 4.5;
PS_MLL.Ep = 18.3;
PS_MLL.Evb = 17;
PS_MLL.Affinity = 3.5;
PS_MLL.Phonon.eps_zero = 2.5; % static dielectric constant
PS_MLL.Phonon.eps_inf = 1.01; % high-frequency dielectric constant (the square of the static refractive index)
PS_MLL.Phonon.eloss = 0.3;
PS_MLL.isMetal = false;

%% Elastic properties
PS_MLL.Elastic.x = zeros(numel(E0),1);
PS_MLL.Elastic.l_el = zeros(numel(E0),1);
PS_MLL.Elastic.l_tr = zeros(numel(E0),1);
PS_MLL.Elastic.x = E0;
PS_MLL.DECS.E0 = E0;
PS_MLL.Composition.Z = [6 1];
PS_MLL.Composition.index = [8 8];

tic;
[data] = ElsepaRunner.RunElsepa(PS_MLL.Composition,E0,0);
toc
PS_MLL.DECS.x = data(1).x;
for i = 1:numel(E0)
    PS_MLL.Elastic.sigma_el(i) = data(i).sigma_el;
    PS_MLL.Elastic.l_el(i) = 1/(data(i).sigma_el*a0^2*PS_MLL.Density);
    PS_MLL.Elastic.l_tr(i) = 1/(data(i).sigma_tr1*a0^2*PS_MLL.Density);
    PS_MLL.DECS.y(:,i) = data(i).y;
end

%% Inelastic properties
osc.model = 'MerminLL';
osc.A = [0.045 0.013 0.011 0.053 0.168 0.093 0.061 0.018 0.066 0.045 0.021 0.01];
osc.G = [0.791 1.931 1.55 3.763 5.415 4.305 4.611 7.793 5.175 5.926 8.526 13.125];
osc.Om = [6.749 11.774 9.924 14.751 21.154 17.898 23.453 37.127 25.588 28.927 32.841 46.14];
osc.u = 1.6415029200190399;
osc.alpha = 1;
osc.beps = 1;
osc.Ef = PS_MLL.Evb; 
osc.qtran = 0.01:0.01:20;
osc.eloss = eps:.1:110;
osc.egap = PS_MLL.Eg;

PS_MLL.ELF = eps_sum_allwq(osc,'bulk');
PS_MLL.eloss = osc.eloss;
PS_MLL.q = osc.qtran;
PS_MLL.DIIMFP = zeros(N,2,numel(E0));
PS_MLL.l_in = zeros(numel(E0),1);
for i = 1:length(E0)
    if E0(i) > 2*PS_MLL.Eg + PS_MLL.Evb
        energy = E0(i) - PS_MLL.Eg - PS_MLL.Evb;
        eloss = PS_MLL.Eg:(energy-PS_MLL.Eg)/(N-1):energy;
        if energy < 10
            osc.eloss = PS_MLL.Eg:0.01:energy;
        else
            osc.eloss = PS_MLL.Eg:0.2:energy;
        end
        PS_MLL.DIIMFP(:,1,i) = eloss;
        [iimfp, diimfp] = ndiimfp(osc,E0(i));
        PS_MLL.DIIMFP(:,2,i) = interp1(osc.eloss,diimfp,eloss);
        PS_MLL.l_in(i) = 1/trapz(osc.eloss/h2ev,iimfp)*a0;
    else
        PS_MLL.l_in(i) = Inf;
    end
end

