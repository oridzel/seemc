   function Au_FPA = Make_Au_FPA

E0 = [1:100 150:50:500 600:100:2500 2750:250:5000]; % 5500:500:30000];
N = 4000;

%% Basic
Au_FPA.Mat = 'Au_FPA';
Au_FPA.M = 196.967;      
Au_FPA.Density = 19.32; %g/cm^3
Au_FPA.Density = Au_FPA.Density*10^-24/Au_FPA.M*6.022*10^23; %#/A^3
Au_FPA.NvTPP = 11;
Au_FPA.Ep = 29.92;
Au_FPA.Ef = 9;
Au_FPA.Wf = 5.2;

%% Elastic properties
%{
Au_FPA.Elastic.x = zeros(numel(E0),1);
Au_FPA.Elastic.l_el = zeros(numel(E0),1);
Au_FPA.Elastic.l_tr = zeros(numel(E0),1);
Au_FPA.Elastic.x = E0;
Au_FPA.DECS.E0 = E0;
Au_FPA.Composition.Z = Au_FPA.Z;
Au_FPA.Composition.index = 1;

tic;
[data] = ElsepaRunner.RunElsepa(Au_FPA.Composition,E0);
toc
Au_FPA.DECS.x = data(1).x;
for i = 1:numel(E0)
    Au_FPA.Elastic.l_el(i) = 1/data(i).sigma_el/Au_FPA.Density;
    Au_FPA.Elastic.l_tr(i) = 1/data(i).sigma_tr1/Au_FPA.Density;
    Au_FPA.DECS.y(:,i) = data(i).y/trapz(data(i).x,data(i).y);
end
%}

%% Inelastic properties

opt_elf = load("../Data/opt/au.diel");
q = 0:1:10;
eloss = 0:.5:110;

tic
[Au_FPA.ELF,~,~] = fpa_vector(q*a0,eloss/h2ev,opt_elf(:,1),opt_elf(:,4));
toc
Au_FPA.eloss = eloss;
Au_FPA.q = q;
Au_FPA.DIIMFP = zeros(N,2,numel(E0));
Au_FPA.l_in = zeros(numel(E0),1);
for i = 1:length(E0)
    if E0(i) > Au_FPA.Ef
        energy = E0(i) - Au_FPA.Ef;
        eloss = eps:(energy-eps)/(N-1):energy;
        omega = eps:0.5:energy;
        Au_FPA.DIIMFP(:,1,i) = eloss;
        [iimfp, diimfp_] = diimfp(E0(i),omega,opt_elf(:,1),opt_elf(:,4));
        Au_FPA.DIIMFP(:,2,i) = interp1(omega,diimfp_,eloss);
        Au_FPA.l_in(i) = 1/trapz(omega/h2ev,iimfp)*a0;
    else
        Au_FPA.l_in(i) = Inf;
    end
 end

%% Ionisation shells
Au_FPA.Shells = {'4S1/2';'4P1/2';'4P3/2';'4D3/2';'4D5/2';'4F5/2';'4F7/2';'5S1/2';'5P1/2';'5P3/2';'5D3/2';'5D5/2';'6S1/2';};
Au_FPA.EB = [763;646;549;356;339;91;87;113;76;60;10;10;9;];

