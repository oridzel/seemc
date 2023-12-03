clear;

N = 1000;
E0 = [50 100 200 500 700 1000];
e = cell(length(E0),1);
matname = 'Si';

tic
for i = 1:length(E0)
    disp(E0(i))
    e{i} = simulateSEE(N,E0(i),matname);
end
toc

%% Ready results
% load Si_DFT_e_array.mat;

%% Histograms
% spec_se = zeros(1);
% spec_pe = zeros(1);
sey = zeros(size(E0));
bse = zeros(size(E0));
for i = 1:length(E0)
    bse(i) = 0;
    sey(i) = 0;
    for j = 1:N
        for k = 1:length(e{i}{j})
            if ~e{i}{j}(k).Inside && ~e{i}{j}(k).Dead
                if e{i}{j}(k).isSecondary                    
                    sey(i) = sey(i) + 1;
                else
                    bse(i) = bse(i) + 1;
                end
            end
        end
    end
end

%% Yields
% {
figure
hold on
box on
plot(E0,(sey+bse)/N,DisplayName='TEY',LineWidth=2)
plot(E0,bse/N,DisplayName='BSE',LineWidth=2)
plot(E0,sey/N,DisplayName='SEY',LineWidth=2)
xlabel('Energy (eV)')
ylabel('Yield')
title(matname)
fontsize(16,"points")
legend
%}

%% Plots
%{
figure
hold on
histogram(spec_pe,500)
histogram(spec_se,500)
xlabel('Electron energy (eV)')
%}

