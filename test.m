clear;

N = 10000;
E0 = [25 50 75 100 200 500 700 1000];
e = cell(length(E0),1);
<<<<<<< HEAD
matname = 'Si_DFT_b1l0';
=======
matname = 'Si_DL';
>>>>>>> 12bbd384c62a2acbfecc98add9e9c72919291e97

tStart = tic;
for i = 1:length(E0)
    disp(E0(i))
    tic
    e{i} = simulateSEE(N,E0(i),matname);
    toc
end
tEnd = toc(tStart);
disp(['Total elapsed time over energies is ' num2str(tEnd) ' seconds.'])

%% Ready results
% load Si_DFT_e_array.mat;

%% Histograms
spec_se = cell(size(E0));
spec_pe = cell(size(E0));
sey = zeros(size(E0));
bse = zeros(size(E0));
for i = 1:length(E0)
    bse(i) = 0;
    sey(i) = 0;
    for j = 1:N
        for k = 1:length(e{i}{j})
            if ~e{i}{j}(k).Inside && ~e{i}{j}(k).Dead
                if e{i}{j}(k).isSecondary
                    spec_se{i}(end+1) = e{i}{j}(k).Energy;
                    sey(i) = sey(i) + 1;
                else
                    spec_pe{i}(end+1) = e{i}{j}(k).Energy;  
                    bse(i) = bse(i) + 1;
                end
            end
        end
    end
    disp([num2str(E0(i)) '  ' num2str((sey(i)+bse(i))/N)])
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
histogram(spec_se{5},"NumBins",100)
histogram(spec_pe{5},"NumBins",100)
xlabel('Electron energy (eV)')
%}

