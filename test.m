clear;

N = 50;
% E0 = [20:10:100 200 500 700 800 1000];
E0 = 500;
e = cell(length(E0),1);
matname = 'Si';

tStart = tic;
for i = 1:length(E0)
    disp(E0(i))
    tic
    e{i} = simulateSEE(N,E0(i),matname,true);
    toc
end
tEnd = toc(tStart);
disp(['Total elapsed time over energies is ' num2str(tEnd) ' seconds.'])

%% Ready results
% load Si_DFT_e_array.mat;

%% Histograms
spec_se = cell(size(E0));
spec_pe = cell(size(E0));
coord_pe = cell(1,1);
coord_se = cell(1,1);
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
                    coord_se{end+1}(:,:) = e{i}{j}(k).coordinates;
                    sey(i) = sey(i) + 1;
                else
                    spec_pe{i}(end+1) = e{i}{j}(k).Energy;
                    coord_pe{end+1}(:,:) = e{i}{j}(k).coordinates;
                    bse(i) = bse(i) + 1;
                end
            end
        end
    end
    disp([num2str(E0(i)) '  ' num2str((sey(i)+bse(i))/N)])
end


%% Trajectories
figure
hold on
for i = 2:length(coord_pe)
    plot3(coord_pe{i}(:,1),coord_pe{i}(:,2),coord_pe{i}(:,3),Color='Blue')
end

for i = 2:length(coord_se)
    plot3(real(coord_se{i}(:,1)),real(coord_se{i}(:,2)),real(coord_se{i}(:,3)),Color='Red')
end

scatter(0,0,150,Color='Black')

%% Yields
%{
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

%% Spectra
%{
figure
hold on
histogram(spec_se{5},"NumBins",100)
histogram(spec_pe{5},"NumBins",100)
xlabel('Electron energy (eV)')
%}

