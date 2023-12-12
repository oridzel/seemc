clear;

run init.m

N = 100;
% E0 = [20:10:100 200 500 700 800 1000];
E0 = 100;
e = cell(length(E0),1);
matname = 'Au';
isMetal = true;
trackTrajectories = false;
coincidences = false;

tStart = tic;
for i = 1:length(E0)
    disp(E0(i))
    tic
    e{i} = simulateSEE(N,E0(i),matname,isMetal,trackTrajectories);
    toc
end
tEnd = toc(tStart);
disp(['Total elapsed time over energies is ' num2str(tEnd) ' seconds.'])

%% Histograms
spec_se = cell(size(E0));
spec_pe = cell(size(E0));
if numel(E0) == 1
    if trackTrajectories
        coord_pe = cell(1,1);
        coord_pe_dead = cell(1,1);
        coord_se = cell(1,1);
        coord_se_dead = cell(1,1);
    end
    if coincidences
        coin_arr = zeros(1,2);
    end
end
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
                    if numel(E0) == 1
                        if trackTrajectories
                            coord_se{end+1}(:,:) = e{i}{j}(k).coordinates;
                        end
                        if coincidences
                            if ~e{i}{j}(e{i}{j}(k).ParentIndex).Dead && ~e{i}{j}(e{i}{j}(k).ParentIndex).Inside && e{i}{j}(k).Energy + e{i}{j}(k).InnerPotential == e{i}{j}(e{i}{j}(k).ParentIndex).EnergyLoss + e{i}{j}(e{i}{j}(k).ParentIndex).EnergySE
                                coin_arr(end+1,:) = [e{i}{j}(e{i}{j}(k).ParentIndex).Energy e{i}{j}(k).Energy];
                            end
                        end
                    end
                    sey(i) = sey(i) + 1;
                else
                    spec_pe{i}(end+1) = e{i}{j}(k).Energy;
                    if numel(E0) == 1 && trackTrajectories
                        coord_pe{end+1}(:,:) = e{i}{j}(k).coordinates;
                    end
                    bse(i) = bse(i) + 1;
                end
            elseif e{i}{j}(k).Inside && e{i}{j}(k).Dead && numel(E0) == 1 && trackTrajectories
                if e{i}{j}(k).isSecondary
                    coord_se_dead{end+1}(:,:) = e{i}{j}(k).coordinates;
                else
                    coord_pe_dead{end+1}(:,:) = e{i}{j}(k).coordinates;
                end
            end
        end
    end
    disp([num2str(E0(i)) '  ' num2str((sey(i)+bse(i))/N)])
end
save([matname '_yields.mat'],"E0","bse","sey")

%% Trajectories
%{
figure
title([matname ' E_0 = ' num2str(E0) ' eV'])
hold on
box on
colormap(turbo)

% 3D
%{
% primaries
for i = 2:length(coord_pe)    
    patch([coord_pe{i}(:,1); nan],[coord_pe{i}(:,2); nan],[coord_pe{i}(:,3); nan],[coord_pe{i}(:,4); nan],'FaceColor','none','EdgeColor','interp')
end
% dead primaries
%{
for i = 2:length(coord_pe_dead)    
    patch([coord_pe_dead{i}(:,1); nan],[coord_pe_dead{i}(:,2); nan],[coord_pe_dead{i}(:,3); nan],[coord_pe_dead{i}(:,4); nan],'FaceColor','none','EdgeColor','interp')
end
%}
% secondaries
for i = 2:length(coord_se)
    patch([coord_se{i}(:,1); nan],[coord_se{i}(:,2); nan],[coord_se{i}(:,3); nan],[coord_se{i}(:,4); nan],'FaceColor','none','EdgeColor','interp')
end
% dead secondaries
%{
for i = 2:length(coord_se_dead)
    patch([coord_se_dead{i}(:,1); nan],[coord_se_dead{i}(:,2); nan],[coord_se_dead{i}(:,3); nan],[coord_se_dead{i}(:,4); nan],'FaceColor','none','EdgeColor','interp')
end
%}

colorbar
view(3)

scatter(0,0,50,'filled',Color='Black')
set(gca,'Zdir','reverse')
%}

% 2D
%{
% primaries
for i = 2:length(coord_pe)    
    patch([coord_pe{i}(:,1); nan],[coord_pe{i}(:,3); nan],[coord_pe{i}(:,4); nan],'FaceColor','none','EdgeColor','interp')
end
% dead primaries
% {
for i = 2:length(coord_pe_dead)    
    patch([coord_pe_dead{i}(:,1); nan],[coord_pe_dead{i}(:,3); nan],[coord_pe_dead{i}(:,4); nan],'FaceColor','none','EdgeColor','interp')
end
%}
% secondaries
for i = 2:length(coord_se)
    patch([coord_se{i}(:,1); nan],[coord_se{i}(:,3); nan],[coord_se{i}(:,4); nan],'FaceColor','none','EdgeColor','interp')
end
% dead secondaries
%{
for i = 2:length(coord_se_dead)
    patch([coord_se_dead{i}(:,1); nan],[coord_se_dead{i}(:,3); nan],[coord_se_dead{i}(:,4); nan],'FaceColor','none','EdgeColor','interp')
end
%}

scatter(0,0,50,'filled',Color='Black')
set(gca,'Ydir','reverse')
c = colorbar;
ylabel(c,'Energy (eV)');
xlabel('X (A)');
ylabel('Z (A)')

xl = xlim;
yl = ylim;
x = [xl(1) xl(2) xl(2) xl(1)];
y = [0 0 yl(2) yl(2)];
fill(x,y,'b','FaceAlpha',0.2)

fontsize(20,'points')
%}
%}

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
histogram(spec_pe{1},"NumBins",200,DisplayName="Primaries")
histogram(spec_se{1},"NumBins",200,DisplayName="Secondaries")
xlabel('Electron energy (eV)')
legend
%}

