clear

inputpar.matName = { 'PS_FPA' };
% inputpar.thickness = 50;
inputpar.numTrajectories = 10000;
inputpar.energy = [10:10:100 200:50:500 600:100:1000];
% inputpar.energy = [2000 3000];

s = SEEMC(inputpar);
s.theta_0 = 0;
s.onlyEscaped = true;
s.cbRef = false;
s.trackTrajectories = false;

%% Run simulation
tic
s.simulate; 
toc

%% Trajectories
%{
s.getTrajectories;
s.plotTrajectories(1);
%}

%% Spectra
%{
s.calculateEnergyHistograms(1);
s.plotEnergyDistribution;
%}

%% Yields
% {
s.calculateYields;
s.plotYields;  
%}

%%
%{
sey = zeros(size(s.energyArray));
tey = zeros(size(s.energyArray));
bse = zeros(size(s.energyArray));
for i = 1:numel(s.energyArray)
    for j = 1:s.numTrajectories
        for k = 1:length(s.statistics{i}{j})
            if ~s.statistics{i}{j}(k).Inside && ~s.statistics{i}{j}(k).Dead
                tey(i) = tey(i) + 1;
                if s.statistics{i}{j}(k).isSecondary || s.statistics{i}{j}(k).Energy <= 50
                    sey(i) = sey(i) + 1;
                else
                    bse(i) = bse(i) + 1;
                end
            end
        end
    end
end
tey = tey/s.numTrajectories;
sey = sey/s.numTrajectories;
bse = bse/s.numTrajectories;
%}