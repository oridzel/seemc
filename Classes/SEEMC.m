classdef SEEMC < handle
    properties
        numTrajectories = 1000
        matName = 'Au'
        isMetal = true
        trackTrajectories = false
        energyArray
        sample
        statistics
        onlyEscaped = true
        bse
        sey
        energyHistogramPE
        energyHistogramSE
        coincidenceHistogram
    end
    methods
        function obj = SEEMC(inputpar)
            obj.matName = inputpar.matName;
            obj.isMetal = inputpar.isMetal;
            obj.numTrajectories = inputpar.numTrajectories;
            obj.trackTrajectories = inputpar.trackTrajectories;
            obj.energyArray = inputpar.energy;
            obj.onlyEscaped = inputpar.onlyEscaped;
            obj.sample = Sample(obj.matName,obj.isMetal);
        end
        function simulate(obj)

            obj.statistics = cell(1);
            n_traj = obj.numTrajectories;
            energy_array = obj.energyArray;
            smpl = obj.sample;
            track = obj.trackTrajectories;
            stat = cell(size(energy_array));
            saveEscaped = obj.onlyEscaped;

            for e = 1:numel(energy_array)
                energy = energy_array(e);
                disp(energy)
                electronData = cell(n_traj,1);            
                parfor i = 1:n_traj
                    e_count = 0;
                    res = Electron.empty;
                    res(end+1) = Electron(energy,smpl,track);
                    while e_count < length(res)
                        e_count = e_count + 1;
                        while res(e_count).Inside && ~res(e_count).Dead
                            res(e_count).travel;
                            res(e_count).escape;
                            if res(e_count).Inside && ~res(e_count).Dead
                                res(e_count).getScatteringType;
                                if res(e_count).scatter
                                    e_se = res(e_count).EnergyLoss + res(e_count).EnergySE;
                                    if e_se > res(e_count).InnerPotential
                                        uvw = [ sin(acos(2*rand-1))*cos(2*rand*pi),...
                                                sin(acos(2*rand-1))*sin(2*rand*pi),...
                                                cos(acos(2*rand-1)) ];
                                        res(end + 1) = Electron(e_se,smpl,track,res(e_count).xyz,uvw,res(e_count).nSecondaries+1,true,e_count);
                                        res(e_count).nSecondaries = res(e_count).nSecondaries + 1;
                                    end
                                end
                            end            
                        end
                        y = struct();
                        if ~res(e_count).Dead && ~res(e_count).Inside && saveEscaped
                            y.Energy = res(e_count).Energy;
                            y.InnerPotential = res(e_count).InnerPotential;
                            y.EnergyLoss = res(e_count).EnergyLoss;
                            y.EnergySE = res(e_count).EnergySE;
                            y.Angles = res(e_count).Angles;
                            y.isSecondary = res(e_count).isSecondary;
                            y.Generation = res(e_count).Generation;
                            y.ParentIndex = res(e_count).ParentIndex;
                            y.Dead = res(e_count).Dead;
                            y.Inside = res(e_count).Inside;
                            if track
                                y.coordinates = res(e_count).coordinates;
                            end
                            electronData{i}(end+1) = y;
                        elseif ~saveEscaped
                            y.Energy = res(e_count).Energy;
                            y.InnerPotential = res(e_count).InnerPotential;
                            y.EnergyLoss = res(e_count).EnergyLoss;
                            y.EnergySE = res(e_count).EnergySE;
                            y.Angles = res(e_count).Angles;
                            y.isSecondary = res(e_count).isSecondary;
                            y.Generation = res(e_count).Generation;
                            y.ParentIndex = res(e_count).ParentIndex;
                            y.Dead = res(e_count).Dead;
                            y.Inside = res(e_count).Inside;
                            if track
                                y.coordinates = res(e_count).coordinates;
                            end
                            electronData{i}(end+1) = y;
                        end
                    end
                end
                stat{e} = electronData;
            end
            obj.statistics = stat;
        end

        function calculateYields(obj)
            obj.sey = zeros(size(obj.energyArray));
            obj.bse = zeros(size(obj.energyArray));
            for i = 1:numel(obj.energyArray)
                for j = 1:obj.numTrajectories
                    for k = 1:length(obj.statistics{i}{j})
                        if ~obj.statistics{i}{j}(k).Inside && ~obj.statistics{i}{j}(k).Dead
                            if obj.statistics{i}{j}(k).isSecondary
                                obj.sey(i) = obj.sey(i) + 1;
                            else
                                obj.bse(i) = obj.bse(i) + 1;
                            end
                        end
                    end
                end
            end
            obj.sey = obj.sey/obj.numTrajectories;
            obj.bse = obj.bse/obj.numTrajectories;
        end

        function calculateEnergyHistograms(obj)
            obj.energyHistogramPE = cell(numel(obj.energyArray),1);
            obj.energyHistogramSE = cell(numel(obj.energyArray),1);
            for i = 1:numel(obj.energyArray)
                for j = 1:obj.numTrajectories
                    for k = 1:length(obj.statistics{i}{j})
                        if ~obj.statistics{i}{j}(k).Inside && ~obj.statistics{i}{j}(k).Dead
                            if obj.statistics{i}{j}(k).isSecondary
                                obj.energyHistogramSE{i}(end+1) = obj.statistics{i}{j}(k).Energy;
                            else
                                obj.energyHistogramPE{i}(end+1) = obj.statistics{i}{j}(k).Energy;
                            end
                        end
                    end
                end
            end
        end

        function plotEnergyDistribution(obj,ind)
            figure
            hold on
            box on
            histogram(obj.energyHistogramPE{ind},200,DisplayName='Primaries')
            histogram(obj.energyHistogramSE{ind},200,DisplayName='Secondaries')
            xlabel('Electron energy (eV)')
            ylabel('Counts')
            fontsize(20,"points")
            title(obj.matName)
            legend
        end

        function calculateCoincidenceHistogram(obj)
            if obj.onlyEscaped
                error('For coincidence histogram all electron trajectories must be stored.')
            end
            obj.coincidenceHistogram = cell(numel(obj.energyArray),1);
            for i = 1:numel(obj.energyArray)
                for j = 1:obj.numTrajectories
                    for k = 1:length(obj.statistics{i}{j})
                        if ~obj.statistics{i}{j}(k).Dead && obj.statistics{i}{j}(k).isSecondary
                            if ~obj.statistics{i}{j}(obj.statistics{i}{j}(k).ParentIndex).Dead && ...
                                    obj.statistics{i}{j}(k).Energy + obj.statistics{i}{j}(k).InnerPotential == ...
                                    obj.statistics{i}{j}(obj.statistics{i}{j}(k).ParentIndex).EnergyLoss + obj.statistics{i}{j}(obj.statistics{i}{j}(k).ParentIndex).EnergySE
                                obj.coincidenceHistogram{i}(end+1,:) = [obj.statistics{i}{j}(obj.statistics{i}{j}(k).ParentIndex).Energy obj.statistics{i}{j}(k).Energy];
                            end
                        end
                    end
                end
            end
        end

        function plotCoincidenceHistogram(obj,ind)
            figure
            hold on
            box on
            histogram2(obj.coincidenceHistogram{ind}(:,1),obj.coincidenceHistogram{ind}(:,2),50,'FaceColor','flat')
            xlabel('Electron energy (eV)')
            ylabel('Electron energy (eV)')
            xlim([0,obj.energyArray(ind)])
            colormap('turbo')
            colorbar
            fontsize(20,"points")
            title(obj.matName)
        end
    end
end