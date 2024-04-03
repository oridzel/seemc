classdef SEEMC < handle
    properties
        numTrajectories = 1000
        matName = 'Au'
        theta_0 = 0
        cbRef = false
        trackTrajectories = false
        energyArray
        layers(1,:) Layer
        statistics
        onlyEscaped = true
        bse
        sey
        energyHistogramSE
        energyHistogramPE
        depthHistogramSE
        coincidenceHistogram
        coordinateArray
    end
    methods
        function obj = SEEMC(inputpar)
            obj.matName = inputpar.matName;
            obj.numTrajectories = inputpar.numTrajectories;
            obj.energyArray = inputpar.energy;

            if length(obj.matName) == 1
                obj.layers = Layer(Sample(obj.matName{1}));
            elseif length(obj.matName) > 1
                if ~isfield(inputpar,'thickness') || length(inputpar.thickness) ~= length(obj.matName)-1
                    error('Set thickness for each layer');
                end
                for i = 1:length(obj.matName)-1
                    obj.layers(i) = Layer(Sample(obj.matName{i}),inputpar.thickness(i));
                end
                obj.layers(length(obj.matName)) = Layer(Sample(obj.matName{i}));
            end
        end
        function simulate(obj)

            obj.statistics = cell(1);
            n_traj = obj.numTrajectories;
            energy_array = obj.energyArray;
            lyrs = obj.layers;
            track = obj.trackTrajectories;
            cb = obj.cbRef;
            stat = cell(size(energy_array));
            saveEscaped = obj.onlyEscaped;
            theta = obj.theta_0;

            for e = 1:numel(energy_array)
                % h = waitbar(0, 'Starting...',Name=['Simulation for E = ',num2str(energy_array(e)),' eV']);
                energy = energy_array(e); 
                disp(energy)
                electronData = cell(n_traj,1);       
                parfor i = 1:n_traj
                    % waitbar(i/n_traj, h, ['Progress: ', num2str(floor(i/n_traj*100)),'%']);
                    e_count = 0;
                    res = Electron.empty;
                    res(end+1) = Electron(energy,lyrs,'cbRef',cb,'trackCoordinates',track,'theta',theta);
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
                                        % uvw = [ sin(acos(2*rand-1))*cos(2*rand*pi),...
                                        %         sin(acos(2*rand-1))*sin(2*rand*pi),...
                                        %         cos(acos(2*rand-1)) ];
                                        uvw = updateDirection(res(e_count).uvw,[asin(cos(res(e_count).Deflection(1))),res(e_count).Deflection(2) + pi],1);
                                        res(end + 1) = Electron(e_se,lyrs,'cbRef',cb,'trackCoordinates',track,'xyz',res(e_count).xyz,...
                                                                'uvw',uvw,'generation',res(e_count).nSecondaries+1,'isSecondary',true,'index',e_count);
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
                            y.isSecondary = res(e_count).isSecondary;
                            y.Generation = res(e_count).Generation;
                            y.ParentIndex = res(e_count).ParentIndex;
                            y.Dead = res(e_count).Dead;
                            y.Inside = res(e_count).Inside;
                            y.InitialDepth = res(e_count).InitialDepth;
                            if track
                                y.Coordinates = res(e_count).coordinates;
                            end
                            electronData{i}(end+1) = y;
                        elseif ~saveEscaped
                            y.Energy = res(e_count).Energy;
                            y.InnerPotential = res(e_count).InnerPotential;
                            y.EnergyLoss = res(e_count).EnergyLoss;
                            y.EnergySE = res(e_count).EnergySE;
                            y.isSecondary = res(e_count).isSecondary;
                            y.Generation = res(e_count).Generation;
                            y.ParentIndex = res(e_count).ParentIndex;
                            y.Dead = res(e_count).Dead;
                            y.Inside = res(e_count).Inside;
                            y.InitialDepth = res(e_count).InitialDepth;
                            if track
                                y.Coordinates = res(e_count).coordinates;
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
                            if obj.statistics{i}{j}(k).isSecondary % || obj.statistics{i}{j}(k).Energy <= 50
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

        function plotYields(obj)
            figure
            hold on
            box on
            plot(obj.energyArray,obj.sey+obj.bse,DisplayName='TEY',LineWidth=2)
            plot(obj.energyArray,obj.sey,DisplayName='SEY',LineWidth=2)
            plot(obj.energyArray,obj.bse,DisplayName='BSE',LineWidth=2)
            xlabel('Energy (eV)')
            ylabel('Yield')
            fontsize(20,"points")
            title(strrep(obj.matName,'_',' '))
            legend
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

        function plotEnergyDistribution(obj,ind,nbins,varargin)
            if nargin < 3
                nbins = 200;
            end
            figure
            hold on
            box on
            histogram(obj.energyHistogramPE{ind},nbins,DisplayName='Primaries')
            histogram(obj.energyHistogramSE{ind},nbins,DisplayName='Secondaries')
            xlabel('Electron energy (eV)')
            ylabel('Counts')
            fontsize(20,"points")
            title(strrep(obj.matName,'_',' '))
            legend
        end

        function calculateDepthHistograms(obj)
            obj.depthHistogramSE = cell(numel(obj.energyArray),1);
            for i = 1:numel(obj.energyArray)
                for j = 1:obj.numTrajectories
                    for k = 1:length(obj.statistics{i}{j})
                        if ~obj.statistics{i}{j}(k).Inside && ~obj.statistics{i}{j}(k).Dead
                            if obj.statistics{i}{j}(k).isSecondary
                                obj.depthHistogramSE{i}(end+1) = obj.statistics{i}{j}(k).InitialDepth;
                            end
                        end
                    end
                end
            end
        end

        function plotDepthDistribution(obj,ind,nbins,varargin)
            if nargin < 3
                nbins = 200;
            end
            figure
            hold on
            box on
            histogram(obj.depthHistogramPE{ind},nbins,DisplayName='Primaries')
            histogram(obj.depthHistogramSE{ind},nbins,DisplayName='Secondaries')
            xlabel('Initial depth (A)')
            ylabel('Counts')
            fontsize(20,"points")
            title(strrep(obj.matName,'_',' '))
            legend
        end

        function calculateCoincidenceHistogram(obj,truePairs)
            if obj.onlyEscaped
                error('For coincidence histogram all electron trajectories must be stored.')
            end
            obj.coincidenceHistogram = cell(numel(obj.energyArray),1);
            for i = 1:numel(obj.energyArray)
                for j = 1:obj.numTrajectories
                    for k = 1:length(obj.statistics{i}{j})
                        if ~obj.statistics{i}{j}(k).Dead && obj.statistics{i}{j}(k).isSecondary
                            if ~obj.statistics{i}{j}(obj.statistics{i}{j}(k).ParentIndex).Dead
                                if ~truePairs
                                    obj.coincidenceHistogram{i}(end+1,:) = [obj.statistics{i}{j}(obj.statistics{i}{j}(k).ParentIndex).Energy obj.statistics{i}{j}(k).Energy];
                                elseif obj.statistics{i}{j}(k).Energy + obj.statistics{i}{j}(k).InnerPotential == ...
                                    obj.statistics{i}{j}(obj.statistics{i}{j}(k).ParentIndex).EnergyLoss + obj.statistics{i}{j}(obj.statistics{i}{j}(k).ParentIndex).EnergySE
                                    obj.coincidenceHistogram{i}(end+1,:) = [obj.statistics{i}{j}(obj.statistics{i}{j}(k).ParentIndex).Energy obj.statistics{i}{j}(k).Energy];
                                end
                            end
                        end
                    end
                end
            end
        end

        function plotCoincidenceHistogram(obj,ind,nbins)
            figure
            hold on
            box on
            histogram2(obj.coincidenceHistogram{ind}(:,1),obj.coincidenceHistogram{ind}(:,2),nbins,'FaceColor','flat')
            xlabel('Electron energy (eV)')
            ylabel('Electron energy (eV)')
            xlim([0,obj.energyArray(ind)])
            colormap('turbo')
            colorbar
            fontsize(20,"points")
            title(strrep(obj.matName,'_',' '))

	        savefig([obj.matName '_' num2str(obj.energyArray(ind)) '_e2e.fig'])
            savepdf(gcf,[obj.matName '_' num2str(obj.energyArray(ind)) '_e2e.pdf'])
        end

        function getTrajectories(obj)
            if ~obj.trackTrajectories
                error('The trajectories were not tracked.')
            end
            obj.coordinateArray = cell(numel(obj.energyArray));
            for i = 1:numel(obj.energyArray)
                temp = cell(1,1);
                for j = 1:obj.numTrajectories                 
                    for k = 1:length(obj.statistics{i}{j})
                        temp{end+1}(:,:) = [obj.statistics{i}{j}(k).Coordinates repmat(obj.statistics{i}{j}(k).Generation,size(obj.statistics{i}{j}(k).Coordinates,1),1)];
                    end
                end
                obj.coordinateArray{i} = temp;
            end
        end

        function plotTrajectories(obj,ind)
            figure
            title([strrep(obj.matName,'_',' ') ' E_0 = ' num2str(obj.energyArray(ind)) ' eV'])
            hold on
            box on
            colormap(turbo)
            for i = 2:length(obj.coordinateArray{ind})    
                patch([obj.coordinateArray{ind}{i}(:,1); nan],[obj.coordinateArray{ind}{i}(:,3); nan],[obj.coordinateArray{ind}{i}(:,4); nan],'FaceColor','none','EdgeColor','flat') %,'Marker','o','MarkerFaceColor','flat')
            end
            scatter(0,0,50,'filled','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0])
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

            figure
            title([strrep(obj.matName,'_',' ') ' E_0 = ' num2str(obj.energyArray(ind)) ' eV'])
            hold on
            box on
            for i = 2:length(obj.coordinateArray{ind})    
                patch([obj.coordinateArray{ind}{i}(:,1); nan],[obj.coordinateArray{ind}{i}(:,3); nan],[obj.coordinateArray{ind}{i}(:,5); nan],'FaceColor','none','EdgeColor','flat') %,'Marker','o','MarkerFaceColor','flat')
            end
            scatter(0,0,50,'filled','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0])
            set(gca,'Ydir','reverse')
            c = colorbar; 
            ylabel(c,'Electron generation');
            cm = lines(c.Limits(2));
            colormap(cm)
            xlabel('X (A)');
            ylabel('Z (A)')
            
            xl = xlim;
            yl = ylim;
            x = [xl(1) xl(2) xl(2) xl(1)];
            y = [0 0 yl(2) yl(2)];
            fill(x,y,'b','FaceAlpha',0.2)
            
            fontsize(20,'points')
        end

    end
end
