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
    end
    methods
        function obj = SEEMC(inputpar)
            obj.matName = inputpar.matName;
            obj.isMetal = inputpar.isMetal;
            obj.numTrajectories = inputpar.numTrajectories;
            obj.trackTrajectories = inputpar.trackTrajectories;
            obj.energyArray = inputpar.energy;

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
    end
end