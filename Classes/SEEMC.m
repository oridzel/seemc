classdef SEEMC < handle
    properties
        numTrajectories = 1000
        matName = 'Au'
        isMetal = true
        trackTrajectories = false
        energy
        sample
        statistics
    end
    methods
        function obj = SEEMC(inputpar)
            obj.matName = inputpar.matName;
            obj.isMetal = inputpar.isMetal;
            obj.numTrajectories = inputpar.numTrajectories;
            obj.trackTrajectories = inputpar.trackTrajectories;
            obj.energy = inputpar.energy;

            obj.sample = Sample(obj.matName,obj.isMetal);
        end
        function simulate(obj)
            obj.statistics = cell(1);
            n = obj.numTrajectories;
            energy_array = obj.energy;
            s = obj.sample;
            track = obj.trackTrajectories;
            stat = cell(numel(energy_array));
            parfor e = 1:numel(energy_array)
                disp(energy_array(e))
                x = cell(1);              
                for i = 1:n
                    e_count = 0;
                    res = Electron.empty;
                    res(end+1) = Electron(energy_array(e),s,track);
                    while e_count < length(res)
                        y = struct;
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
                                        res(end + 1) = Electron(e_se,s,track,res(e_count).xyz,uvw,res(e_count).nSecondaries+1,true,e_count);
                                        res(e_count).nSecondaries = res(e_count).nSecondaries + 1;
                                    end
                                end
                            end            
                        end
                        if ~res(e_count).Dead && ~res(e_count).Inside
                            y.Energy = res(e_count).Energy;
                            y.EnergyLoss = res(e_count).EnergyLoss;
                            y.EnergySE = res(e_count).EnergySE;
                            y.Angles = res(e_count).Angles;
                            y.isSecondary = res(e_count).isSecondary;
                            y.Generation = res(e_count).Generation;
                            y.ParentIndex = res(e_count).ParentIndex;
                            if track
                                y.coordinates = res(e_count).coordinates;
                            end
                            x{end+1} = y;
                        end
                    end
                end
                stat{e} = x;
            end
            obj.statistics = stat;
        end
    end
end