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
            for e = 1:numel(obj.energy)
                disp(obj.energy(e))
                tic
                x = cell(1);              
                for i = 1:obj.numTrajectories
                    e_count = 0;
                    res = Electron.empty;
                    res(end+1) = Electron(obj.energy(e),obj.sample,obj.trackTrajectories);
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
                                        res(end + 1) = Electron(e_se,obj.sample,obj.trackTrajectories,res(e_count).xyz,uvw,res(e_count).nSecondaries+1,true,e_count);
                                        res(e_count).nSecondaries = res(e_count).nSecondaries + 1;
                                    end
                                end
                            end            
                        end
                        if ~res(e_count).Dead && ~res(e_count).Inside
                            x{end+1}.Energy = res(e_count).Energy;
                            x{end+1}.EnergyLoss = res(e_count).EnergyLoss;
                            x{end+1}.EnergySE = res(e_count).EnergySE;
                            x{end+1}.Angles = res(e_count).Angles;
                            x{end+1}.isSecondary = res(e_count).isSecondary;
                            x{end+1}.Generation = res(e_count).Generation;
                            x{end+1}.ParentIndex = res(e_count).ParentIndex;
                            if obj.trackTrajectories
                                x{end+1}.coordinates = res(e_count).coordinates;
                            end
                        end
                    end
                end
                obj.statistics{e} = x;
                toc
            end
        end
    end
end