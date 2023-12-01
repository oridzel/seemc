function e = simulateSEE(N_traj,E0,MatName)

    e_count = 0;
    Mat = Sample(MatName,false);
    e = Electron.empty;
    
    for i = 1:N_traj
        e(end + 1) = Electron(E0,Mat);
        while e_count < length(e)
            e_count = e_count + 1;
            % disp(e_count);
            while e(e_count).Inside && ~e(e_count).Dead
                e(e_count).travel;
                e(e_count).escape;
                if e(e_count).Inside && ~e(e_count).Dead
                    e(e_count).getScatteringType;
                    if e(e_count).scatter
                        % if e(j).EnergyLoss + e(j).Material.MaterialData.Evb > e(j).InnerPotential
                        if e(e_count).EnergyLoss - e(e_count).Material.MaterialData.Eg - e(e_count).EnergySE > e(e_count).Material.MaterialData.Affinity
                            % e_se = e(j).EnergyLoss + e(j).Material.MaterialData.Evb;
                            e_se = e(e_count).EnergyLoss - e(e_count).Material.MaterialData.Eg - e(e_count).EnergySE;
                            uvw = [ sin(acos(2*rand-1))*cos(2*rand*pi),...
                                    sin(acos(2*rand-1))*sin(2*rand*pi),...
                                    cos(acos(2*rand-1)) ];
                            e(end + 1) = Electron(e_se,Mat,e(e_count).xyz,uvw,e(e_count).nSecondaries+1,true);
                            e(e_count).nSecondaries = e(e_count).nSecondaries + 1;
                        end
                    end
                end
            end
        end
    end

end

