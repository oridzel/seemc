function e = simulateSEE(N_traj,E0,MatName,isMetal,trackTrajectories)

    Mat = Sample(MatName,isMetal);
    e = cell(N_traj,1);
    % parpool(2)

    for i = 1:N_traj
        e{i}(1) = Electron(E0,Mat,trackTrajectories);
        e_count = 0;
        while e_count < length(e{i})

            e_count = e_count + 1;
            while e{i}(e_count).Inside && ~e{i}(e_count).Dead
                e{i}(e_count).travel;
                e{i}(e_count).escape;
                if e{i}(e_count).Inside && ~e{i}(e_count).Dead
                    e{i}(e_count).getScatteringType;
                    if e{i}(e_count).scatter
                        if e{i}(e_count).EnergyLoss + e{i}(e_count).EnergySE > e{i}(e_count).InnerPotential
                            e_se = e{i}(e_count).EnergyLoss + e{i}(e_count).EnergySE;
                            uvw = [ sin(acos(2*rand-1))*cos(2*rand*pi),...
                                    sin(acos(2*rand-1))*sin(2*rand*pi),...
                                    cos(acos(2*rand-1)) ];
                            e{i}(end + 1) = Electron(e_se,Mat,trackTrajectories,e{i}(e_count).xyz,uvw,e{i}(e_count).nSecondaries+1,true);
                            e{i}(e_count).nSecondaries = e{i}(e_count).nSecondaries + 1;
                        end
                    end
                end

            end
        end
    end

end

