classdef Sample
    properties
        Name
    end
    properties (Constant)
        T = 300 % K
		k_B = 8.617e-5 % eV/K
    end    
    properties
        MaterialData
    end
    methods
        function obj = Sample(MatName)
            ps = inputParser;
            ps.FunctionName = 'Sample';
            MatNameValidation = @(x) (iscell(x) && isscalar(x) && ischar(x{1})) || ischar(x);
            ps.addRequired('MatName', MatNameValidation);
            ps.parse(MatName);

            obj.Name = char(ps.Results.MatName);
        end

        function decs = getDECS(obj,energy)
            decs = interp2(obj.MaterialData.DECS.E0,obj.MaterialData.DECS.x,obj.MaterialData.DECS.y,energy,obj.MaterialData.DECS.x);
        end

        function [eloss,diimfp] = getDIIMFP(obj,energy)
            diimfp = interp1(obj.MaterialData.DECS.E0,squeeze(obj.MaterialData.DIIMFP(:,2,:))',energy);
            if obj.MaterialData.isMetal
                energy = energy - obj.MaterialData.Ef;
                eloss = eps:(energy-eps)/(length(obj.MaterialData.DIIMFP(:,1,1))-1):energy;
            else
                energy = energy - obj.MaterialData.Eg - obj.MaterialData.Evb;
                eloss = obj.MaterialData.Eg:(energy-obj.MaterialData.Eg)/(length(obj.MaterialData.DIIMFP(:,1,1))-1):energy;
            end
        end

        function [theta,dist] = getAngularIIMFP(obj,energy,eloss)
            theta = linspace(0,pi/2,100);            
            eloss = eloss/h2ev;
            if obj.MaterialData.isMetal
                energy = energy/h2ev;
            else
                energy = (energy - obj.MaterialData.Eg)/h2ev;
            end
            q_squared = 4*energy - 2*eloss - 4*sqrt(energy*(energy - eloss))*cos(theta);
            elf_int = interp2(obj.MaterialData.q*a0,obj.MaterialData.eloss/h2ev,obj.MaterialData.ELF,sqrt(q_squared),eloss);
            dist = 1./(pi^2*q_squared).*sqrt(1 - eloss/energy).*elf_int;
        end

        function imfp = getIMFP(obj,energy)
            imfp = interp1(obj.MaterialData.DECS.E0,obj.MaterialData.l_in,energy);
        end

        function emfp = getEMFP(obj,energy)
            emfp = interp1(obj.MaterialData.DECS.E0,obj.MaterialData.Elastic.l_el,energy);
        end

        function iphmfp = getIPHMFP(obj,energy)
            if isfield(obj.MaterialData,'Phonon')
                de_over_e = obj.MaterialData.Phonon.eloss/energy;
                n_lo = 1./(exp(obj.MaterialData.Phonon.eloss/(obj.T*obj.k_B)) - 1);
                ln_plus = (1 + sqrt(abs(1 - de_over_e))) / (1 - sqrt(abs(1 - de_over_e)));
                iphmfp = 1/a0 * (1/obj.MaterialData.Phonon.eps_inf - 1/obj.MaterialData.Phonon.eps_zero) .* ( n_lo + 1 ) .* de_over_e/2 * log(ln_plus);
            else
                iphmfp = 0;
            end
        end

        function plotEnergyMomentumELF(obj)
            figure
            box on
            h = pcolor(obj.MaterialData.q,obj.MaterialData.eloss,obj.MaterialData.ELF);
            h.EdgeColor = 'None';
            xlabel("Momentum transfer (A^{-1})");
            ylabel("Energy loss (eV)")
            title(strrep(obj.Name,'_',' '))
            fontsize(18,"points")
        end

        function val = get.MaterialData(obj)
            val = MaterialDatabase.getData(obj.Name);            
        end
    end

    methods (Static = true)
        function MaterialList = GetMaterialList()
            Data = load(Sample.databaseFileName);            
            MaterialList = fields(Data.MaterialData);
        end
    end
end