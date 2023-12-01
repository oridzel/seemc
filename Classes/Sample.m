classdef Sample
    properties
        Name
        isMetal    
    end
    properties
        MaterialData
    end
    methods
        function obj = Sample(MatName,isMetal)
            ps = inputParser;
            ps.FunctionName = 'Sample';
            MatNameValidation = @(x) (iscell(x) && isscalar(x) && ischar(x{1})) || ischar(x);
            ps.addRequired('MatName', MatNameValidation);
            ps.parse(MatName);

            obj.Name = char(ps.Results.MatName);
            obj.isMetal = isMetal;
        end

        function decs = getDECS(obj,energy)
            decs = interp2(obj.MaterialData.DECS.E0,obj.MaterialData.DECS.x,obj.MaterialData.DECS.y,energy,obj.MaterialData.DECS.x);
        end

        function [eloss,diimfp] = getDIIMFP(obj,energy)
            idx = interp1(obj.MaterialData.DECS.E0,1:length(obj.MaterialData.DECS.E0),energy,'nearest');
            if energy > obj.MaterialData.DECS.E0(idx)
                prob = log(energy/obj.MaterialData.DECS.E0(idx)) / log(obj.MaterialData.DECS.E0(idx+1)/obj.MaterialData.DECS.E0(idx));
                if rand < prob
                    idx = idx + 1;
                end
            elseif energy < obj.MaterialData.DECS.E0(idx)
                prob = log(obj.MaterialData.DECS.E0(idx)/energy) / log(obj.MaterialData.DECS.E0(idx)/obj.MaterialData.DECS.E0(idx-1));
                if rand < prob
                    idx = idx - 1;
                end
            end
            eloss = obj.MaterialData.DIIMFP(:,1,idx);
            diimfp = obj.MaterialData.DIIMFP(:,2,idx);
        end

        function [theta,dist] = getAngularIIMFP(obj,energy,eloss)
            theta = linspace(0,pi/2,100);
            energy = energy/h2ev;
            eloss = eloss/h2ev;
            q_squared = 4*(energy - obj.MaterialData.Eg/h2ev) - 2*eloss - ...
                4*sqrt((energy - obj.MaterialData.Eg/h2ev)*(energy - obj.MaterialData.Eg/h2ev - eloss))*cos(theta);
            elf_int = interp2(obj.MaterialData.q*a0,obj.MaterialData.eloss/h2ev,obj.MaterialData.ELF,sqrt(q_squared),eloss);
            dist = 1./(pi^2*q_squared).*sqrt(1 - eloss/(energy - obj.MaterialData.Eg/h2ev)).*elf_int;
        end

        function imfp = getIMFP(obj,energy)
            imfp = interp1(obj.MaterialData.DECS.E0,obj.MaterialData.l_in,energy);
        end

        function emfp = getEMFP(obj,energy)
            emfp = interp1(obj.MaterialData.DECS.E0,obj.MaterialData.Elastic.l_el,energy);
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