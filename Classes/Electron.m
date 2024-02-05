classdef Electron < handle
    properties
        Energy double {mustBeNonnegative}
        EnergyLoss double {mustBeNonnegative}
        InnerPotential double {mustBeNonnegative}
        PathLength double {mustBeNonnegative} = 0
        Layers(1,:) Layer
        currentLayer = 1
        multiLayer = false
        Deflection(1,2) double {mustBeReal} = [0,0]
        xyz(1,3) double {mustBeReal} = [0,0,0]
        uvw(1,3) double {mustBeReal} = [0,0,1]
        coordinates {mustBeReal}
        trackCoordinates = false
        isSecondary = false
        Generation = 0 % Primary electron
        ParentIndex = 0 % Primary electron
        Inside = true
        Dead = false
        ScatteringType
        nSecondaries = 0
        EnergySE
        ConductionBandreference = false
        IncidentAngle = 0
    end
    properties (Dependent)
        IIMFP
        IEMFP
        IPHMFP
        ITMFP
    end
    methods
        function obj = Electron(energy,layers,varargin)
            p = inputParser;
            validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
            addRequired(p,'energy',validScalarPosNum);
            addRequired(p,'layers');
            addOptional(p,'theta',obj.IncidentAngle,@(x) isnumeric(x) && isscalar(x) && (x >= 0) && (x <= pi));
            addOptional(p,'cbRef',obj.ConductionBandreference);
            addOptional(p,'trackCoordinates',obj.trackCoordinates);
            addOptional(p,'xyz',obj.xyz);
            addOptional(p,'uvw',obj.uvw);
            addOptional(p,'generation',obj.Generation);
            addOptional(p,'isSecondary',obj.isSecondary);
            addOptional(p,'index',obj.ParentIndex);
            parse(p,energy,layers,varargin{:});

            obj.IncidentAngle = p.Results.theta;
            obj.Layers = p.Results.layers;
            if length(layers) == 2
                obj.multiLayer = true;
            end
            obj.ConductionBandreference = p.Results.cbRef;
            if obj.Layers(1).Material.MaterialData.isMetal
                obj.InnerPotential = obj.Layers(1).Material.MaterialData.Ef + obj.Layers(1).Material.MaterialData.Wf;
            else
                obj.InnerPotential = obj.Layers(1).Material.MaterialData.Affinity + obj.Layers(1).Material.MaterialData.Evb + obj.Layers(1).Material.MaterialData.Eg;
            end
            obj.isSecondary = p.Results.isSecondary;
            obj.xyz = p.Results.xyz;
            obj.Generation = p.Results.generation;
            obj.ParentIndex = p.Results.index;
            if obj.isSecondary
                obj.uvw = p.Results.uvw;
                obj.Energy = p.Results.energy;
            else
                obj.Energy = p.Results.energy + obj.InnerPotential;
                theta_0 = asin(sin(p.Results.theta)*sqrt(obj.Energy/(obj.Energy + obj.InnerPotential)));
                obj.uvw = [sin(theta_0),0,cos(theta_0)];
            end
            obj.trackCoordinates = p.Results.trackCoordinates;
            if obj.trackCoordinates
                obj.coordinates(1,:) = [obj.xyz,obj.Energy];
            end
        end
        function travel(obj)
            s = -(1/obj.ITMFP)*log(rand);

            if obj.xyz(3) + obj.uvw(3)*s < 0 % crosses the boundary solid -> vacuum
                s = abs(obj.xyz(3)/obj.uvw(3)) + 0.0001;
                obj.currentLayer = 1;
            elseif obj.multiLayer
                % check if crosses the boundary between layers
                for i = 1:length(obj.Layers)-1
                    if obj.xyz(3) < obj.Layers(i).Thickness && obj.xyz(3) + obj.uvw(3)*s > obj.Layers(i).Thickness
                        s = abs( (obj.Layers(i).Thickness - obj.xyz(3))/obj.uvw(3) ) + 0.0001;
                        obj.currentLayer = i + 1;
                        break
                    elseif obj.xyz(3) > obj.Layers(i).Thickness && obj.xyz(3) + obj.uvw(3)*s < obj.Layers(i).Thickness
                        s = abs( (obj.Layers(i).Thickness - obj.xyz(3))/obj.uvw(3) ) + 0.0001;
                        obj.currentLayer = i;
                        break
                    end
                end
            end

            obj.PathLength = obj.PathLength + s;
            obj.xyz(1) = obj.xyz(1) + obj.uvw(1)*s;
            obj.xyz(2) = obj.xyz(2) + obj.uvw(2)*s;
            obj.xyz(3) = obj.xyz(3) + obj.uvw(3)*s;
            if obj.trackCoordinates
                obj.coordinates(end+1,:) = [obj.xyz,obj.Energy];
            end
        end
        function getScatteringType(obj)
            rn = rand;
            if rn < obj.IEMFP/obj.ITMFP
                obj.ScatteringType = 0;
            elseif rn < (obj.IEMFP + obj.IIMFP)/obj.ITMFP
                obj.ScatteringType = 1;
            else
                obj.ScatteringType = 2;
            end
        end
        function loss = scatter(obj)
            obj.Deflection(2) = rand*2*pi;
            [~,energy_index] = min(abs(obj.Layers(obj.currentLayer).Material.MaterialData.DECS.E0 - obj.Energy));
            if obj.ScatteringType == 0
                loss = false;
                % decs = obj.Layers(obj.currentLayer).Material.getDECS(obj.Energy);
                decs = obj.Layers(obj.currentLayer).Material.MaterialData.DECS.y(:,energy_index);
                cumsigma = cumtrapz( obj.Layers(obj.currentLayer).Material.MaterialData.DECS.x,2*pi*decs.*sin(obj.Layers(obj.currentLayer).Material.MaterialData.DECS.x) );
                % cumsigma = (cumsigma - cumsigma(1))/(cumsigma(end) - cumsigma(1));
                obj.Deflection(1) = interp1(cumsigma,obj.Layers(obj.currentLayer).Material.MaterialData.DECS.x,rand*cumsigma(end));
                obj.uvw = updateDirection(obj.uvw,obj.Deflection,1);
            elseif obj.ScatteringType == 1
                % [eloss,diimfp] = obj.Layers(obj.currentLayer).Material.getDIIMFP(obj.Energy);
                eloss = obj.Layers(obj.currentLayer).Material.MaterialData.DIIMFP(:,1,energy_index);
                diimfp = obj.Layers(obj.currentLayer).Material.MaterialData.DIIMFP(:,2,energy_index);
                cumdiimfp = cumtrapz(eloss,diimfp);
                % cumdiimfp = (cumdiimfp - cumdiimfp(1))/(cumdiimfp(end)-cumdiimfp(1));
                 while true
                    obj.EnergyLoss = interp1(cumdiimfp,eloss,rand*cumdiimfp(end));
                    if obj.EnergyLoss < obj.Energy
                        break;
                    end
                end
                loss = true;
                if obj.Layers(obj.currentLayer).Material.MaterialData.isMetal
                    obj.EnergySE = fegdos(obj.EnergyLoss,obj.Layers(obj.currentLayer).Material.MaterialData.Ef);
                else
                    obj.EnergySE = fegdos(obj.EnergyLoss,obj.Layers(obj.currentLayer).Material.MaterialData.Evb);
                end
                obj.Energy = obj.Energy - obj.EnergyLoss;
                obj.died;
                if obj.Layers(obj.currentLayer).Material.MaterialData.isMetal
                    min_energy = 1;
                else
                    min_energy = obj.Layers(obj.currentLayer).Material.MaterialData.Eg;
                end
                if ~obj.Dead && obj.Energy > min_energy
                    [theta, angdist] = obj.Layers(obj.currentLayer).Material.getAngularIIMFP(obj.Energy+obj.EnergyLoss,obj.EnergyLoss);
                    cumang = cumtrapz(theta, angdist.*sin(theta));
                    cumang = (cumang - cumang(1))/(cumang(end)-cumang(1));
                    if isfinite(cumang)
                        obj.Deflection(1) = interp1(cumang,theta,rand);
                        obj.uvw = updateDirection(obj.uvw,obj.Deflection,1);
                    end
                end
            else
                rn = rand;
                e = (obj.Energy - obj.Layers(obj.currentLayer).Material.MaterialData.Eg - obj.Layers(obj.currentLayer).Material.MaterialData.Evb)/h2ev;
                de = obj.Layers(obj.currentLayer).Material.MaterialData.Phonon.eloss/h2ev;
                if e - de > 0
                    bph = (e + e - de + 2*sqrt(e*(e - de))) / (e + e - de - 2*sqrt(e*(e - de)));
                    obj.Deflection(1) = acos( (e + e - de)/(2*sqrt(e*(e - de)))*(1 - bph^rn) + bph^rn );
                    obj.Energy = obj.Energy - obj.Layers(obj.currentLayer).Material.MaterialData.Phonon.eloss;
                    obj.died;
                    if ~obj.Dead && isreal(obj.Deflection(1))
                        obj.uvw = updateDirection(obj.uvw,obj.Deflection,1);
                    end
                    loss = false;
                else
                    obj.Dead = true;
                end
            end
        end
        function escape(obj)
            if ~obj.Dead
                if obj.xyz(end) < 0
                    if ~obj.Layers(1).Material.MaterialData.isMetal && obj.ConductionBandreference
                        ecos = (obj.Energy - obj.Layers(1).Material.MaterialData.Eg - obj.Layers(1).Material.MaterialData.Evb)*obj.uvw(3)^2;
                        ui = obj.Layers(1).Material.MaterialData.Affinity;
                    else
                        ecos = obj.Energy*obj.uvw(3)^2;
                        ui = obj.InnerPotential;
                    end
                    if ecos > ui
                        t = 4*sqrt(1 - ui/ecos)/(1 + sqrt(1 - ui/ecos))^2;
                    else
                        t = 0;
                    end
                    if rand < t
                        obj.Inside = false;
                        old_xyd = sqrt(1 - obj.uvw(3)^2);
                        obj.uvw(3) = -sqrt( (ecos - ui)/(obj.Energy - obj.InnerPotential) );
                        xyd = sqrt(1 - obj.uvw(3)^2);
                        if abs(obj.uvw(3)) ~= 1
                            obj.uvw(1) = obj.uvw(1)*xyd/old_xyd;
                            obj.uvw(2) = obj.uvw(2)*xyd/old_xyd;
                        end
                        obj.Energy = obj.Energy - obj.InnerPotential;
                        if obj.trackCoordinates
                            x = obj.xyz(1) + 100*obj.uvw(1);
                            y = obj.xyz(2) + 100*obj.uvw(2);
                            z = obj.xyz(3) + 100*obj.uvw(3);
                            obj.coordinates(end+1,:) = [x,y,z,obj.Energy];
                        end
                    else
                        obj.uvw(end) = -1*obj.uvw(end);
                        obj.xyz(end) = -1*obj.xyz(end);
                        if obj.trackCoordinates
                            obj.coordinates(end,:) = [obj.xyz(1),obj.xyz(2),-obj.xyz(3),obj.Energy];
                        end
                    end
                end
            end
        end
        function died(obj)
            if obj.Energy < obj.InnerPotential
                obj.Dead = true;
            end
        end
        function testDECSsampling(obj)
            n = 10000;
            angle = zeros(n);
            decs = obj.Layers(obj.currentLayer).Material.getDECS(obj.Energy);
            cumdecs = cumtrapz(obj.Layers(obj.currentLayer).Material.MaterialData.DECS.x,decs);
            for i = 1:n
                angle(i) = interp1(cumdecs,obj.Layers(obj.currentLayer).Material.MaterialData.DECS.x,rand);
            end
            figure
            hold on
            h = histogram(angle);
            plot(obj.Layers(obj.currentLayer).Material.MaterialData.DECS.x,decs/max(decs)*max(h.Values))
            xlabel('Angle (rad)')
        end
        function testDIIMFPsampling(obj)
            n = 10000;
            loss = zeros(n);
            [eloss,diimfp] = obj.Layers(obj.currentLayer).Material.getDIIMFP(obj.Energy);
            cumdiimfp = cumtrapz(eloss,diimfp);
            cumdiimfp = (cumdiimfp - cumdiimfp(1))/(cumdiimfp(end)-cumdiimfp(1));
            for i = 1:n
                loss(i) = interp1(cumdiimfp,eloss,rand);
            end
            figure
            hold on
            h = histogram(loss);
            plot(eloss,diimfp/max(diimfp)*max(h.Values))
            xlabel('Energy loss (eV)')
        end
        function testAngIIMFPsampling(obj)
            n = 10000;
            ang = zeros(n);
            for i = 1:n
                [theta, angdist] = obj.Layers(obj.currentLayer).Material.getAngularIIMFP(500,20);
                cumang = cumtrapz(theta, angdist);
                cumang = (cumang - cumang(1))/(cumang(end)-cumang(1));
                ang(i) = interp1(cumang,theta,rand);
            end
            figure
            hold on
            h = histogram(ang);
            plot(theta, angdist/max(angdist)*max(h.Values))
            xlabel('Theta (rad)')
        end
        function set.Energy(obj,val)
            obj.Energy = val;
        end
        function val = get.IIMFP(obj)
            val = 1/obj.Layers(obj.currentLayer).Material.getIMFP(obj.Energy);
        end
        function val = get.IEMFP(obj)
            if obj.Layers(obj.currentLayer).Material.MaterialData.isMetal
                val = 1/obj.Layers(obj.currentLayer).Material.getEMFP(obj.Energy);
            else
                val = 1/obj.Layers(obj.currentLayer).Material.getEMFP(obj.Energy - ...
                    obj.Layers(obj.currentLayer).Material.MaterialData.Eg - obj.Layers(obj.currentLayer).Material.MaterialData.Evb);
            end
        end
        function val = get.IPHMFP(obj)
            if ~obj.Layers(obj.currentLayer).Material.MaterialData.isMetal
                val = obj.Layers(obj.currentLayer).Material.getIPHMFP(obj.Energy - ...
                    obj.Layers(obj.currentLayer).Material.MaterialData.Eg - obj.Layers(obj.currentLayer).Material.MaterialData.Evb);
            else
                val = 0;
            end
        end
        function val = get.ITMFP(obj)
            val = obj.IEMFP + obj.IIMFP + obj.IPHMFP;
        end
    end
end