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
        Angles(1,2) double
        xyz(1,3) double {mustBeReal} = [0,0,0]
        uvw(1,3) double {mustBeReal} = [0,0,1]
        coordinates {mustBeReal}
        saveCoordinates = false
        isSecondary = false
        Generation = 0 % Primary electron
        ParentIndex = 0 % Primary electron
        Inside = true
        Dead = false
        ScatteringType
        nSecondaries = 0
        EnergySE
        ConductionBandreference = false
    end
    properties (Dependent)
        IIMFP
        IEMFP
        IPHMFP
        ITMFP
    end
    methods
        function obj = Electron(e,layers,cbRef,saveCoord,xyz,uvw,gen,se,ind)
            obj.Layers = layers;
            if length(layers) == 2
                obj.multiLayer = true;
            end
            obj.ConductionBandreference = cbRef;
            if obj.Layers(1).Material.isMetal
                obj.InnerPotential = obj.Layers(1).Material.MaterialData.Ef + obj.Layers(1).Material.MaterialData.Wf;
            else
                obj.InnerPotential = obj.Layers(1).Material.MaterialData.Affinity + obj.Layers(1).Material.MaterialData.Evb + obj.Layers(1).Material.MaterialData.Eg;
            end
            if nargin > 4
                obj.xyz = xyz;
                obj.uvw = uvw;
                obj.Generation = gen;
                obj.isSecondary = se;
                obj.ParentIndex = ind;
            end
            
            if ~obj.isSecondary
                obj.Energy = e + obj.InnerPotential;
            else
                obj.Energy = e;
            end

            if saveCoord
                obj.saveCoordinates = true;
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
            if obj.saveCoordinates
                obj.coordinates(end+1,:) = [obj.xyz,obj.Energy];
            end
        end
        function dircos2ang(obj)
            obj.Angles(1) = acos(obj.uvw(3));
            obj.Angles(2) = atan2(obj.uvw(2),obj.uvw(1));
        end
        function getScatteringType(obj)
            rn = rand;
            if rn < obj.IEMFP/obj.ITMFP
                obj.ScatteringType = 0;
            elseif rn < (obj.IEMFP+obj.IIMFP)/obj.ITMFP
                obj.ScatteringType = 1;
            else
                obj.ScatteringType = 2;
            end
        end
        function loss = scatter(obj)
            obj.Deflection(2) = rand*2*pi;
            if obj.ScatteringType == 0
                loss = false;
                decs = obj.Layers(obj.currentLayer).Material.getDECS(obj.Energy);
                cumdecs = cumtrapz(obj.Layers(obj.currentLayer).Material.MaterialData.DECS.x,decs);
                obj.Deflection(1) = interp1(cumdecs,obj.Layers(obj.currentLayer).Material.MaterialData.DECS.x,rand);
                obj.uvw = updateDirection(obj.uvw,obj.Deflection,1);
            elseif obj.ScatteringType == 1
                [eloss,diimfp] = obj.Layers(obj.currentLayer).Material.getDIIMFP(obj.Energy);
                cumdiimfp = cumtrapz(eloss,diimfp);
                cumdiimfp = (cumdiimfp - cumdiimfp(1))/(cumdiimfp(end)-cumdiimfp(1));
                while true
                    obj.EnergyLoss = interp1(cumdiimfp,eloss,rand);
                    if obj.EnergyLoss < obj.Energy
                        break;
                    end
                end
                loss = true;
                if obj.Layers(obj.currentLayer).Material.isMetal
                    obj.EnergySE = fegdos(obj.EnergyLoss,obj.Layers(obj.currentLayer).Material.MaterialData.Ef);
                else
                    obj.EnergySE = fegdos(obj.EnergyLoss,obj.Layers(obj.currentLayer).Material.MaterialData.Evb);
                end
                obj.Energy = obj.Energy - obj.EnergyLoss;
                obj.died;
                if obj.Layers(obj.currentLayer).Material.isMetal
                    min_energy = 1;
                else
                    min_energy = obj.Layers(obj.currentLayer).Material.MaterialData.Eg;
                end
                if ~obj.Dead && obj.Energy > min_energy
                    [theta, angdist] = obj.Layers(obj.currentLayer).Material.getAngularIIMFP(obj.Energy+obj.EnergyLoss,obj.EnergyLoss);
                    cumang = cumtrapz(theta, angdist);
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
            obj.dircos2ang;
            if ~obj.Dead
                if obj.xyz(end) < 0
                    beta = pi - obj.Angles(1);
                    if ~obj.Layers(1).Material.isMetal && obj.ConductionBandreference
                        ecos = (obj.Energy - obj.Layers(1).Material.MaterialData.Eg - obj.Layers(1).Material.MaterialData.Evb)*cos(beta)^2;
                        ui = obj.Layers(1).Material.MaterialData.Affinity;
                    else
                        ecos = obj.Energy*cos(beta)^2;
                        ui = obj.InnerPotential;
                    end
                    if ecos > ui
                        t = 4*sqrt(1 - ui/ecos)/(1 + sqrt(1 - ui/ecos))^2;
                    else
                        t = 0;
                    end
                    if rand < t
                        obj.Inside = false;
                        obj.xyz(1) = obj.xyz(1) + sin(beta)*cos(obj.Angles(2))*obj.xyz(end)/cos(beta);
                        obj.xyz(2) = obj.xyz(2) + sin(beta)*sin(obj.Angles(2))*obj.xyz(end)/cos(beta);
                        obj.xyz(3) = 0;
                        obj.Angles(1) = pi - asin(sin(beta)*sqrt(obj.Energy/(obj.Energy-obj.InnerPotential)));
                        if obj.saveCoordinates
                            obj.coordinates(end,:) = [obj.xyz,obj.Energy];
                        end
                        obj.Energy = obj.Energy - obj.InnerPotential;
                        if obj.saveCoordinates
                            x = obj.xyz(1) + 100*sin(obj.Angles(1))*cos(obj.Angles(2));
                            y = obj.xyz(2) + 100*sin(obj.Angles(1))*sin(obj.Angles(2));
                            z = obj.xyz(3) + 100*cos(obj.Angles(1));
                            obj.coordinates(end+1,:) = [x,y,z,obj.Energy];
                        end
                    else
                        obj.uvw(end) = -1*obj.uvw(end);
                        obj.xyz(end) = -1*obj.xyz(end);
                        if obj.saveCoordinates
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
            if obj.Layers(obj.currentLayer).Material.isMetal
                val = 1/obj.Layers(obj.currentLayer).Material.getEMFP(obj.Energy);
            else
                val = 1/obj.Layers(obj.currentLayer).Material.getEMFP(obj.Energy - ...
                    obj.Layers(obj.currentLayer).Material.MaterialData.Eg - obj.Layers(obj.currentLayer).Material.MaterialData.Evb);
            end
        end
        function val = get.IPHMFP(obj)
            if ~obj.Layers(obj.currentLayer).Material.isMetal
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