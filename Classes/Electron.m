classdef Electron < handle
    properties
        Energy double {mustBeNonnegative}
        EnergyLoss double {mustBeNonnegative}
        InnerPotential double {mustBeNonnegative}
        PathLength double {mustBeNonnegative} = 0
        Material Sample
        Deflection(1,2) double = [0,0]
        Angles(1,2) double
        xyz(1,3) double = [0,0,0]
        uvw(1,3) double = [0,0,1]
        coordinates
        saveCoordinates = false
        isSecondary = false
        Generation = 1 % Primary electron
        ParentIndex = 1 % Primary electron
        Inside = true
        Dead = false
        ScatteringType
        nSecondaries = 0
        EnergySE
    end
    properties (Dependent)
        IIMFP
        IEMFP
        IPHMFP
        ITMFP
    end
    methods
        function obj = Electron(e,mat,saveCoord,xyz,uvw,gen,se,ind)
            obj.Material = mat;
            if obj.Material.isMetal
                obj.InnerPotential = obj.Material.MaterialData.Ef + obj.Material.MaterialData.Wf;
            else
                obj.InnerPotential = obj.Material.MaterialData.Affinity + obj.Material.MaterialData.Evb + obj.Material.MaterialData.Eg;
            end
            if nargin > 3
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
            if obj.xyz(3) + obj.uvw(3)*s < 0
                s = abs(obj.xyz(3)/obj.uvw(3))+0.0001;
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
        function updateDirection(obj)
            theta0 = acos(obj.uvw(end));
            phi0 = atan2(obj.uvw(2),obj.uvw(1));
            theta = acos(cos(theta0)*cos(obj.Deflection(1)) - sin(theta0)*sin(obj.Deflection(1))*cos(obj.Deflection(2)));
            phi = asin( sin(obj.Deflection(1))*sin(obj.Deflection(2))/sin(theta) ) + phi0;

            obj.uvw(1) = sin(theta)*cos(phi);
            obj.uvw(2) = sin(theta)*sin(phi);
            obj.uvw(3) = cos(theta);
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
                decs = obj.Material.getDECS(obj.Energy);
                cumdecs = cumtrapz(obj.Material.MaterialData.DECS.x,decs);
                obj.Deflection(1) = interp1(cumdecs,obj.Material.MaterialData.DECS.x,rand);
                obj.updateDirection;
            elseif obj.ScatteringType == 1
                [eloss,diimfp] = obj.Material.getDIIMFP(obj.Energy);
                cumdiimfp = cumtrapz(eloss,diimfp);
                cumdiimfp = (cumdiimfp - cumdiimfp(1))/(cumdiimfp(end)-cumdiimfp(1));
                while true
                    obj.EnergyLoss = interp1(cumdiimfp,eloss,rand);
                    if obj.EnergyLoss < obj.Energy
                        break;
                    end
                end
                loss = true;
                if obj.Material.isMetal
                    obj.EnergySE = fegdos(obj.EnergyLoss,obj.Material.MaterialData.Ef);
                else
                    obj.EnergySE = fegdos(obj.EnergyLoss,obj.Material.MaterialData.Evb);
                end
                obj.Energy = obj.Energy - obj.EnergyLoss;
                obj.died;
                if obj.Material.isMetal
                    min_energy = 1;
                else
                    min_energy = obj.Material.MaterialData.Eg;
                end
                if ~obj.Dead && obj.Energy > min_energy
                    [theta, angdist] = obj.Material.getAngularIIMFP(obj.Energy+obj.EnergyLoss,obj.EnergyLoss);
                    cumang = cumtrapz(theta, angdist);
                    cumang = (cumang - cumang(1))/(cumang(end)-cumang(1));
                    if isfinite(cumang)
                        obj.Deflection(1) = interp1(cumang,theta,rand);
                        obj.updateDirection;
                    end
                end
            else
                rn = rand;
                e = obj.Energy/h2ev;
                de = obj.Material.MaterialData.Phonon.eloss/h2ev;
                bph = (e + e - de + 2*sqrt(e*(e - de))) / (e + e - de - 2*sqrt(e*(e - de)));
                obj.Deflection(1) = acos( (e + e - de)/(2*sqrt(e*(e - de)))*(1 - bph^rn) + bph^rn );
                obj.Energy = obj.Energy - obj.Material.MaterialData.Phonon.eloss;
                obj.died;
                obj.updateDirection;
                loss = false;
            end
        end
        function escape(obj)
            obj.dircos2ang;
            if obj.xyz(end) < 0
                beta = pi - obj.Angles(1);
                ecos = obj.Energy*cos(beta)^2;             
                if ecos > obj.InnerPotential
                    t = 4*sqrt(1 - obj.InnerPotential/ecos)/(1 + sqrt(1 - obj.InnerPotential/ecos))^2;
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
        function died(obj)
            if obj.Energy < obj.InnerPotential
                obj.Dead = true;
            end
        end
        function testDECSsampling(obj)
            decs = obj.Material.getDECS(obj.Energy);
            cumdecs = cumtrapz(obj.Material.MaterialData.DECS.x,decs);
            for i = 1:100000
                angle(i) = interp1(cumdecs,obj.Material.MaterialData.DECS.x,rand);
            end
            figure
            hold on
            h = histogram(angle);
            plot(obj.Material.MaterialData.DECS.x,decs/max(decs)*max(h.Values))
            xlabel('Angle (rad)')
        end
        function testDIIMFPsampling(obj)
            [eloss,diimfp] = obj.Material.getDIIMFP(obj.Energy);
            cumdiimfp = cumtrapz(eloss,diimfp);
            cumdiimfp = (cumdiimfp - cumdiimfp(1))/(cumdiimfp(end)-cumdiimfp(1));
            for i = 1:100000
                loss(i) = interp1(cumdiimfp,eloss,rand);
            end
            figure
            hold on
            h = histogram(loss);
            plot(eloss,diimfp/max(diimfp)*max(h.Values))
            xlabel('Energy loss (eV)')
        end
        function testAngIIMFPsampling(obj)
            for i = 1:10000
                [theta, angdist] = obj.Material.getAngularIIMFP(500,20);
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
            val = 1/obj.Material.getIMFP(obj.Energy);
        end
        function val = get.IEMFP(obj)
            if obj.Material.isMetal
                val = 1/obj.Material.getEMFP(obj.Energy);
            else
                val = 1/obj.Material.getEMFP(obj.Energy - obj.Material.MaterialData.Eg - obj.Material.MaterialData.Evb);
            end
        end
        function val = get.IPHMFP(obj)
            if ~obj.Material.isMetal
                val = obj.Material.getIPHMFP(obj.Energy - obj.Material.MaterialData.Eg - obj.Material.MaterialData.Evb);
            else
                val = 0;
            end
        end
        function val = get.ITMFP(obj)
            val = obj.IEMFP + obj.IIMFP + obj.IPHMFP;
        end
    end
end