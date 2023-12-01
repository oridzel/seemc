classdef Barrier
    properties
        Height = 1;
        Width = 0;
    end
    properties (Dependent)
        Transmission
    end
    methods
        function obj = Barrier(h,w)
            if nargin > 1
                obj.Width = w;
            else
                obj.Height = h;
                obj.Width = 0;
            end
        end
        function obj = set.Height(obj,h)
            if h < 0
                error('Barrier height must be positive')
            end
            obj.Height = h;
        end
        function val = get.Transmission(obj)
            val = 1;
        end
    end
    methods (Static)
        function disp(obj)
            h = obj.Height;
            disp(['Barrier with height: ',num2str(h)])
        end
    end
end