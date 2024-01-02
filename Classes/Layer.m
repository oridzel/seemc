classdef Layer < handle
    
    properties (SetObservable = true)
        Material Sample
        Thickness double {mustBeNonnegative} = Inf
    end
    
    properties (Access = private, Constant = true)
        defaultLayerThickness = Inf;
    end
    
    methods
        function obj = Layer(Material,thickness,varargin)
            obj.Material = Material;
            if nargin == 2
                obj.Thickness = thickness;
            end
        end       
        function set.Thickness(obj, val)
            obj.Thickness = val;
        end
    end
    
end

