classdef ElemObj
    properties
        D (1,:) double = [1]
        L (1,:) double = [0]
        F (1,:) double = [0]
        sol {mustBeNumeric} = 0
        x (1,2) double
        id {mustBeNumeric} = 0
        t {mustBeNumeric} = 0 % Current time for mesh
    end
    properties (Dependent)
        E % error
    end
    methods
        function obj = ElemObj(id,x0,x1,basis)
            obj.id = id;
            switch basis
                case 1
                    obj.x = [x0 x1];
                case 2
                    obj.x = [x0 (x0+x1)/2 x1];
            end
        end
        
    end
end