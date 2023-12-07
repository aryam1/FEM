classdef MeshObj
    properties
        t {mustBeNumeric} = 0 % Current time for mesh
        elemN {mustBeNumeric} % Number of elements in mesh
        minX {mustBeGreaterThanOrEqual(minX,0)} % Starting x position of mesh
        maxX {mustBeGreaterThan(maxX,0)} % End x position of mesh
        lBoundary
        lBoundaryType
        rBoundary
        rBoundaryType
        basisType
    end
    properties (Dependent)
        dx % Spatial width of linear nodes
        L2 % L2 Norm error for the mesh
        globalN % Number of nodes in mesh
        nVec % Position of every global node
    end
    methods
        function xStep = get.dx(obj)
            xStep = (obj.maxX - obj.minX)/obj.elemN; % Computing linear spatial step
        end
        function gNodes = get.nVec(obj)
            gNodes = obj.minX:obj.dx:obj.maxX; % Change to Array of element objects
        end
        function error = get.L2(obj)
            error =1; % Change to sum e of each element
        end
        function nodes = get.globalN(obj)
            nodes = obj.elemN+1; % Number of nodes is 1 greater than elements
        end
        
    end
end