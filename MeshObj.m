classdef MeshObj
    properties
        t {mustBeNumeric} % Current time for mesh
        dt {mustBeNumeric}
        elemN {mustBeNumeric} % Number of elements in mesh
        minX {mustBeGreaterThanOrEqual(minX,0)} % Starting x position of mesh
        maxX {mustBeGreaterThan(maxX,0)} % End x position of mesh
        lBoundary
        lBoundaryType
        rBoundary
        rBoundaryType
        basisType
        Dvec (2,:) double % spatially varying diffusion term 
        Lvec (2,:) double % spatially varying reaction term
        Fvec (1,:) double % spatially varying source term
        elems = ElemObj.empty
    end
    properties (Dependent)
        dx % Spatial width of linear nodes
        L2 % L2 Norm error for the mesh
        globalN % Number of nodes in mesh
        nVec % Position of every global node
    end
    methods
        function obj = MeshObj(minX, maxX, elemN, boundaries, basis, t)
             obj.minX = minX;
             obj.maxX = maxX;
             obj.elemN = elemN;
             
             [obj.lBoundary,obj.lBoundaryType,obj.rBoundary,obj.rBoundaryType] = boundaries{:};
             obj.basisType = basis;
             obj.t = t;
             obj = obj.elemSetup();
        end

        function obj = elemSetup(obj)
            for i = 1:obj.elemN
                obj.elems(i) = ElemObj(i,obj.nVec(i),obj.nVec(i+1),obj.basisType);
            end
        end

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
            nodes = (obj.elemN*obj.basisType)+1; % Number of nodes is 1 greater than elements times basis order
        end
        
    end
end