classdef MeshObj
    properties
        t {mustBeNumeric} % Current time
        dt {mustBeNumeric} % Time step
        elemN {mustBeNumeric} % Number of elements
        lBoundary % Left boundary value
        lBoundaryType % Left boundary type (d or n)
        rBoundary % Right boundary value
        rBoundaryType % Right boundary type (d or n)
        basisType {mustBeNumeric} % Basis function order
        Dvec (2,:) double % spatially varying diffusion term
        Lvec (2,:) double % spatially varying reaction term
        Fvec (2,:) double % spatially varying source term
        NBc % Neumann boundary conditions
        elems = ElemObj.empty % Array of element objects
        nVec % Position of every global node
        transient = false % Transient flag for time varying parameters
        KGlobal % Global stiffness matrix
        MGlobal % Global mass matrix
        FGlobal % Global source vector
        solution % Solution vector
        L2 % L2 Norm Error of mesh
    end
    properties (Dependent)
        globalN % Number of nodes in mesh 
    end
    methods
        function obj = MeshObj(x0,elemN, elemDist, boundaries, basis, t, dt) % Constructor
            scalar = round(elemN / sum(elemDist(2,:))); % Scaling factor to ensure correct number of elements
            obj.nVec = x0; % Initialize position vector with starting x 
            for i = 1:size(elemDist,2) % Loop through element distributions
                xpoints = linspace(obj.nVec(end),elemDist(1,i),basis*scalar*elemDist(2,i)+1); % Create equally spaced points, scaled by basis order and distribution
                obj.nVec = [obj.nVec,xpoints(2:end)]; % Append points to position vector
            end
            obj.elemN = (length(obj.nVec)-1)/basis; % Calculate true number of elements 

            [obj.lBoundary,obj.lBoundaryType,obj.rBoundary,obj.rBoundaryType] = boundaries{:}; % Unpack boundary conditions
            obj.NBc = zeros(obj.globalN,1); % Initialize Neumann boundary conditions
            obj.NBc(1) = obj.lBoundary * (obj.lBoundaryType == "n"); % Set Neumann boundary conditions
            obj.NBc(end) = obj.rBoundary * (obj.rBoundaryType == "n"); 

            obj.basisType = basis; % Set basis function order
            obj.t = t; % Set current time
            obj.dt = dt; % Set time step
        end

        function self = ElemSetup(self) % Create element objects
            for i = 1:self.elemN % Loop through elements
                I = self.basisType * (i-1) + 1; % Calculate insertion index for all basis types
                self.elems(i) = ElemObj(i,self.nVec(I:I+self.basisType),self.basisType); % Create element object
            end
        end

        function self = SetParams(self,dVec,lVec,fVec) % Set spatially varying parameters
            if any([~isnumeric(dVec),~isnumeric(lVec),~isnumeric(fVec)]) % Check if parameters are numeric arrays
                causeException = MException('MeshObj:SetParams',"Invalid parameters");
                throw(causeException);
            end
            self = self.ElemSetup(); % Create element objects

            self.Dvec = dVec;
            self.Lvec = lVec;
            self.Fvec = fVec;

            dPointer = 1; % Initialize pointers for parameter distributions
            lPointer = 1; 
            fPointer = 1;
            for i = 1:self.elemN % Loop through elements
                for j = 1:self.basisType+1 % Loop through nodes
                    while self.elems(i).x(j) > dVec(1,dPointer) % Find parameter values at node
                        dPointer = dPointer + 1; 
                    end
                    self.elems(i).D(j) = dVec(2,dPointer); % Set parameter values

                    while self.elems(i).x(j) > lVec(1,lPointer) % Find parameter values at node
                        lPointer = lPointer + 1;
                    end
                    self.elems(i).L(j) = lVec(2,lPointer); % Set parameter values

                    while self.elems(i).x(j) > fVec(1,fPointer) % Find parameter values at node
                        fPointer = fPointer + 1;
                    end
                    self.elems(i).F(j) = fVec(2,fPointer); % Set parameter values
                end
            end
        end

        function self = GlobalSetup(self) % Create global matrices
            KG = zeros(self.globalN); % Initialize global matrices
            MG = zeros(self.globalN); 
            FG = zeros(self.globalN,1);

            [psi,dpsi] = self.BasisFunc(self.basisType); % Get basis functions
            [p,w]= self.GQScheme(4); % Get quadrature scheme

            for i = 1:self.elemN % Loop through elements
                [K,M,F] = self.elems(i).LocalMatrix(psi,dpsi,p,w); % Get local matrices
                I = self.basisType*(i-1)+1; % Calculate insertion index for all basis types
                I2 = I+self.basisType; 

                KG(I:I2,I:I2) = KG(I:I2,I:I2) + K; % Insert local matrices into global matrices
                MG(I:I2,I:I2) = MG(I:I2,I:I2) + M;

                FG(I:I2,1) = FG(I:I2,1) + F; % Insert local source vector into global source vector
            end
            
            self.KGlobal = KG; % Set global matrices
            self.MGlobal = MG;
            self.FGlobal = FG;
        end

        function self = Solve(self,theta,previous) % Solve system
            if ~ismember(theta,[0 0.5 1]) % Check if theta is valid
                causeException = MException('MeshObj:Solve',"Invalid theta scheme");
                throw(causeException);
            end
            if ~isa(previous,'MeshObj') % Check if previous mesh is valid
                causeException = MException('MeshObj:Solve',"Invalid mesh object");
                throw(causeException);
            end
            % If not transient, copy global matrices, 
            % elements and materials distributions from previous mesh
            if ~self.transient 
                self.KGlobal = previous.KGlobal;
                self.MGlobal = previous.MGlobal;
                self.FGlobal = previous.FGlobal;
                self.elems = previous.elems;
                self.Dvec = previous.Dvec;
                self.Lvec = previous.Lvec;
                self.Fvec = previous.Fvec;
            else
                self = self.GlobalSetup(); % Otherwise, create global matrices
            end

            res = (self.MGlobal-(1-theta).*self.dt.*self.KGlobal)*previous.solution; % Calculate residual vector
            res = res + self.dt*theta*(self.FGlobal + self.NBc); % Add source vector and Neumann boundary conditions
            res = res + self.dt*(1-theta)*(previous.FGlobal + previous.NBc); % Add source vector and Neumann boundary conditions from previous time step

            gm = self.MGlobal+(theta.*self.dt.*self.KGlobal); % Calculate global matrix

            if self.lBoundaryType == "d" % Apply Dirichlet boundary conditions on left
                gm(1,:) = [1,zeros(1,self.globalN-1)];
                res(1)=self.lBoundary;
            end
            if self.rBoundaryType == "d" % Apply Dirichlet boundary conditions on right
                gm(end,:) = [zeros(1,self.globalN-1),1];
                res(end)=self.rBoundary;
            end

            self.solution = gm\res; % Solve system

        end

        function self = L2Norm(self)
            if any([self.Dvec ~= [1;1], self.Lvec ~= [1;0], self.Fvec ~= [1;0]]) % Check if analytical solution is known
                causeException = MException('MeshObj:L2Norm',"Analytical Solution Not Known For This System");
                throw(causeException);
            end
            [p,w] = self.GQScheme(4); % Get quadrature scheme
            [psi,~] = self.BasisFunc(self.basisType); % Get basis functions
            xVals = zeros(self.elemN,self.basisType+1);
            solVals = zeros(self.elemN,self.basisType+1);
            Jvec = zeros(self.elemN,1);
            
            for i = 1:self.elemN
                I = self.basisType * (i-1) + 1;
                I2 = I + self.basisType;
                solVals(i,:) = self.solution(I:I2); % Get solution values
                xVals(i,:) = self.nVec(I:I2); % Get x values
                Jvec(i) = self.elems(i).J; % Get Jacobian
            end
            E = zeros(self.elemN,4);
            for g = 1:4 % Loop through quadrature points
                xInterp = sum(xVals.*(cellfun(@(fun) fun(p(g)), psi)),2); % Interpolate x values
                analytic = self.TransientAnalyticalSol(xInterp,self.t); % Get analytical solution
                solInterp = sum(solVals.*(cellfun(@(fun) fun(p(g)), psi)),2); % Interpolate solution
                E(:,g) = w(g)*Jvec.*(analytic-solInterp).^2; % Calculate L2 norm error
            end
            self.L2 = sum(E,2); % Sum L2 norm error
        end

        function nodes = get.globalN(obj) % Get number of nodes
            nodes = (obj.elemN*obj.basisType)+1; % Number of nodes is 1 greater than elements times basis order
        end

        function ind = xInd(self,x)
            ind = find(self.nVec >= x, 1); % Find index of x value
        end

    end

    methods (Static)
        function [psi,dpsi] = BasisFunc(basis) % Get basis functions
            funcs = {
                {@(x)(1-x)/2, @(x)(1+x)/2} ...
                {@(x)-0.5,@(x)0.5};
                {@(x)0.5.*x.*(x-1), @(x)1-x.^2, @(x)0.5.*x.*(1+x)} ...
                {@(x)x-0.5,@(x)-2*x,@(x)x+0.5}
                }; % Cell array of basis functions
            try
                psi = funcs{basis}; % Get basis functions
                dpsi = funcs{basis+2}; % Get basis function derivatives
            catch
                causeException = MException('MeshObj:BasisFunc',"Invalid basis function order");
                throw(causeException);
            end
        end

        function [p,w] = GQScheme(order) % Get Gauss quadrature scheme
            data = {
                [0; 2]
                [-1/sqrt(3) 1/sqrt(3); 1 1]
                [-sqrt(3/5) 0 sqrt(3/5); 5/9 8/9 5/9]
                [-sqrt(3/7 + (2/7*sqrt(1.2))) ...
                -sqrt(3/7 - (2/7*sqrt(1.2))) ...
                sqrt(3/7 - (2/7*sqrt(1.2))) ...
                sqrt(3/7 + (2/7*sqrt(1.2)));
                (18-sqrt(30))/36 (18+sqrt(30))/36 ...
                (18+sqrt(30))/36 (18-sqrt(30))/36]
                }; % Cell array of Gauss quadrature schemes

            try
                gq = data{order}; % Get Gauss quadrature scheme
                p= gq(1,:); % Get quadrature points
                w= gq(2,:); % Get quadrature weights
            catch
                causeException = MException('MeshObj:GQ',"Invalid quadrature order");
                throw(causeException);
            end
        end

        function c = TransientAnalyticalSol(xVec,t)
            trans = 0;
            for k=1:1000
                trans = trans + ((((-1)^k)/k) * exp(-k^2*pi^2*t)*sin(xVec.*k*pi)); 
            end
            c = xVec + (2/pi)*trans;
        end

    end

end