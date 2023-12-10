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
        Fvec (2,:) double % spatially varying source term
        NBc
        elems = ElemObj.empty
        nVec % Position of every global node
        transient = false
        KGlobal
        MGlobal
        FGlobal
        solution
    end
    properties (Dependent)
        globalN % Number of nodes in mesh
    end
    methods
        function obj = MeshObj(minX, maxX, elemN, boundaries, basis, t, dt)
            obj.minX = minX;
            obj.maxX = maxX;
            obj.elemN = elemN;
            obj.nVec = linspace(minX,maxX,(elemN*basis)+1);

            [obj.lBoundary,obj.lBoundaryType,obj.rBoundary,obj.rBoundaryType] = boundaries{:};
            obj.NBc = zeros(obj.globalN,1);
            obj.NBc(1) = obj.lBoundary * (obj.lBoundaryType == "n");
            obj.NBc(end) = obj.rBoundary * (obj.rBoundaryType == "n");

            obj.basisType = basis;
            obj.t = t;
            obj.dt = dt;
        end

        function self = ElemSetup(self)
            for i = 1:self.elemN
                I = self.basisType * (i-1) + 1;
                self.elems(i) = ElemObj(i,self.nVec(I:I+self.basisType),self.basisType);
            end
        end

        function self = SetParams(self,dVec,lVec,fVec)
            self = self.ElemSetup();

            dPointer = 1;
            lPointer = 1;
            fPointer = 1;
            for i = 1:self.elemN
                for j = 1:self.basisType+1
                    while self.elems(i).x(j) > dVec(1,dPointer)
                        dPointer = dPointer + 1;
                    end
                    self.elems(i).D(j) = dVec(2,dPointer);

                    while self.elems(i).x(j) > lVec(1,lPointer)
                        lPointer = lPointer + 1;
                    end
                    self.elems(i).L(j) = lVec(2,lPointer);

                    while self.elems(i).x(j) > fVec(1,fPointer)
                        fPointer = fPointer + 1;
                    end
                    self.elems(i).F(j) = fVec(2,fPointer);
                end
            end
        end

        function self = GlobalSetup(self)
            KG = zeros(self.globalN);
            MG = zeros(self.globalN);
            FG = zeros(self.globalN,1);

            [psi,dpsi] = self.BasisFunc(self.basisType);

            for i = 1:self.elemN
                [K,M,F] = self.elems(i).LocalMatrix(psi,dpsi,4);
                I = self.basisType*(i-1)+1; %adjusting the insertion index for all basis type
                I2 = I+self.basisType;
                % Insert local matrices into global matrices
                KG(I:I2,I:I2) = KG(I:I2,I:I2) + K;
                MG(I:I2,I:I2) = MG(I:I2,I:I2) + M;

                FG(I:I2,1) = FG(I:I2,1) + F;
            end
            self.KGlobal = KG;
            self.MGlobal = MG;
            self.FGlobal = FG;
        end

        function self = Solve(self,theta,previous)
            if ~self.transient
                self.KGlobal = previous.KGlobal;
                self.MGlobal = previous.MGlobal;
                self.FGlobal = previous.FGlobal;
                self.elems = previous.elems;
            else
                self = self.GlobalSetup();
            end
            res = (self.MGlobal-(1-theta).*self.dt.*self.KGlobal)*previous.solution;
            res = res + self.dt*theta*(self.FGlobal + self.NBc);
            res = res + self.dt*(1-theta)*(previous.FGlobal + previous.NBc);

            gm = self.MGlobal+(theta.*self.dt.*self.KGlobal);

            if self.lBoundaryType == "d"
                gm(1,:) = [1,zeros(1,self.globalN-1)];
                res(1)=self.lBoundary;
            end
            if self.rBoundaryType == "d"
                gm(end,:) = [zeros(1,self.globalN-1),1];
                res(end)=self.rBoundary;
            end

            self.solution = gm\res;

        end

        function nodes = get.globalN(obj)
            nodes = (obj.elemN*obj.basisType)+1; % Number of nodes is 1 greater than elements times basis order
        end



    end

    methods (Static)
        function [psi,dpsi] = BasisFunc(basis)
            funcs = {
                {@(x)(1-x)/2, @(x)(1+x)/2} ...
                {@(x)-0.5,@(x)0.5};
                {@(x)0.5.*x.*(x-1), @(x)1-x.^2, @(x)0.5.*x.*(1+x)} ...
                {@(x)x-0.5,@(x)-2*x,@(x)x+0.5}
                };
            try
                psi = funcs{basis};
                dpsi = funcs{basis+2};
            catch
                causeException = MException('MeshObj:BasisFunc',"Invalid basis function order");
                throw(causeException);
            end
        end

        function [p,w] = GQScheme(order)
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
                };

            try
                gq = data{order};
                p= gq(1,:);
                w= gq(2,:);
            catch
                causeException = MException('MeshObj:GQ',"Invalid quadrature order");
                throw(causeException);
            end
        end

    end

end