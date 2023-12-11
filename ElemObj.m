classdef ElemObj
    properties
        D (1,:) double
        L (1,:) double
        F (1,:) double
        x (1,:) double
        id {mustBeNumeric} = 0
        t {mustBeNumeric} = 0 % Current time for mesh
        basis = 1
        J % Jacobian
    end
    properties (Dependent)
    end
    methods
        function obj = ElemObj(id,x,basis)
            obj.id = id;
            obj.x = x;
            obj.basis = basis;
            obj.J = (x(end)-x(1))/2;
        end

        function [K,M,F] = LocalMatrix(self,psi,dpsi,p,w)
            order = self.basis+1;
            K = zeros(order);
            M = zeros(order);
            F = zeros(order,1);

            for i = 1:length(p)
                xipt = p(i);
                gw = w(i);
                for j = 1:order
                    d(j)=self.D(j)*gw*(psi{j}(xipt));
                    l(j)=self.D(j)*gw*(psi{j}(xipt));
                    f(j)=self.F(j)*gw*(psi{j}(xipt));
                end
                totalD = sum(d);
                totalL = sum(l);
                totalF = sum(f);
                for m = 1:order
                    F(m,1) = totalF*self.J;
                    for n = 1:order
                        M1 = self.J*gw*(psi{n}(xipt)*psi{m}(xipt));
                        M(n,m) = M(n,m) + M1;
                        K(n,m) = K(n,m) + ((totalD/self.J)*gw*dpsi{n}(xipt)*dpsi{m}(xipt))-totalL*M1;
                    end
                end
            end
        end

        
    end
end