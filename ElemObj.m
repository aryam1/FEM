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

                psiVals = cellfun(@(fun) fun(xipt), psi);
                dpsiVals = cellfun(@(fun) fun(xipt), dpsi);

                psiMatrix = psiVals.*psiVals';
                dpsiMatrix = dpsiVals.*dpsiVals';
                
                totalD = psiVals*(self.D)';
                totalL = psiVals*(self.L)';
                totalF = psiVals*(self.F)';

                F = F + gw*totalF*self.J;
                M1 = self.J*gw*psiMatrix;
                M = M + M1;
                K = K + (totalD/self.J)*gw*dpsiMatrix - totalL*M1;
            end
        end

        
    end
end