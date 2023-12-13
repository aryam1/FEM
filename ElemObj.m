classdef ElemObj
    properties
        D (1,:) double % Diffusion vector
        L (1,:) double % Reaction vector
        F (1,:) double % Source vector
        x (1,:) double % Node coordinates
        id {mustBeNumeric} = 0 % Element ID
        t {mustBeNumeric} = 0 % Current time for mesh
        basis = 1 % Basis function order
        J % Jacobian
    end
    properties (Dependent)
    end
    methods
        function obj = ElemObj(id,x,basis)
            obj.id = id;
            obj.x = x;
            obj.basis = basis;
            obj.J = (x(end)-x(1))/2; % Jacobian calculation
        end

        function [K,M,F] = LocalMatrix(self,psi,dpsi,p,w) % Local matrix calculation
            order = self.basis+1; 
            K = zeros(order); 
            M = zeros(order);
            F = zeros(order,1);

            for i = 1:length(p)
                xipt = p(i); % Gauss point
                gw = w(i); % Gauss weight

                % Evaluate basis functions for each Gauss point
                psiVals = cellfun(@(fun) fun(xipt), psi); 
                % Evaluate basis function derivatives for each Gauss point
                dpsiVals = cellfun(@(fun) fun(xipt), dpsi); 

                psiMatrix = psiVals.*psiVals'; % Basis function product matrix
                dpsiMatrix = dpsiVals.*dpsiVals'; % Basis function derivative product matrix
                
                totalD = psiVals*(self.D)'; % Varying Diffusion coefficient interpolation
                totalL = psiVals*(self.L)'; % Varying Reaction coefficient interpolation
                totalF = psiVals*(self.F)'; % Varying Source coefficient interpolation

                F = F + gw*totalF*self.J; % Source vector
                M1 = self.J*gw*psiMatrix; % Mass matrix
                M = M + M1; 
                K = K + (totalD/self.J)*gw*dpsiMatrix - totalL*M1; % Stiffness matrix
            end
        end

        
    end
end