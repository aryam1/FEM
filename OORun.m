function [matrix] = OORun(basis,xElems,xDist,x0,t0,t1,dt,IC,LB,RB,LBType,RBType,DDist,LDist,FDist)

Nt = 1+t1/dt;   % Number of time steps

tvals = t0:dt:t1;  % Time values

matrix = MeshObj.empty(0,Nt); % Empty array of MeshObj objects

for i = 1:t1/dt + 1
    % Create MeshObj object for each time step
    matrix(i) = MeshObj(x0, xElems, xDist, {LB,LBType,RB,RBType}, basis, tvals(i), dt);
    % Set parameters for first MeshObj object
    if i == 1
        xvals = matrix(i).nVec; % x values for first MeshObj object
        matrix(i) = matrix(i).SetParams(DDist,LDist,FDist); % Set parameters 
        matrix(i).solution = IC*ones(length(xvals),1); % Set initial condition
        matrix(i) = matrix(i).GlobalSetup(); % Set up global matrix
        continue % Skip to next iteration
    end
    % Solve and calulate L2 norm for all other MeshObj objects
    matrix(i) = matrix(i).Solve(0.5,matrix(i-1));
    matrix(i) = matrix(i).L2Norm();
end
[sol] = [matrix.solution]; % Solution vector
end