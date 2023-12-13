basis = 1;      % 1 = linear, 2 = quadratic

xElems = 10;   % Number of elements in x direction
xDist = [1;1];  % Spatial distribution of elements
xDist = [0.2 0.9 1; 5 2 3];
x0 = 0;         % Left boundary

t0 = 0;         % Initial time
t1 = 1;         % Final time
dt = 0.01;      % Time step
Nt = 1+t1/dt;   % Number of time steps

tvals = t0:dt:t1;  % Time values

matrix = MeshObj.empty(0,Nt); % Empty array of MeshObj objects

IC = 0; % Initial condition
LB = 0; % Left boundary condition
RB = 1; % Right boundary condition

DDist = [1;1]; % Diffusion coefficient distribution
LDist = [1;0]; % Reaction coefficient distribution
FDist = [1;0]; % Source coefficient distribution

for i = 1:t1/dt + 1
    % Create MeshObj object for each time step
    matrix(i) = MeshObj(x0, xElems, xDist, {LB,'d',RB,'d'}, basis, tvals(i), dt);
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
