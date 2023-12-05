function [cout] = FESolver(mesh,D,L,f,bc0,bc1,type0,type1)

mat = zeros(mesh.ngn); % Initialize the global stiffness matrix.
res = zeros(mesh.ngn,1); % Initialize the solution vector.

for i = 1:mesh.ne
    % Compute local matrices.
    diffmat = LaplaceElemMatrix(D,i,mesh);
    reactmat = LaplaceReactMatrix(L,i,mesh);
    sourcemat = LaplaceSourceMatrix(f,i,mesh);

    % Insert local matrices into global matrices
    mat(i:i+1,i:i+1) = mat(i:i+1,i:i+1) + diffmat - reactmat;
    res(i:i+1) = res(i:i+1) + sourcemat;
end

% Apply boundary conditions at the left end.
switch type0
    case 'n' %Neumann
        res(1)=res(1)-(bc0*D);    
    case 'd' %Dirichlet
        mat(1,:) = [1,zeros(1,mesh.ngn-1)];
        res(1)=bc0;
    otherwise
        error('FESolver:invalidBoundary','Invalid boundary condition type entered');
end

% Apply boundary conditions at the right end.
switch type1
    case 'n' %Neumann
        res(mesh.ngn)=res(mesh.ngn)+(bc1*D);    
    case 'd' %Dirichlet
        mat(mesh.ngn,:) = [zeros(1,mesh.ngn-1),1];
        res(mesh.ngn)=bc1;
    otherwise
        error('FESolver:invalidBoundary','Invalid boundary condition type entered')
end

% Solve the linear system of equations using backslash operator.
cout = mat\res;
end
