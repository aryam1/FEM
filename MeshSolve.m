function [solVector] = MeshSolve(mesh,cur,dt,theta)

gm = zeros(mesh.ngn);
Kglobal = zeros(mesh.ngn); % Initialize the global stiffness matrix.
Mglobal = zeros(mesh.ngn); % Initialize the global mass matrix.
res = zeros(mesh.ngn); % Initialize the solution vector.

for i = 1:mesh.ne
    % Compute local matrices.
    [K,M] = LocalMatrices(i,mesh);

    % Insert local matrices into global matrices
    Kglobal(i:i+1,i:i+1) = Kglobal(i:i+1,i:i+1) + K;
    Mglobal(i:i+1,i:i+1) = Mglobal(i:i+1,i:i+1) + M;
end

res = (Mglobal-(1-theta).*dt.*Kglobal)*cur;
gm = Mglobal+(theta.*dt.*Kglobal);

gm(1,:) = [1,zeros(1,mesh.ngn-1)];
res(1)=0;
gm(mesh.ngn,:) = [zeros(1,mesh.ngn-1),1];
res(mesh.ngn)=1;


% Solve the linear system of equations using backslash operator.
solVector = gm\res;

end