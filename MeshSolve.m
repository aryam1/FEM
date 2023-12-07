function [solVector] = MeshSolve(mesh,prev,dt,theta,basis)

Ne = basis*mesh.ne+1;
Kglobal = zeros(Ne); % Initialize the global stiffness matrix.
Mglobal = zeros(Ne); % Initialize the global mass matrix.

for i = 1:mesh.ne
    % Compute local matrices.
    [K,M] = LocalMatrices(i,mesh,basis);
    I = basis*(i-1)+1; %adjusting the insertion index for all basis type

    % Insert local matrices into global matrices
    Kglobal(I:I+basis,I:I+basis) = Kglobal(I:I+basis,I:I+basis) + K;
    Mglobal(I:I+basis,I:I+basis) = Mglobal(I:I+basis,I:I+basis) + M;
end

res = (Mglobal-(1-theta).*dt.*Kglobal)*prev.sol;
%res = res + dt.*theta.*(mesh.F + mesh.NBc);
%res = res + dt.*(1-theta).*(prev.F + prev.NBc); doesnt scale with quad
%basis

gm = Mglobal+(theta.*dt.*Kglobal);

gm(1,:) = [1,zeros(1,Ne-1)];
res(1)=prev.sol(1);
gm(Ne,:) = [zeros(1,Ne-1),1];
res(Ne)=prev.sol(Ne);

% Solve the linear system of equations using backslash operator.
solVector = gm\res;

end