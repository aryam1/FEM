function [sol] = Run(dx,dt,tmax,theta,basis)

Nt = tmax/dt;
Nx = 1/dx;

tVal= 0:dt:tmax;

mesh = OneDimLinearMeshGen(0,1,Nx);

for t = 1:Nt+1
    mesh.t = tVal(t);
    mesh.sol(1)=0;
    mesh.sol(end)=1;
    meshMatrix(t) = mesh;
end

sol = zeros(basis*mesh.ne+1,Nt+1);

sol(end,:)=1; % all time at end space equals 1 (last row)
sol(1,:)=0; % all time at 0 space equals 0 (first row)
sol(:,1)=0; % all space at 0 time equals 0 (first column)
meshMatrix(1).sol = sol(:,1);

for i = 2:Nt+1
    newSol = MeshSolve(meshMatrix(i),meshMatrix(i-1),0,1,dt,theta,basis);
    meshMatrix(i).sol = newSol;
    sol(:,i) = newSol;
end

end
