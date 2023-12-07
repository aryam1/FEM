function [sol] = Run(dx,dt,tmax,theta,basis)

Nt = tmax/dt;
Nx = 1/dx;

tVal= 0:dt:tmax;

mesh = OneDimLinearMeshGen(0,1,Nx);

for t = 1:Nt+1
    mesh.t = tVal(t);
    meshMatrix(t) = mesh;
end

sol = zeros(basis*mesh.ne+1,Nt+1);

sol(:,1)=0; % all space at 0 time equals 0 (first column)
sol(1,:)=0; % all time at 0 space equals 0 (first row)
sol(basis*mesh.ne+1,:)=1; % all time at end space equals 1 (last row)

for i = 2:Nt+1
    sol(:,i) = MeshSolve(meshMatrix(i),sol(:,i-1),dt,theta,basis);
end

end
