function [sol] = Run()

Nx = 10;
Nt = 10;

xmax = 1;
tmax = 10;

D = 1; L = 1; F = 1;

theta = 0;

tVal= 0:tmax/Nt:tmax;

mesh = OneDimLinearMeshGen(0,xmax,Nx);

for t = 1:Nt+1
    mesh.t = tVal(t);
    mesh.D = D;
    mesh.L = L;
    mesh.F = F;
    meshMatrix(t) = mesh;
end

sol = zeros(Nx+1,Nt+1);

sol(:,1)=0; % all space at 0 time equals 0 (first column)
sol(1,:)=0; % all time at 0 space equals 0 (first row)
sol(Nx+1,:)=1; % all time at end space equals 1 (last row)

for i = 2:Nt+1
    M=1;
    J=1;
    sol(:,i) = MeshSolve(meshMatrix(i),tVal(i));
end

end
