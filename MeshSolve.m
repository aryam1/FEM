function [solVector] = MeshSolve(mesh,t)
M=1;
K=1;
solVector = (mesh.nvec+t)*t;
end