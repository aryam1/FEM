basis = 2;

xElems = 1000;
xDist = [1;1];
x0 = 0;

t0 = 0;
t1 = 1;
dt = 0.01;
Nt = 1+t1/dt;

tvals = t0:dt:t1;

matrix = MeshObj.empty(0,Nt);

IC = 0;
LB = 0;
RB = 1;

DDist = [1;1];
LDist = [1;0];
FDist = [1;0];

for i = 1:t1/dt + 1
    matrix(i) = MeshObj(x0, xElems, xDist, {LB,'d',RB,'d'}, basis, tvals(i), dt);
    if i == 1
        xvals = matrix(i).nVec;
        matrix(i) = matrix(i).SetParams(DDist,LDist,FDist);
        matrix(i).solution = IC*ones(length(xvals),1);
        matrix(i) = matrix(i).GlobalSetup();
        continue
    end
    matrix(i) = matrix(i).Solve(0.5,matrix(i-1));
    matrix(i) = matrix(i).L2Norm();
end
[sol] = [matrix.solution];
