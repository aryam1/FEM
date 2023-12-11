basis = 1;

xElems = 100;
xDist = [1;1];
x0 = 0;

t0 = 0;
t1 = 30;
dt = 0.0001;
Nt = 1+t1/dt;

tvals = t0:dt:t1;

matrix = MeshObj.empty(0,Nt);

IC = 30;
LB = 0;
RB = 0;

DDist = [5/3000 5/1000 0.01; 4e-6 5e-6 2e-6];
LBDist = [5/3000 5/1000 0.01; 0 0.01 0.01]; % Beta
LGDist = [5/3000 5/1000 0.01; 0.02 0.02 0.02]; % Gamma
LDist = -(LBDist + LGDist);
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
end
[sol] = [matrix.solution];

Dindex = 0;