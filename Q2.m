basis = 2;

xElems = 100;
xDist = [0.01;1];
%xDist = [5/3000 5/1000 0.01; 2 8 5];
x0 = 0;

t0 = 0;
t1 = 30;
dt = 0.01;
Nt = 1+t1/dt;

tvals = t0:dt:t1;

IC = 0;
RB = 0;
LB = 70.01856;

DDist = [5/3000 5/1000 0.01; 4e-6 5e-6 2e-6];
LBDist = [5/3000 5/1000 0.01; 0 0.01 0.01]; % Beta
LGDist = [5/3000 5/1000 0.01; 0.02 0.02 0.02]; % Gamma
LDist = [LBDist(1,:);-(LBDist(2,:) + LGDist(2,:))];
FDist = [1;0];

matrix = MeshObj.empty(0,Nt);

for i = 1:t1/dt + 1
    if mod(i,100) == 0
        disp(i);
    end
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

DInd = matrix(1).xInd(0.005);
cDose = sol(DInd,:);
cDoseEff = cDose(cDose > 40);
K = trapz(cDoseEff)*dt;
