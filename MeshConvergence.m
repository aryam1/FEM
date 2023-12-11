basis = 1;
xDist = [1;1];
x0 = 0;

t0 = 0;
t1 = 1;

IC = 0;
LB = 0;
RB = 1;

DDist = [1;1];
LDist = [1;0];
FDist = [1;0];

Nx = 10:100:10000;
Nt = 10:100:10000;

E = zeros(length(Nt),length(Nx));
iterations = size(E);
totalIter = prod(iterations);

ppm = ParforProgressbar(totalIter);

parfor ix = 1:totalIter
    matrix = MeshObj.empty();
    [row,col] = ind2sub(iterations,ix);
    xelem = Nx(col);
    telem = Nt(row);
    
    tvals = linspace(t0,t1,telem+1)
    dt = t1/telem;

    for i = 1:telem+1
        matrix(i) = MeshObj(x0, xelem, xDist, {LB,'d',RB,'d'}, basis, tvals(i), dt);
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
    L2 = sum([matrix.L2]);
    E(ix) = dt*trapz(L2);
    pause(100/totalIter);
    ppm.increment();
end
E= reshape(E,length(Nx),length(Nt));
delete(ppm);