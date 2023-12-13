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

Nx = [[10:10:100],[200:100:1000]]; 
Nt = [[10:10:100],[200:100:1000]];

for theta = [0 0.5 1] 
for basis = [1,2] 
    E = zeros(length(Nt),length(Nx)); 
    iterations = size(E); 
    totalIter = prod(iterations); 
    
    parfor ix = 1:totalIter %parallelized for loop
        matrix = MeshObj.empty(); %creates empty array of MeshObj
        [row,col] = ind2sub(iterations,ix); %converts linear index to subscript index
        xelem = Nx(col); %number of elements in x direction
        telem = Nt(row); %number of elements in t direction
        
        tvals = linspace(t0,t1,telem+1) %creates array of t values
        dt = t1/telem;
        % Solves system
        for i = 1:telem+1 
            matrix(i) = MeshObj(x0, xelem, xDist, {LB,'d',RB,'d'}, basis, tvals(i), dt);
            if i == 1 
                xvals = matrix(i).nVec;
                matrix(i) = matrix(i).SetParams(DDist,LDist,FDist);
                matrix(i).solution = IC*ones(length(xvals),1);
                matrix(i) = matrix(i).GlobalSetup();
                continue
            end
            matrix(i) = matrix(i).Solve(theta,matrix(i-1)); 
            matrix(i) = matrix(i).L2Norm();
        end
        L2 = sqrt(sum([matrix.L2])); % L2 norm of solution
        matrix = MeshObj.empty();
        E(ix) = dt*trapz(L2); % Integrated L2 norm of solution
    end
    E= reshape(E,length(Nx),length(Nt)); %reshapes E into a matrix
    name = "%dbasis_%d_theta_E.mat"; %saves E as a .mat file
    save(sprintf(name,basis,theta),'E','Nx','Nt') 
    clear E;
end
end
