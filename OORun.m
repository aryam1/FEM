basis = 1;

x0 = 0;
x1 = 1;
dx = 0.01;
Nx = 1+(basis * x1/dx);

t0 = 0;
t1 = 2;
dt = 0.01;
Nt = 1+t1/dt;

xvals = x0:dx:x1;
tvals = t0:dt:t1;

%matrix = MeshObj.empty(0,Nt);

IC = 0;
LB = 0;
RB = 1;

sol = zeros(Nx,Nt);
sol(:,1) = IC;


% for i = 1:t1/dt + 1
%     matrix(i) = MeshObj(x0,x1,x1/dx,{LB,'d',RB,'d'},basis,tvals(i),dt);
%     if i == 1
%         matrix(i) = matrix(i).SetParams([1;1],[1;0],[1;0]);
%         matrix(i).solution = sol(:,1);
%         matrix(i) = matrix(i).GlobalSetup();
%         continue
%     end
%     matrix(i) = matrix(i).Solve(0.5,matrix(i-1));
% end

previous = MeshObj.empty;
for i = 1:t1/dt + 1
    mesh = MeshObj(x0,x1,x1/dx,{LB,'d',RB,'d'},basis,tvals(i),dt);
    if i == 1
        mesh = mesh.SetParams([1;1],[1;0],[1;0]);
        mesh.solution = sol(:,1);
        mesh = mesh.GlobalSetup();
        previous = mesh;
        continue
    end
    mesh = mesh.Solve(0.5,previous);
    previous = mesh;
    sol(:,i) = mesh.solution;
end