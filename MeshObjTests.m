tol = 1e-10;

xDist = [1;1];
xElems = 10;
x0 = 0;
t0 = 0;
t1 = 1;
dt = 0.01;

basis = 1;

tvals = t0:dt:t1;

IC = 0; 
LB = 0; 
RB = 1;

DDist = [1;1];
LDist = [1;0];
FDist = [1;0];

msh = MeshObj(x0,xElems,xDist,{LB,'d',RB,'d'}, basis, tvals(1), dt);
xvals = msh.nVec;
msh = msh.SetParams(DDist,LDist,FDist); 
msh.solution = IC*ones(length(xvals),1);
msh = msh.GlobalSetup();

%% Test 1: Check symmetry of Mass matrix
M = msh.MGlobal;

assert(abs(M(1, 2) - M(2, 1)) <= tol)

%% Test 2: Check symmetry of Stiffness matrix

K = msh.KGlobal;

assert(abs(K(1, 2) - K(2, 1)) <= tol)


%% Test 3: Check against computed stiffness matrix

K = msh.KGlobal;
error = rmse(K(1:2,1:2),[10 -10; -10 20]);

assert(error(1) <= tol)


%% Test 4: Check against computed mass matrix

M = msh.MGlobal;

error = rmse(M(1:2,1:2),[10/300	5/300; 5/300 20/300]);

assert(error(1) <= tol)

%% Test 5: Check element distribution 

assert(sum(xvals - linspace(0,1,11)) <= tol)


%% Test 6: Check quadratic matrix symmetry
basis = 2;
xElems = 5;
xDist = [0.5 1; 2 3];
msh = MeshObj(x0,xElems,xDist,{LB,'d',RB,'d'}, basis, tvals(1), dt);
xvals = msh.nVec;
msh = msh.SetParams(DDist,LDist,FDist); 
msh.solution = IC*ones(length(xvals),1);
msh = msh.GlobalSetup();

M = msh.MGlobal;

assert(sum(M(1:3,1) - M(1,1:3)') <= tol)


%% Test 7: Check erroneous material parameter 
try
    msh.SetParams(DDist,LDist,"erroneous")
catch e
    assert(e.identifier == "MeshObj:SetParams")
end

%% Test 8: Check erroneous theta scheme parameter 
try
    msh.Solve(5,msh)
catch e
    assert(e.identifier == "MeshObj:Solve")
end

%% Test 9: Check erroneous previous mesh parameter 
try
    msh.Solve(0.5,"Hello World")
catch e
    assert(e.identifier == "MeshObj:Solve")
end

%% Test 10: Try compute L2 Norm for a different equation
try
    DDist = [1;5];
    msh = MeshObj(x0,xElems,xDist,{LB,'d',RB,'d'}, basis, tvals(1), dt);
    xvals = msh.nVec;
    msh = msh.SetParams(DDist,LDist,FDist); 
    msh.solution = IC*ones(length(xvals),1);
    msh = msh.GlobalSetup();
    msh = msh.Solve(0.5,msh);
    msh.L2Norm();
catch e
    assert(e.identifier == "MeshObj:L2Norm")
end

%% Test 11: Try pull cubic basis functions
try
    MeshObj.BasisFunc(3);
catch e
    assert(e.identifier == "MeshObj:BasisFunc")
end

%% Test 11: Try pull 6th order quadrature scheme
try
    MeshObj.GQScheme(6);
catch e
    assert(e.identifier == "MeshObj:GQ")
end