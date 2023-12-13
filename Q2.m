basis = 2;

xElems = 50;
xDist = [5/3000 5/1000 0.01; 5 10 5];
x0 = 0;

t0 = 0;
t1 = 30;
dt = 0.01;

IC = 0;
RB = 0;
RBType = 'd';
LB = 70.245;
LBType = 'd';

DDist = [5/3000 5/1000 0.01; 4e-6 5e-6 2e-6];
LBDist = [5/3000 5/1000 0.01; 0 0.01 0.01]; % Beta
LGDist = [5/3000 5/1000 0.01; 0.02 0.02 0.02]; % Gamma
LDist = [LBDist(1,:);-(LBDist(2,:) + LGDist(2,:))];
FDist = [1;0];

Ksol = []; % Empty array for K values
dif = (max(DDist(2,:))-min(DDist(2,:)))/2; % Difference between max and min D values
adjust = linspace(-dif,dif,11); % Array of D values to adjust by
for i = 1:11
    adjusted=DDist(2,:)+adjust(i); % Adjusted D values
    DDist = [DDist(1,:);adjusted]; % Adjusted DDist
    matrix = OORun(basis,xElems,xDist,x0,t0,t1,dt,IC,LB,RB,LBType,RBType,DDist,LDist,FDist);
    [sol] = [matrix.solution]; % Solution vector
    
    DInd = matrix(1).xInd(0.005); % Index where x = 0.005
    cDose = sol(DInd,:); % Concentration at x = 0.005
    cDoseEff = cDose(cDose > 40); % Concentration at x = 0.005 where cDose > 40
    K = trapz(cDoseEff)*dt; % K value
    Ksol(i) = K; % K values
end 