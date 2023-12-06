function Analysis(dx,dt,tmax,theta)

x=linspace(0,1,1+1/dx);
t=linspace(0,tmax,1+tmax/dt);

tran = zeros(1+1/dx,1+tmax/dt);
num = Run(dx,dt,tmax,theta);

for j = 1:tmax/dt+1
    for i = 1:1/dx+1
        tran(i,j) = TransientAnalyticSoln(x(i),t(j));
        E(i) = tran(i,j) - num(i,j);
    end
end

err = immse(tran,num)
end
