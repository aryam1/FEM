function Analysis(dx,dt,tmax,theta)

x=linspace(0,1,1+1/dx);
t=linspace(0,tmax,1+tmax/dt);

tran = zeros(1+1/dx,1+tmax/dt);
num = zeros(1+1/dx,1+tmax/dt);


for i = 1:1/dx
    for j = 1:tmax/dt
        tran(i,j) = TransientAnalyticSoln(x(i),t(j));
    end
end

num = Run(dx,dt,tmax,theta);
err = immse(tran,num)
end
