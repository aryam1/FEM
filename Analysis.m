function Analysis(dx,dt,tmax,theta,basis)

x=linspace(0,1,1+1/dx);
t=linspace(0,tmax,1+tmax/dt);

tran = zeros(1+1/dx,1+tmax/dt);
num = Run(dx,dt,tmax,theta,basis);
E = zeros(1,1+tmax/dt);
[p,w] = GQScheme(3);

for j = 1:tmax/dt+1
    for i = 1:1/dx
        c0 = num(i,j); c1 = num(i+1,j);
        x0 = x(i); x1 = x(i+1);
        J=1;
        tran(i,j) = TransientAnalyticSoln(x(i),t(j));
        E1 = J*w.*(TransientAnalyticSoln(0.5*x0*(1-p)+0.5*x1*(1+p),t(j)) - 0.5*c0*(1-p) - 0.5*c1*(1+p)).^2;
        E(j)=E(j)+sum(E1);
    end
end
plot(E)
mean(E)
err = immse(tran,num);
end
