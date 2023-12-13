t = tiledlayout(2,2);
xlabel(t,'Number of elements')
ylabel(t,"L2 Norm")
nexttile
leg = [];
hold on
for i = 1:4:19
plot(Nx,E1B(i,:));
leg = [leg,Nx(i)];
end
title('Effect of increasing elements on L2 Norm of Linear Backwards Differencing')
legend(cellstr(num2str(leg', '%-d timesteps')));
nexttile

hold on
for i = 1:4:19
plot(Nx,E2B(i,:));
end
legend(cellstr(num2str(leg', '%-d timesteps')));
title('Effect of increasing elements on L2 Norm of Quadratic Backwards Differencing')
nexttile

hold on
for i = 1:4:19
plot(Nx,E1CN(i,:));
end
legend(cellstr(num2str(leg', '%-d timesteps')));
title('Effect of increasing elements on L2 Norm of Linear Crank Nicolson')
nexttile

hold on
for i = 1:4:19
plot(Nx,E2CN(i,:));
end
legend(cellstr(num2str(leg', '%-d timesteps')));
title('Effect of increasing elements on L2 Norm of Quadratic Crank Nicolson')

clear t
figure

t = tiledlayout(2,2);
xlabel(t,'Number of timesteps')
ylabel(t,"L2 Norm")
nexttile
leg = [];
hold on
for i = 1:4:19
plot(Nx,E1B(:,i));
leg = [leg,Nx(i)];
end
title('Effect of increasing timestep on L2 Norm of Linear Backwards Differencing')
legend(cellstr(num2str(leg', '%-d elements')));
nexttile

hold on
for i = 1:4:19
plot(Nx,E2B(:,i));
end
legend(cellstr(num2str(leg', '%-d elements')));
title('Effect of increasing timestep on L2 Norm of Quadratic Backwards Differencing')
nexttile

hold on
for i = 1:4:19
plot(Nx,E1CN(:,i));
end
legend(cellstr(num2str(leg', '%-d elements')));
title('Effect of increasing timestep on L2 Norm of Linear Crank Nicolson')
nexttile

hold on
for i = 1:4:19
plot(Nx,E2CN(:,i));
end
legend(cellstr(num2str(leg', '%-d elements')));
title('Effect of increasing timestep on L2 Norm of Quadratic Crank Nicolson')

clear t
figure

plot(Nx,E2CN(:,14));
%legend(cellstr(num2str(leg', '%-d elements')));
title('Effect of increasing timestep on L2 Norm of Quadratic Crank Nicolson')









