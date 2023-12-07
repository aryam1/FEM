function [psi,dpsi] = BasisFunc(order)
switch order
    case 1
        psi = {@(x)(1-x)/2, @(x)(1+x)/2};
        dpsi = {@(x)-0.5,@(x)0.5};
    case 2
        psi = {@(x)x.*(x-1)/2,@(x)1-x.^2,@(x)x.*(1+x)/2};
        dpsi = {@(x)x-0.5,@(x)-2.*x,@(x)x+0.5};
end
end