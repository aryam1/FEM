function [ psi ] = BasisFunc(order)
switch order
    case 1
        psi = {@(x)(1-x)/2,@(x)(1+x)/2};
    case 2
        psi = {@(x)x*(x-1)/2,@(x)1-x^2,@(x)x*(1+x)/2};
end
end