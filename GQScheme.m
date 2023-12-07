function [p,w] = GQScheme(order)
data = {
    [0; 2]
    [-1/sqrt(3) 1/sqrt(3); 1 1]
    [-sqrt(3/5) 0 sqrt(3/5); 5/9 8/9 5/9]
    [-sqrt(3/7 + (2/7*sqrt(1.2))) -sqrt(3/7 - (2/7*sqrt(1.2))) sqrt(3/7 - (2/7*sqrt(1.2))) sqrt(3/7 + (2/7*sqrt(1.2))); 
    (18-sqrt(30))/36 (18+sqrt(30))/36 (18+sqrt(30))/36 (18-sqrt(30))/36]
    };

try
    gq = data{order};
    p= gq(1,:);
    w= gq(2,:);
catch
    causeException = MException('MATLAB:GQ',"Invalid quadrature order");
    throw(causeException);
end

%int = gauss(2,:)*integrand(gauss(1,:))';
end
