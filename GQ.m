function [int] = GQ(func,order)
data = {
    [0; 2]
    [-1/sqrt(3) 1/sqrt(3); 1 1] 
    [-sqrt(3/5) 0 sqrt(3/5); 5/9 8/9 5/9]
    [-sqrt(3/7 + (2/7*sqrt(1.2))) -sqrt(3/7 - (2/7*sqrt(1.2))) sqrt(3/7 - (2/7*sqrt(1.2))) sqrt(3/7 + (2/7*sqrt(1.2))); 
    (18-sqrt(30))/36 (18+sqrt(30))/36 (18+sqrt(30))/36 (18-sqrt(30))/36]
    };

gauss= data{order};
int = gauss(2,:)*func(gauss(1,:))';
end
