function [K,M] = LocalMatrices2(eID, mesh,basis)

    order = basis+1;    

    J = mesh.elem(eID).J;
    L = mesh.elem(eID).L;
    D = mesh.elem(eID).D;

    [p,w]= GQScheme(4);
    [psi,dpsi] = BasisFunc(basis);
    K = zeros(order);
    M = zeros(order);

    for i = 1:4
        xipt = p(i);
        gw = w(i);

        for m = 1:order
            for n = 1:order
                M1 = J*gw*(psi{n}(xipt)*psi{m}(xipt));
                M(n,m) = M(n,m) + M1;
                K(n,m) = K(n,m) + ((D/J)*gw*dpsi{n}(xipt)*dpsi{m}(xipt))-L*M1;
            end
        end
    end

end
