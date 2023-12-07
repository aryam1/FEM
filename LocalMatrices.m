function [K,M] = LocalMatrices(eID, mesh,basis)

    J = mesh.elem(eID).J;
    L = mesh.elem(eID).L;
    D = mesh.elem(eID).D;
    [p,w]= GQScheme(2);
    psi = BasisFunc(basis);

    switch basis
        case 1
            psi0 = psi{1}(p);
            psi1 = psi{2}(p);
            M = J * [w*(psi0.*psi0)' w*(psi1.*psi0)';
                     w*(psi0.*psi1)' w*(psi1.*psi1)'];
            
            K= ((D/J)/2 * [1 -1; -1 1]) - (L*M);
        
        case 2
            psi0 = psi{1}(p);
            psi1 = psi{2}(p);
            psi2 = psi{3}(p);
            M = J * [w*(psi0.*psi0)' w*(psi1.*psi0)' w*(psi2.*psi0)';
                     w*(psi0.*psi1)' w*(psi1.*psi1)' w*(psi2.*psi1)';
                     w*(psi0.*psi2)' w*(psi1.*psi2)' w*(psi2.*psi2)'];
            
            K= ((D/J)/2 * [1 -1; -1 1]) - (L*M);
    end

    
    %M = J/4 * [w*(p.^2-p.*2+1)', w*(1-p.^2)';w*(1-p.^2)', w*(1+2.*p+p.^2)'];

end
