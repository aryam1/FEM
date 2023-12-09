function [K,M] = LocalMatrices(eID, mesh,basis)

    J = mesh.elem(eID).J;
    L = mesh.elem(eID).L;
    D = mesh.elem(eID).D;
    [p,w]= GQScheme(4);
    [psi,dpsi] = BasisFunc(basis);
   
    switch basis
        case 1
            psi0 = psi{1}(p); psi1 = psi{2}(p);
            dpsi0 = dpsi{1}(p); dpsi1 = dpsi{2}(p);

            M = J * [w*(psi0.*psi0)' w*(psi1.*psi0)';
                     w*(psi0.*psi1)' w*(psi1.*psi1)'];
            K1 = (D/J) * [dpsi0*dpsi0 dpsi1*dpsi0;
                          dpsi0*dpsi1 dpsi1*dpsi1];
            % Need to multiply by 2 becuase function handle is only
            % returning single value not 1x2 matrix for GQ, so Matlab cant
            % sum it

            K = 2*K1 - (L*M);
        
        case 2
            psi0 = psi{1}(p); psi1 = psi{2}(p); psi2 = psi{3}(p);
            dpsi0 = dpsi{1}(p); dpsi1 = dpsi{2}(p); dpsi2 = dpsi{3}(p);
            
            M = J * [w*(psi0.*psi0)' w*(psi1.*psi0)' w*(psi2.*psi0)';
                     w*(psi0.*psi1)' w*(psi1.*psi1)' w*(psi2.*psi1)';
                     w*(psi0.*psi2)' w*(psi1.*psi2)' w*(psi2.*psi2)'];
            
            K1 = (D/J) * [w*(dpsi0.*dpsi0)' w*(dpsi1.*dpsi0)' w*(dpsi2.*dpsi0)';
                        w*(dpsi0.*dpsi1)' w*(dpsi1.*dpsi1)' w*(dpsi2.*dpsi1)';
                        w*(dpsi0.*dpsi2)' w*(dpsi1.*dpsi2)' w*(dpsi2.*dpsi2)'];

            K = K1 - (L*M);
    end
end
