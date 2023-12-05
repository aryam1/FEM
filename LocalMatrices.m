function [K,M] = LocalMatrices(eID, mesh)

    J = mesh.elem(eID).J;
    L = mesh.elem(eID).L;
    D = mesh.elem(eID).D;
    
    M = J/3 * [2 1; 1 2];
    K= ((D/J)/2 * [1 -1; -1 1]) - (L*M);
end
