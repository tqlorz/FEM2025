function [A, b] = assembling(node, elem, edge, elem2dof, f)
    %% Calculate parameters 
    NV = size(node, 1);
    NE = size(edge, 1);
    NT = size(elem, 1);
    Nu = NV + NE;

    %% Assembling the algebra equation
    A = sparse(Nu, Nu);
    b = zeros(Nu, 1);
    [D, area] = gradbasis(node, elem);
    fmid = f((node(elem(:, 1), :) ...
                + node(elem(:, 2), :) ...
                + node(elem(:, 3), :)) / 3);
    [lambda, weight] = quadpts(2);
    coefs(:, :, :) = getcoefs(D);
    for nt = 1 : NT 
        for p = 1 : length(weight)
            
            for nv = 1 : 6
                Hphi(nt, :, :, nv) = 2 * coefs(nt, 1, nv) * D(nt, :, 1).' * D(nt, :, 1) ...
                                   + 2 * coefs(nt, 2, nv) * D(nt, :, 2).' * D(nt, :, 2) ...
                                   + 2 * coefs(nt, 3, nv) * D(nt, :, 3).' * D(nt, :, 3) ...
                                   + coefs(nt, 4, nv) * (D(nt, :, 1).' * D(nt, :, 2) + D(nt, :, 2).' * D(nt, :, 1)) ...
                                   + coefs(nt, 5, nv) * (D(nt, :, 1).' * D(nt, :, 3) + D(nt, :, 3).' * D(nt, :, 1)) ...
                                   + coefs(nt, 6, nv) * (D(nt, :, 2).' * D(nt, :, 3) + D(nt, :, 3).' * D(nt, :, 2));
            end
            for i = 1 : 6 
                for j = 1 : 6
                    A_ij = weight(p) * area(nt) ...
                            .* sum(Hphi(nt, :, :, i) .* Hphi(nt, :, :, j), 'all');
                    A = A + sparse(elem2dof(nt, i), elem2dof(nt, j), A_ij, Nu, Nu);
                end
            end 
            
            for nv = 1 : 6
                phi(nt, nv) = coefs(nt, 1, nv) * lambda(p,1) * lambda(p,1) ...
                            + coefs(nt, 2, nv) * lambda(p,2) * lambda(p,2) ...
                            + coefs(nt, 3, nv) * lambda(p,3) * lambda(p,3) ...
                            + coefs(nt, 4, nv) * lambda(p,1) * lambda(p,2) ...
                            + coefs(nt, 5, nv) * lambda(p,1) * lambda(p,3) ...
                            + coefs(nt, 6, nv) * lambda(p,2) * lambda(p,3);
            end
            for k = 1 : 6
                b_k = weight(p) .* area(nt) .* fmid(nt) .* phi(nt, k);
                b(elem2dof(nt, k)) = b(elem2dof(nt, k)) + b_k; 
            end
        end
    end
end