function [A, b] = assembling_P3(node, elem, edge, elem2dof, f)
    %% Calculate parameters 
    NV = size(node, 1);
    NE = size(edge, 1);
    NT = size(elem, 1);
    Nu = NV + 2 * NE + NT;

    %% Assembling the algebra equation
    A = sparse(Nu, Nu);
    bt = zeros(NT, 10);
    [D, area] = gradbasis(node, elem);
    fmid = f((node(elem(:, 1), :) ...
                + node(elem(:, 2), :) ...
                + node(elem(:, 3), :)) / 3);
    [lambda, weight] = quadpts(3);
    for p = 1 : length(weight)
        dphi(:, :, 1) = (13.5 * lambda(p, 1) * lambda(p, 1) - 9 * lambda(p, 1) + 1) * D(:, :, 1);
        dphi(:, :, 2) = (13.5 * lambda(p, 2) * lambda(p, 2) - 9 * lambda(p, 2) + 1) * D(:, :, 2);
        dphi(:, :, 3) = (13.5 * lambda(p, 3) * lambda(p, 3) - 9 * lambda(p, 3) + 1) * D(:, :, 3);
        dphi(:, :, 4) = 4.5 * ((3 * lambda(p, 2) * lambda(p, 2) - lambda(p, 2)) * D(:, :, 3) ...
                             + lambda(p, 3) * (6 * lambda(p, 2) - 1) * D(:, :, 2));
        dphi(:, :, 5) = 4.5 * ((3 * lambda(p, 3) * lambda(p, 3) - lambda(p, 3)) * D(:, :, 2) ...
                             + lambda(p, 2) * (6 * lambda(p, 3) - 1) * D(:, :, 3));
        dphi(:, :, 6) = 4.5 * ((3 * lambda(p, 3) * lambda(p, 3) - lambda(p, 3)) * D(:, :, 1) ...
                             + lambda(p, 1) * (6 * lambda(p, 3) - 1) * D(:, :, 3));
        dphi(:, :, 7) = 4.5 * ((3 * lambda(p, 1) * lambda(p, 1) - lambda(p, 1)) * D(:, :, 3) ...
                             + lambda(p, 3) * (6 * lambda(p, 1) - 1) * D(:, :, 1));
        dphi(:, :, 8) = 4.5 * ((3 * lambda(p, 1) * lambda(p, 1) - lambda(p, 1)) * D(:, :, 2) ...
                             + lambda(p, 2) * (6 * lambda(p, 1) - 1) * D(:, :, 1));
        dphi(:, :, 9) = 4.5 * ((3 * lambda(p, 2) * lambda(p, 2) - lambda(p, 2)) * D(:, :, 1) ...
                             + lambda(p, 1) * (6 * lambda(p, 2) - 1) * D(:, :, 2));
        dphi(:, :, 10) = 27 * (lambda(p, 1) * lambda(p, 2) * D(:, :, 3) ...
                                + lambda(p, 1) * lambda(p, 3) * D(:, :, 2)  ...
                                + lambda(p, 3) * lambda(p, 2) * D(:, :, 1));
        for i = 1 : 10 
            for j = 1 : 10
                A_ij = weight(p) * area ...
                        .* dot(dphi(:, :, i), dphi(:, :, j), 2);
                A = A + sparse(elem2dof(:, i), elem2dof(:, j), A_ij, Nu, Nu);
            end
        end 
        phi(:, 1) = 0.5 * (3 * lambda(p, 1) - 1) * (3 * lambda(p, 1) - 2) * lambda(p, 1);
        phi(:, 2) = 0.5 * (3 * lambda(p, 2) - 1) * (3 * lambda(p, 2) - 2) * lambda(p, 2);
        phi(:, 3) = 0.5 * (3 * lambda(p, 3) - 1) * (3 * lambda(p, 3) - 2) * lambda(p, 3);
        phi(:, 4) = 4.5 * lambda(p, 3) * lambda(p, 2) * (3 * lambda(p, 2) - 1);
        phi(:, 5) = 4.5 * lambda(p, 3) * lambda(p, 2) * (3 * lambda(p, 3) - 1);
        phi(:, 6) = 4.5 * lambda(p, 1) * lambda(p, 3) * (3 * lambda(p, 3) - 1);
        phi(:, 7) = 4.5 * lambda(p, 1) * lambda(p, 3) * (3 * lambda(p, 1) - 1);
        phi(:, 8) = 4.5 * lambda(p, 1) * lambda(p, 2) * (3 * lambda(p, 1) - 1);
        phi(:, 9) = 4.5 * lambda(p, 1) * lambda(p, 2) * (3 * lambda(p, 2) - 1);
        phi(:, 10) = 27 * lambda(p, 1) * lambda(p, 2) * lambda(p, 3);
        for k = 1 : 10
            bt(:, k) = bt(:, k) + weight(p) .* area .* fmid .* phi(:, k);
        end
    end
    b = accumarray(elem2dof(:), bt(:), [Nu, 1]);
end