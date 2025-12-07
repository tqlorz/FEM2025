clc;
clear;

%% Set parameters 
h = pow2(-5);
f = @(x) 24*(x(:,1).^2 .* (1-x(:,1)).^2 + x(:,2).^2 .* (1-x(:,2)).^2) ...
       + 8*(1 - 6*x(:,1) + 6*x(:,1).^2) .* (1 - 6*x(:,2) + 6*x(:,2).^2);

%% Solve FEM problem
[node, elem, edge, elem2dof, u_h] = test_solveFEM(h, f);
allnode = getallnode(node, edge, elem2dof);

%% Draw picture of mesh
if h >= pow2(-3)
    figure(1);
    hold on;
    axis on;
    showmesh(node, elem);
    findelem(node, elem);       % plot indices of all triangles
    findnode(node);             % plot indices of all vertices
    findedgedof(node, edge);    % plot indices of all edge
end

%% Draw picture of FEM
plotsolution(node, elem, elem2dof, u_h)

function [node, elem, edge, elem2dof, u_h] = test_solveFEM(h, f)
    %% Generate mesh for the unit square
    [node, elem] = squaremesh([0, 1, 0, 1], h);
    [~, edge, bdDof] = dofP2(elem);

    %% Calculate parameters
    NV = size(node, 1);
    NE = size(edge, 1);
    Nu = NV + NE;

    %% Generate auxiliary datas
    T = auxstructure(elem);
    elem2dof = [elem, NV + T.elem2edge];

    %% Assembling the algebra equation
    [A, b] = test_assembling(node, elem, edge, elem2dof, f);
    
    %% Set boundary conditions
    isBdDof = false(Nu, 1);
    isBdDof(bdDof(:)) = 1;

    %% Cut off matrix
    A_sub = A(logical(~isBdDof), logical(~isBdDof));
    b_sub = b(logical(~isBdDof));
    
    %% Solve the algebraic equation
    u_sub = A_sub \ b_sub;
    
    u_h = zeros(Nu, 1);
    u_h(~isBdDof) = u_sub;  
    u_h(isBdDof) = 0;     
end

function [A, b] = test_assembling(node, elem, edge, elem2dof, f)
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

