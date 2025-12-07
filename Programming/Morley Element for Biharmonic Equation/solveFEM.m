function [node, elem, edge, elem2dof, u_h] = solveFEM(h, f)
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
    [A, b] = assembling(node, elem, edge, elem2dof, f);
    
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