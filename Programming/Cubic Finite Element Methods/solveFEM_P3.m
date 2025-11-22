function [node, elem, edge, elem2dof, u_h] = solveFEM_P3(h, f)
    %% Generate mesh for the unit square
    [node, elem] = squaremesh([0, 1, 0, 1], h);
    [elem2dof, elem2edge, edge, bdDof] = dofP3(elem);

    %% Calculate parameters
    NV = size(node, 1);
    NE = size(edge, 1);
    NT = size(elem, 1);
    Nu = NV + 2 * NE + NT;

    %% Assembling the algebra equation
    [A, b] = assembling_P3(node, elem, edge, elem2dof, f);
    
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