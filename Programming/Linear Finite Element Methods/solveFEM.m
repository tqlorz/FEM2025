function [node, elem, u_h] = solveFEM(h, f)
    %% Generate mesh for the unit square
    [node, elem] = squaremesh([0, 1, 0, 1], h);
    
    %% Assembling the matrix
    [A, area] = assemblingstandard(node, elem);
    % [A, area] = assemblingsparse(node, elem);
    
    %% Right hand side
    b = righthandside(node, elem, area, f);
    
    %% Set boundary conditions
    bdFlag = setboundary(node, elem, 'Dirichlet', 'x == 0 | x == 1 | y == 0 | y == 1');
    [bdNode, bdEdge, isBdNode] = findboundary(elem, bdFlag);
    
    %% Cut off matrix
    A_sub = A(logical(~isBdNode), logical(~isBdNode));
    b_sub = b(logical(~isBdNode));
    
    %% Solve the algebraic equation
    u_sub = A_sub \ b_sub;
    
    u_h = zeros(size(node, 1), 1);
    u_h(~isBdNode) = u_sub;  
    u_h(isBdNode) = 0;     
end