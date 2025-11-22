function allnode = getallnode(node, edge, elem2dof)
    NV = size(node, 1);
    NE = size(edge, 1);
    Nu = NV + NE;
    allnode = zeros(Nu, 2); 
    allnode(1 : NV, :) = node;
    allnode(elem2dof(:, 4), :) = 0.5 * (node(elem2dof(:, 2), :) + node(elem2dof(:, 3), :));
    allnode(elem2dof(:, 5), :) = 0.5 * (node(elem2dof(:, 1), :) + node(elem2dof(:, 3), :));
    allnode(elem2dof(:, 6), :) = 0.5 * (node(elem2dof(:, 1), :) + node(elem2dof(:, 2), :));
end