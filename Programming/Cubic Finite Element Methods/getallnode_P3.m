function allnode = getallnode_P3(node, elem, edge, elem2dof)
    NV = size(node, 1);
    NE = size(edge, 1);
    NT = size(elem, 1);
    Nu = NV + 2 * NE + NT;
    
    allnode = zeros(Nu, 2);
    allnode(1 : NV, :) = node;    
    allnode(elem2dof(:, 4), :) = 2 / 3 * node(elem2dof(:, 2), :) +  1 / 3 * node(elem2dof(:, 3), :);
    allnode(elem2dof(:, 5), :) = 1 / 3 * node(elem2dof(:, 2), :) +  2 / 3 * node(elem2dof(:, 3), :);
    allnode(elem2dof(:, 6), :) = 1 / 3 * node(elem2dof(:, 1), :) +  2 / 3 * node(elem2dof(:, 3), :);
    allnode(elem2dof(:, 7), :) = 2 / 3 * node(elem2dof(:, 1), :) +  1 / 3 * node(elem2dof(:, 3), :);
    allnode(elem2dof(:, 8), :) = 2 / 3 * node(elem2dof(:, 1), :) +  1 / 3 * node(elem2dof(:, 2), :);
    allnode(elem2dof(:, 9), :) = 1 / 3 * node(elem2dof(:, 1), :) +  2 / 3 * node(elem2dof(:, 2), :);
    allnode(elem2dof(:, 10), :) = 1 / 3 * node(elem2dof(:, 1), :) + ...
                                  1 / 3 * node(elem2dof(:, 2), :) + ...
                                  1 / 3 * node(elem2dof(:, 3), :);
end