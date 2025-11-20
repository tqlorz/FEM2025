clear;
clc;

%% Set parameters 
h = pow2(-2);
f = @(x) 2 * pi^2 .* sin(pi * x(:, 1)) .* sin(pi * x(:, 2));

%% Solve FEM problem
[node, elem, edge, elem2dof, u_h] = solveFEM(h, f);

NV = size(node, 1);
NE = size(edge, 1);
Nu = NV + NE;
allnode = zeros(Nu, 2); 
allnode(1 : NV, :) = node;

allnode(elem2dof(:, 4), :) = 0.5 * (node(elem2dof(:, 2), :) + node(elem2dof(:, 3), :));
allnode(elem2dof(:, 5), :) = 0.5 * (node(elem2dof(:, 1), :) + node(elem2dof(:, 3), :));
allnode(elem2dof(:, 6), :) = 0.5 * (node(elem2dof(:, 1), :) + node(elem2dof(:, 2), :));
