clc;
clear;

h = power(2, -1);
[node, elem] = squaremesh([0, 1, 0, 1], h);
[~, edge, bdDof] = dofP2(elem);

[node_refined, elem_refined] = uniformrefine(node, elem);

figure(1);
hold on;
axis on;
showmesh(node, elem);
findelem(node, elem);       % plot indices of all triangles
findnode(node);             % plot indices of all vertices
findedgedof(node, edge);    % plot indices of all edge

figure(2);
hold on;
axis on;
showmesh(node_refined, elem_refined);
findelem(node_refined, elem_refined);       % plot indices of all triangles
findnode(node_refined);             % plot indices of all vertices
% findedgedof(node, edge);    % plot indices of all edge