clc;
clear;

%% Set parameters 
h = pow2(-2);
f = @(x) 2 * pi^2 .* sin(pi * x(:, 1)) .* sin(pi * x(:, 2));

%% Solve FEM problem
[node, elem, edge, elem2dof, u_h] = solveFEM_P3(h, f);
allnode = getallnode_P3(node, elem, edge, elem2dof);

%% Draw picture of mesh
figure(1);
hold on;
axis on;
showmesh(node, elem);
findelem(node, elem);       % plot indices of all triangles
findnode(allnode);             % plot indices of all vertices

%% Draw picture of FEM
[Xi, Yi] = meshgrid(linspace(0, 1, 100), ...
                    linspace(0, 1, 100));

Zi = griddata(allnode(:,1), allnode(:,2), u_h, Xi, Yi, 'cubic');

figure(2);
contourf(Xi, Yi, Zi, 20);
colorbar;
title('Interpolated Solution');
xlabel('X'); ylabel('Y');

figure(3);
surf(Xi, Yi, Zi, 'EdgeColor', 'none');
colorbar;
title('3D Surface Plot of u_h');
xlabel('X'); ylabel('Y'); zlabel('u_h');

%% Verify convergence in H1-norm
[node, elem] = squaremesh([0, 1, 0, 1], 0.25);
for k = 1 : 4
    exactu = inline('sin(pi*pxy(:,1)).*sin(pi*pxy(:,2))','pxy');
    Du = inline('[pi*cos(pi*pxy(:,1)).*sin(pi*pxy(:,2)) pi*sin(pi*pxy(:,1)).*cos(pi*pxy(:,2))]','pxy');
    [elem2dof, elem2edge, edge, bdDof] = dofP3(elem);
    allnode = getallnode_P3(node, elem, edge, elem2dof);
    uI = exactu(allnode);
    N(k) = size(node, 1);
    err(k) = getH1error(node, elem, Du, uI);
    [node, elem] = uniformrefine(node, elem);
end
figure(4);
showrate(N, err);

%% Verify convergence in L2-norm
exactu = @(x) sin(pi * x(:, 1)) .* sin(pi * x(:, 2));
[node, elem] = squaremesh([0, 1, 0, 1], 0.25);
err = zeros(4, 1); N = zeros(4, 1);
for k = 1 : 4
    [node, elem] = uniformrefine(node, elem);
    [elem2dof, elem2edge, edge, bdDof] = dofP3(elem);
    allnode = getallnode_P3(node, elem, edge, elem2dof);
    uI = exactu(allnode);
    N(k) = 1 / sqrt(size(elem, 1));
    err(k) = getL2error(node, elem, exactu, uI);
end
figure(5);
showrateh(N, err);


