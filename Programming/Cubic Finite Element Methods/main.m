clc;
clear;

%% Set parameters 
h = pow2(-2);
f = @(x) 2 * pi^2 .* sin(pi * x(:, 1)) .* sin(pi * x(:, 2));

%% Solve FEM problem
[node, elem, edge, elem2dof, u_h] = solveFEM(h, f);
allnode = getallnode(node, edge, elem2dof);

%% Draw picture of mesh
figure(1);
hold on;
axis on;
showmesh(node, elem);
findelem(node, elem);       % plot indices of all triangles
findnode(node);             % plot indices of all vertices
findedgedof(node, edge);    % plot indices of all edge

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
% h = pow2(-2);
% err = zeros(4, 1);
% N = zeros(4, 1);
% Du = @(x) [pi * cos(pi * x(:, 1)) .* sin(pi * x(:, 2))
%            pi * sin(pi * x(:, 1)) .* cos(pi * x(:, 2))];
% for k = 1 : 4
%     [node, elem, u_h] = solveFEM(h, f);
%     N(k) = size(node, 1);
%     err(k) = getH1error(node, elem, Du, u_h);
%     h = h / 2;
% end
% figure(3);
% showrate(N,err);

%% Verify convergence in L2-norm
% h = pow2(-1);
% err = zeros(4,1);
% N = zeros(4,1);
% exactu = @(x) sin(pi * x(:, 1)) .* sin(pi * x(:, 2));
% for k = 1 : 4 
%     [node, elem, edge, elem2dof, u_h] = solveFEM(h, f);
%     uI = exactu([node; (node(edge(:,1),:)+node(edge(:,2),:))/2]);
%     N(k) = 1/sqrt(size(elem, 1));
%     err(k) = getL2error(node, elem, exactu, u_h);
%     h = h / 2;
% end
% figure(4);
% showrateh(N, err);

% pde = Stokesdata2;
% [node,elem] = squaremesh([0,1,0,1],0.25);
% err = zeros(4,1); h = zeros(4,1);
% for k = 1:4
%     [node,elem] = uniformrefine(node,elem);
%     [elem2edge,edge] = dofedge(elem);
%     uI = pde.exactu([node; (node(edge(:,1),:)+node(edge(:,2),:))/2]);
%     h(k) = 1/sqrt(size(elem,1));
%     err(k) = getL2error(node,elem,pde.exactu,uI);
% end
% showrateh(h,err);

