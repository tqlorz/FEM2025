clc;
clear;

%% Set parameters 
h = pow2(-4);
f = @(x) 24*(x(:,1).^2 .* (1-x(:,1)).^2 + x(:,2).^2 .* (1-x(:,2)).^2) ...
       + 8*(1 - 6*x(:,1) + 6*x(:,1).^2) .* (1 - 6*x(:,2) + 6*x(:,2).^2);

%% Solve FEM problem
[node, elem, edge, elem2dof, u_h] = solveFEM(h, f);
allnode = getallnode(node, edge, elem2dof);

%% Draw picture of mesh
% if h >= pow2(-2)
%     figure(1);
%     hold on;
%     axis on;
%     showmesh(node, elem);
%     findelem(node, elem);       % plot indices of all triangles
%     findnode(node);             % plot indices of all vertices
%     findedgedof(node, edge);    % plot indices of all edge
% end
%% Draw picture of FEM
plotsolution(node, elem, elem2dof, u_h);

%% Verify convergence in H1-norm
% [node, elem] = squaremesh([0, 1, 0, 1], 0.25);
% for k = 1 : 4
%     exactu = inline('sin(pi*pxy(:,1)).*sin(pi*pxy(:,2))','pxy');
%     Du = inline('[pi*cos(pi*pxy(:,1)).*sin(pi*pxy(:,2)) pi*sin(pi*pxy(:,1)).*cos(pi*pxy(:,2))]','pxy');
%     [elem2dof, elem2edge, edge, bdDof] = dofP2(elem);
%     allnode = getallnode_P2(node, elem, edge, elem2dof);
%     uI = exactu(allnode);
%     N(k) = size(node, 1);
%     err(k) = getH1error(node, elem, Du, uI);
%     [node, elem] = uniformrefine(node, elem);
% end
% figure(4);
% showrate(N, err);

%% Verify convergence in L2-norm
% exactu = @(x) sin(pi * x(:, 1)) .* sin(pi * x(:, 2));
% [node, elem] = squaremesh([0, 1, 0, 1], 0.25);
% err = zeros(4, 1); N = zeros(4, 1);
% for k = 1 : 4
%     [node, elem] = uniformrefine(node, elem);
%     [elem2dof, elem2edge, edge, bdDof] = dofP2(elem);
%     allnode = getallnode_P2(node, elem, edge, elem2dof);
%     uI = exactu(allnode);
%     N(k) = 1 / sqrt(size(elem, 1));
%     err(k) = getL2error(node, elem, exactu, uI);
% end
% figure(5);
% showrateh(N, err);