clc;
clear;

%% Set parameters 
h = pow2(-2);
f = @(x) (2 * pi^2 + 1) .* sin(pi * x(:, 1)) .* sin(pi * x(:, 2));

%% Solve FEM problem
[node, elem, u_h] = solveFEM(h, f);

%% Draw picture of answer
figure(1);
hold on;
axis on;
showmesh(node, elem);
findelem(node, elem);  % plot indices of all triangles
findnode(node);       % plot indices of all vertices

figure(2);
showresult(node, elem, u_h);

%% Verify convergence in H1-norm
h = pow2(-2);
err = zeros(4, 1);
N = zeros(4, 1);
Du = @(x) [pi*cos(pi * x(:, 1)) .* sin(pi * x(:,2)) pi * sin(pi * x(:, 1)) .* cos(pi * x(:, 2))];
for k = 1 : 4
    [node, elem, u_h] = solveFEM(h, f);
    N(k) = size(node, 1);
    err(k) = getH1error(node, elem, Du, u_h);
    h = h / 2;
end
figure(3);
showrate(N,err);

%% Verify convergence in L2-norm
h = pow2(-1);
err = zeros(4,1);
N = zeros(4,1);
exactu = @(x) sin(pi * x(:, 1)) .* sin(pi * x(:, 2));
for k = 1 : 4 
    [node, elem, u_h] = solveFEM(h, f);
    N(k) = 1/sqrt(size(elem, 1));
    err(k) = getL2error(node, elem, exactu, u_h);
    h = h / 2;
end
figure(4);
showrateh(N, err);