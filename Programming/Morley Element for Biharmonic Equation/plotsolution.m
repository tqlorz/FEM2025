function plotsolution(node, elem, elem2dof, u_h, coefs, n_refine)
    NV = size(node, 1);
    NT = size(elem, 1);


    %% n_refine
    if nargin < 6
        n_refine = 6 - 0.5 * log2(NT/2);
    end
    
    %% calculate the gradbasis
    [D, ~] = gradbasis(node, elem);
    if nargin < 5 || isempty(coefs)
        coefs = getcoefs(D);
    end
    
    %% refine mesh
    node_refined = node;
    elem_refined = elem;

    for k = 1 : n_refine
        [node_refined, elem_refined] = uniformrefine(node_refined, elem_refined);
    end
        
    %% calculate numerical solution at each point
    u_refined = zeros(size(node_refined, 1), 1);
    point_in_triangle = zeros(size(node_refined, 1), 1);  % record which triangle this point belong to 
    
    for i = 1 : size(node_refined, 1)
        x = node_refined(i, :);
        found = false;
        
        for nt = 1 : NT
            % check whether belong to this triangle
            v1 = node(elem(nt, 1), :);
            v2 = node(elem(nt, 2), :);
            v3 = node(elem(nt, 3), :);
            
            [lambda, in_triangle] = point_in_triangle_barycentric(x, v1, v2, v3);
            
            if in_triangle
                point_in_triangle(i) = nt;
                
                u_val = 0;
                for k = 1 : 6
                    % value of basis function
                    phi = coefs(nt, 1, k) * lambda(1)^2 ...
                        + coefs(nt, 2, k) * lambda(2)^2 ...
                        + coefs(nt, 3, k) * lambda(3)^2 ...
                        + coefs(nt, 4, k) * lambda(1) * lambda(2) ...
                        + coefs(nt, 5, k) * lambda(1) * lambda(3) ...
                        + coefs(nt, 6, k) * lambda(2) * lambda(3);
                    
                    % value of point
                    dof_idx = elem2dof(nt, k);
                    u_val = u_val + u_h(dof_idx) * phi;
                end
                
                u_refined(i) = u_val;
                found = true;
                break;
            end
        end
        
        if ~found
            u_refined(k) = NaN;
        end
    end
    
    %% draw pictuce
    figure;
    trisurf(elem_refined, node_refined(:,1), node_refined(:,2), u_refined, ...
            'EdgeColor', 'none', 'FaceColor', 'interp');
    trimesh(elem_refined, node_refined(:,1), node_refined(:,2), u_refined, ...
            'EdgeColor', 'k', 'FaceColor', 'interp');
    xlabel('x'); ylabel('y'); zlabel('u');
    title(['Morley Element Solution (Refined, n=', num2str(n_refine), ')']);
    colorbar;
    
    %% draw picture
    figure;
    trimesh(elem_refined, node_refined(:,1), node_refined(:,2), u_refined, ...
            'EdgeColor', 'k', 'FaceColor', 'interp');
    hold on;
    xlabel('x'); ylabel('y'); zlabel('u');
    title('Solution with Contours');
    colorbar;
    view(2);
end

function [lambda, in_triangle] = point_in_triangle_barycentric(p, v1, v2, v3)
    % calculate in axis of barycentric
    A = [v2-v1; v3-v1]';
    b = (p - v1)';
    
    if rank(A) < 2
        lambda = [1/3, 1/3, 1/3];
        in_triangle = false;
        return;
    end
    
    % solve equation
    lambda23 = A \ b;
    lambda2 = lambda23(1);
    lambda3 = lambda23(2);
    lambda1 = 1 - lambda2 - lambda3;
    
    lambda = [lambda1, lambda2, lambda3];
    
    % check
    in_triangle = (lambda1 >= -1e-10) && (lambda2 >= -1e-10) && (lambda3 >= -1e-10);
end