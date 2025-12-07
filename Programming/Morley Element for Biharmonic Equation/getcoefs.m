function coefs = getcoefs(D)
    NT = size(D, 1);

    D1(:, :) = D(:, :, 1);
    D2(:, :) = D(:, :, 2);
    D3(:, :) = D(:, :, 3);
    
    for nt = 1 : NT
        d1(nt) = 1 ./ norm(D1(nt, :));
        d2(nt) = 1 ./ norm(D2(nt, :));
        d3(nt) = 1 ./ norm(D3(nt, :));

        n1(nt, :) = -d1(nt) .* D1(nt, :);
        n2(nt, :) = -d2(nt) .* D2(nt, :);
        n3(nt, :) = -d3(nt) .* D3(nt, :);    
    end
    
    for nt = 1 : NT
        A = [1,0,0,0,0,0;
             0,1,0,0,0,0;
             0,0,1,0,0,0;
             0,dot(n1(nt,:),D2(nt,:)),dot(n1(nt,:),D3(nt,:)), 0.5*dot(n1(nt,:),D1(nt,:)),0.5*dot(n1(nt,:),D1(nt,:)),0.5*dot(n1(nt,:),D2(nt,:)+D3(nt,:));
             dot(n2(nt,:),D1(nt,:)),0,dot(n2(nt,:),D3(nt,:)), 0.5*dot(n2(nt,:),D2(nt,:)),0.5*dot(n2(nt,:),D1(nt,:)+D3(nt,:)),0.5*dot(n2(nt,:),D2(nt,:));
             dot(n3(nt,:),D1(nt,:)),dot(n3(nt,:),D2(nt,:)),0, 0.5*dot(n3(nt,:),D1(nt,:)+D2(nt,:)),0.5*dot(n3(nt,:),D3(nt,:)),0.5*dot(n3(nt,:),D3(nt,:))];
        coefs(nt, :, :) = inv(A);
    end
end