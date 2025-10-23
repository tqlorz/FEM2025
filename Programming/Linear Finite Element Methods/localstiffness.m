function [At, area] = localstiffness(p)
    A1t = zeros(3, 3);
    A2t = zeros(3,3);
    B = [p(1, :) - p(3, :); p(2, :) - p(3, :)];
    G = [[1, 0]', ...
         [0, 1]', ...
         [-1, -1]'];
    area = 0.5 * abs(det(B));
    
    for i = 1 : 3
        for j = 1 : 3
            A1t(i, j) = area * ((B \ G(:, i))' * (B \ G(:, j)));
        end
    end
    for i = 1 : 3
        for j = 1 : 3
            if i == j
                A2t(i, j) = area * (0.5^2 + 0.5^2) / 3;    
            else
                A2t(i, j) = area * (0.5^2) / 3;    
            end
        end
    end
    At = A1t + A2t;
end