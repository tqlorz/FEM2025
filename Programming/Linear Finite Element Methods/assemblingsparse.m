function [A, area] = assemblingsparse(node, elem)
    N = size(node, 1); 
    NT = size(elem, 1);
    i = zeros(9 * NT, 1); 
    j = zeros(9 * NT, 1);
    s = zeros(9 * NT, 1);
    
    index=0;
    for t = 1 : NT
    [At, area(t)] = localstiffness(node(elem(t,:),:));
        for ti = 1 : 3
            for tj = 1 : 3
            index = index + 1;
            i(index) = elem(t, ti);
            j(index) = elem(t, tj);
            s(index) = At(ti, tj);
            end
        end
    end
    A = sparse(i, j, s, N, N);
    area = area';
end