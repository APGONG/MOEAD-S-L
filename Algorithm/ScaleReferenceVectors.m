function adjusted_vectors = ScaleReferenceVectors(N, W, Vc, k)
    % 初始化新的参考向量集
    adjusted_vectors = zeros(size(W));
    %计算每个参考向量和中心参考向量的相似度
    P = zeros(N, 1);
    max_max_diff = 0;
    max_diff = 0;

    display(W);
    display(Vc);
    for i = 1 : N
        max_diff = 0;
        for j = 1 : length(Vc)
            max_diff = max(max_diff, W(i, j) - Vc(j));
        end
        P(i) = max_diff;
        max_max_diff = max(max_diff,max_max_diff);
    end
    for i = 1 : N
        P(i) = P(i) / max_max_diff;
    end
    % 遍历每一个参考向量
    for i = 1:N
        % 计算当前参考向量与中心向量的相似度
        pi = P(i);
        
        % 计算ri，这里假设ri直接等于pi，根据实际需求可能需要调整
        ri = pi^k;
        
        % 调用缩放函数SF，这里将ri传递给SF
        adjusted_vectors(i, :) = SF(W(i, :), Vc, ri);
    end
end
