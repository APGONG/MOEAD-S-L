function Wi_new = SF(Wi, Vc, ri)
    % 计算哈达玛积和线性组合
    result = Wi * ri + Vc * (1 - ri);
    % 归一化
    Wi_new = result / norm(result);
end

