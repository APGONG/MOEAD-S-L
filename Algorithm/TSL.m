function W_new = TSL(Problem, pop, N)
    %归一化：方便eps设置

    PA = pop;
    Dist = Calculate_Minimum_Distance(PA, 1)
    dmin = 0;
    % 删除旧解
    for i = 1 : N / 10
        %注意这里矩阵对角线是无穷大
        [dist_min, idx] = min(Dist(:));
        [row, col] = ind2sub(size(Dist), idx); % a=row,b=col
        DistA = inf;
        DistB = inf;
        for i = 1 : length(Dist)
            if i~= col
                DistA = min(DistA, Dist(row, i));
            end
            if i~= row
                DistB = min(DistB, Dist(col, i));
            end
        end

        if DistA < DistB
            PA(row) = [];
        else
            PA(col) = [];
        end
        Dist = Calculate_Minimum_Distance(PA, 1);  % 更新距离矩阵
    end
    Dist = Calculate_Minimum_Distance(PA, 0);
    % 添加新解
    sub = N - length(PA);
    it = 1;
    while length(PA) < N
        [dmax, idx] = max(Dist(:));
        [row, col] = ind2sub(size(Dist), idx);
        pos = 0;
        if it <= sub / 2
            pos = row;
        else
            pos = col;
        end
        x = rand();
        % display(PA(row).dec);
        % averDec = [x; f(x)];
        averDec = (PA(row).dec + PA(col).dec) / 2;
        newObj = Problem.CalObj(averDec);  % 计算新的目标函数值
        newCon = Problem.CalCon(averDec);  % 计算新的约束违反度
        newSol = SOLUTION(averDec, newObj, newCon);

        PA = [PA, newSol];
        Dist = Calculate_Minimum_Distance(PA, 0);  % 更新距离矩阵  
        it = it + 1;
    end
    W_new = zeros(size(PA, 1), Problem.M);
    for i = 1 : N
        W_new(i, :) = PA(i).objs / norm(PA(i).objs);
    end
end

function Dist = Calculate_Minimum_Distance(PA, flag)
    % 计算解决方案间的欧氏距离矩阵
    numSolutions = length(PA);
    Dist = zeros(numSolutions);
    for i = 1:numSolutions
        for j = 1:numSolutions
            if i ~= j
                Dist(i, j) = norm(PA(i).dec - PA(j).dec);
            else
                if flag == 1
                    Dist(i, j) = inf;
                else
                    Dist(i, j) = 0;
                end
            end
        end
    end
end

function y = f(x)
    y = 1 - sqrt(x);
end
