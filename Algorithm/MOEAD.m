classdef MOEAD < ALGORITHM
% <multi/many> <real/integer/label/binary/permutation>
% Multiobjective evolutionary algorithm based on decomposition
% type --- 1 --- The type of aggregation function

%------------------------------- Reference --------------------------------
% Q. Zhang and H. Li, MOEA/D: A multiobjective evolutionary algorithm based
% on decomposition, IEEE Transactions on Evolutionary Computation, 2007,
% 11(6): 712-731.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2024 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        % function print_res(Algorithm, Problem)
        %     % 检查result属性中是否存在可用的种群
        %     if isempty(Algorithm.result)
        %         disp('没有可用的种群数据');
        %         return;
        %     end
        % 
        %     % 假设最后一个保存的种群代表最新的解集
        %     lastPopulation = Algorithm.result{end};
        % 
        %     % 获取目标函数值
        %     objectiveValues = lastPopulation.objs;
        % 
        %     % 绘制Pareto前沿
        %     figure; % 创建一个新的图形窗口
        %     plot(objectiveValues(:,1), objectiveValues(:,2), 'o', 'MarkerSize', 5, 'MarkerFaceColor', 'b');
        %     xlabel('f_1'); % x轴标签，第一个目标函数
        %     ylabel('f_2'); % y轴标签，第二个目标函数
        %     title('MOEA/D-DE Pareto Front'); % 图标题
        %     grid on; % 添加网格
        % end

        function main(Algorithm,Problem)
            % Algorithm.outputFcn = @print_res;
            %% Parameter setting
            type = Algorithm.ParameterSet(1);

            %% Generate the weight vectors
            [W,Problem.N] = UniformPoint    (Problem.N,Problem.M);
            T = ceil(Problem.N/10);

            %% Detect the neighbours of each solution
            B = pdist2(W,W);
            [~,B] = sort(B,2);
            B = B(:,1:T);

            %% Generate random population
            Population = Problem.Initialization();
            Z = min(Population.objs,[],1);
            iteration = 1;
            
            index = find_index(Population, Problem.N);
            center_vector = Population(index).objs;
            % 调用函数调整参考向量
            W = ScaleReferenceVectors(Problem.N, W, center_vector, 0.2);
            % 可能需要重新计算邻居索引
            B = pdist2(W, W);
            [~, B] = sort(B, 2);
            B = B(:, 1:T);

            %% Optimization
            while Algorithm.NotTerminated(Population)
                iteration = iteration + 1;
                
                if  iteration == 50
                    W = TSL(Problem, Population, Problem.N);
                    %% Detect the neighbours of each solution
                    B = pdist2(W,W);
                    [~,B] = sort(B,2);
                    B = B(:,1:T);
                end
                
                % 每20次迭代调整一次参考向量

                % if iteration == 75%|| iteration == 150
                %     index = find_index(Population, Problem.N);
                %     center_vector = Population(index).dec;
                % 
                %     % 调用函数调整参考向量
                %     W = ScaleReferenceVectors(Problem.N, W, center_vector, 0.5);
                % 
                %     % 可能需要重新计算邻居索引
                %     B = pdist2(W, W);
                %     [~, B] = sort(B, 2);
                %     B = B(:, 1:T);
                % end


                for i = 1 : Problem.N
                    % Choose the parents
                    P = B(i,randperm(size(B,2)));

                    % Generate an offspring
                    Offspring = OperatorGAhalf(Problem,Population(P(1:2)));

                    % Update the ideal point
                    Z = min(Z,Offspring.obj);

                    % Update the neighbours
                    switch type
                        case 1
                            % PBI approach
                            normW   = sqrt(sum(W(P,:).^2,2));
                            normP   = sqrt(sum((Population(P).objs-repmat(Z,T,1)).^2,2));
                            normO   = sqrt(sum((Offspring.obj-Z).^2,2));
                            CosineP = sum((Population(P).objs-repmat(Z,T,1)).*W(P,:),2)./normW./normP;
                            CosineO = sum(repmat(Offspring.obj-Z,T,1).*W(P,:),2)./normW./normO;
                            g_old   = normP.*CosineP + 5*normP.*sqrt(1-CosineP.^2);
                            g_new   = normO.*CosineO + 5*normO.*sqrt(1-CosineO.^2);
                        case 2
                            % Tchebycheff approach
                            g_old = max(abs(Population(P).objs-repmat(Z,T,1)).*W(P,:),[],2);
                            g_new = max(repmat(abs(Offspring.obj-Z),T,1).*W(P,:),[],2);
                        case 3
                            % Tchebycheff approach with normalization
                            Zmax  = max(Population.objs,[],1);
                            g_old = max(abs(Population(P).objs-repmat(Z,T,1))./repmat(Zmax-Z,T,1).*W(P,:),[],2);
                            g_new = max(repmat(abs(Offspring.obj-Z)./(Zmax-Z),T,1).*W(P,:),[],2);
                        case 4
                            % Modified Tchebycheff approach
                            g_old = max(abs(Population(P).objs-repmat(Z,T,1))./W(P,:),[],2);
                            g_new = max(repmat(abs(Offspring.obj-Z),T,1)./W(P,:),[],2);
                    end
                    Population(P(g_old>=g_new)) = Offspring;
                end
            end
        end
    end
end
