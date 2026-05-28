function run_Weight_GridSearch(measData, ctx, test_angles)
    % DE-CMA-ES 混合算法目标函数权重网格搜索 (Grid Search)
    % 这是一个独立函数，不会污染主程序的工作区
    
    % --- 1. 基础参数设置 (保持与你主程序一致) ---
    dim = 4; % 优化维度
    opts.PopSize_Total = 100;
    opts.PopSize_DE    = 50;
    opts.PopSize_CMA   = 50;
    opts.MaxIter       = 100;
    opts.Eigen_Prob    = 0.3;
    opts.DE_F          = 0.5;
    opts.DE_CR         = 0.8;
    opts.ExchangeFreq  = 10;
    opts.Alpha         = 0.1;

    % --- 2. 权重网格定义 ---
    w_mag_list = 0 : 0.1 : 1.0;
    num_tests = length(w_mag_list);
    results = zeros(num_tests, 6); % [w_mag, w_phs, eps_r, d_m, RMSE_mag, RMSE_phs]

    fprintf('\n==================================================\n');
    fprintf('开始目标函数权重网格搜索验证 (Pareto Front)...\n');
    
    for i = 1:num_tests
        w_mag = w_mag_list(i);
        w_phs = 1.0 - w_mag;
        
        fprintf('进度 [%2d/%d] : 评估 w_mag = %.1f, w_phs = %.1f ... ', i, num_tests, w_mag, w_phs);
        
        % 构造适应度函数句柄
        fhandle = @(x) compute_normalized_fitness(x, test_angles, measData, ctx, w_mag, w_phs);
        
        % 调用核心优化器
        [bestX, ~, ~, ~] = Hybrid_DE_CMA_ES_Paper_Solver(fhandle, dim, opts);
        
        % 计算独立的误差用于画图
        [~, rmse_mag, rmse_phs] = compute_normalized_fitness(bestX, test_angles, measData, ctx, w_mag, w_phs);
        
        eps_inv = bestX(1);
        d_inv   = bestX(4);
        
        results(i, :) = [w_mag, w_phs, eps_inv, d_inv, rmse_mag, rmse_phs];
        fprintf('完成! (eps: %.2f, d: %.2f)\n', eps_inv, d_inv);
    end

    % --- 3. 绘制帕累托前沿图 ---
    figure('Color', 'w', 'Position', [100, 100, 900, 400], 'Name', 'Weight Grid Search');

    subplot(1,2,1); hold on; grid on; box on;
    plot(results(:,1), results(:,3), '-s', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', [61, 103, 160]/255);
    xlabel('幅度权重 w_{mag}', 'FontName', 'Times New Roman', 'FontSize', 12); 
    ylabel('反演得到的介电常数 \epsilon_r', 'FontName', 'Times New Roman', 'FontSize', 12);
    title('权重对反演结果的敏感度分析', 'FontName', 'SimHei', 'FontSize', 14);

    subplot(1,2,2); hold on; grid on; box on;
    scatter(results(:,5), results(:,6), 120, results(:,1), 'filled', 'MarkerEdgeColor', 'k');
    colormap('jet'); cb = colorbar; ylabel(cb, '幅度权重 w_{mag}', 'FontName', 'Times New Roman');
    xlabel('幅度误差 RMSE_{mag} (归一化)', 'FontName', 'Times New Roman', 'FontSize', 12); 
    ylabel('相位误差 RMSE_{phs} (归一化)', 'FontName', 'Times New Roman', 'FontSize', 12);
    title('Pareto Front (寻找左下角最优点)', 'FontName', 'SimHei', 'FontSize', 14);
    
    for i = 1:num_tests
        text(results(i,5), results(i,6)+0.02, sprintf('w_m=%.1f', results(i,1)), ...
            'FontSize', 9, 'HorizontalAlignment', 'center', 'FontName', 'Times New Roman');
    end
end

%% =========================================================
%% 附属函数：包含归一化的适应度计算器
%% =========================================================
function [J_total, rmse_mag, rmse_phs] = compute_normalized_fitness(x, angles, measData, ctx, w_mag, w_phs)
    % 映射回物理范围 (假设你主程序里有 map_norm2phy)
    % x_phy = map_norm2phy(x);  % 如果你在外层已经处理了映射，这里可以省略
    
    % 调用正向模型计算 S11
    simRes = calculate_S11_response(x, ctx, angles); 
    
    L_mag_all = [];
    L_phs_all = [];
    
    for i = 1:length(angles)
        angStr = sprintf('deg%d', angles(i));
        
        sim_mag = simRes.(angStr); 
        meas_mag = measData.(angStr).mag; 
        
        % 归一化幅度误差
        err_mag = abs(sim_mag - meas_mag) ./ (abs(meas_mag) + 1e-6); 
        L_mag_all = [L_mag_all; err_mag(:)];
        
        % 归一化相位误差 (如果你的数据里有相位的话，没有的话可以注释掉这两行)
        sim_phs = simRes.([angStr '_phs']);
        meas_phs = measData.(angStr).phs;
        err_phs = abs(sim_phs - meas_phs) ./ (abs(meas_phs) + 1e-6);
        L_phs_all = [L_phs_all; err_phs(:)];
    end
    
    rmse_mag = mean(L_mag_all);
    rmse_phs = mean(L_phs_all); % 如果没有相位数据，这里直接写 rmse_phs = 0;
    
    J_total = w_mag * rmse_mag + w_phs * rmse_phs;
end