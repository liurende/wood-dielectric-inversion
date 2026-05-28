function plot_wood_comparison()
%% 木材纹理对比图: 垂直 vs 平行纹理
clc; close all;

resDir = 'D:\woods\lzmwoods\results_Wood';
algos = {'Proposed', 'DE', 'CMA', 'GA', 'PSO'};
algo_labels = {'DE-CMA-ES', 'DE', 'CMA-ES', 'GA', 'PSO'};
nAlgo = length(algos);

% 样品定义与纹理配对
pairs = { ...
    {'LargeWood_Vertical', 'LargeWood_Parallel'}, ...
    {'SmallWood_Vertical', 'SmallWood_Parallel'} ...
    };
pair_labels = {'大木方', '小木块'};

% 字体配置
fn = 'Microsoft YaHei'; fs_label = 16; fs_axes = 13; fs_legend = 12;

% --- 学术配色 ---
c_meas = [0.20, 0.20, 0.20];
c_fit  = [217, 48, 38] / 255;    % 砖红
c_vert = [74, 125, 186] / 255;   % 蓝 - 垂直纹理
c_para = [220, 165, 60] / 255;   % 金黄 - 平行纹理

% --- 加载所有数据 ---
allData = struct();
for p = 1:length(pairs)
    for s = 1:2
        sn = pairs{p}{s};
        snFld = matlab.lang.makeValidName(sn);
        for a = 1:nAlgo
            fpath = fullfile(resDir, sprintf('Data_%s_%s.mat', sn, algos{a}));
            if exist(fpath, 'file')
                ld = load(fpath, 'saveData');
                allData.(snFld).(algos{a}) = ld.saveData;
            end
        end
    end
end

%% ========================================================================
%% Figure 1: 透射拟合对比 (2x2 布局)
%% ========================================================================
fig1 = figure('Color','w', 'Position',[50 100 1200 850]);

for p = 1:2
    snV = pairs{p}{1}; snP = pairs{p}{2};
    snVFld = matlab.lang.makeValidName(snV);
    snPFld = matlab.lang.makeValidName(snP);

    % 用 Proposed 算法结果 (索引1)
    dV = allData.(snVFld).(algos{1});
    dP = allData.(snPFld).(algos{1});

    % --- 子图: 垂直纹理 ---
    subplot(2, 2, (p-1)*2 + 1);
    hold on;
    fGHz = dV.fGHz;
    h_meas = plot(fGHz, dV.ym_raw, 'o', 'Color', c_meas, ...
        'MarkerSize', 4, 'LineWidth', 1.0, 'MarkerFaceColor', 'none');
    h_fit  = plot(fGHz, dV.fit_mag, '-', 'Color', c_fit, 'LineWidth', 1.8);

    set(gca, 'FontName', fn, 'FontSize', fs_axes, 'Box', 'on', 'LineWidth', 1.0);
    xlabel('Frequency (GHz)', 'FontName', fn, 'FontSize', fs_label);
    ylabel('|S_{21}|', 'FontName', fn, 'FontSize', fs_label);
    xlim([18, 44]); grid on;

    title(sprintf('%s - 垂直纹理  (RMSE=%.3f dB)', pair_labels{p}, dV.rmse), ...
        'FontName', fn, 'FontSize', fs_label-1, 'FontWeight', 'bold');
    legend([h_meas, h_fit], {'Measured', 'Fit (DE-CMA-ES)'}, ...
        'FontName', fn, 'FontSize', fs_legend, 'Location', 'southwest');

    % 标注提取参数
    txt = sprintf('\\epsilon_r''=%.2f  tan\\delta=%.3f  d=%.0fmm', ...
        dV.params(1), 10^dV.params(2), dV.params(4)*1000);
    text(0.98, 0.12, txt, 'Units', 'normalized', 'FontName', fn, ...
        'FontSize', fs_axes-1, 'HorizontalAlignment', 'right', ...
        'BackgroundColor', [1 1 1 0.7], 'EdgeColor', [0.5 0.5 0.5]);

    % --- 子图: 平行纹理 ---
    subplot(2, 2, (p-1)*2 + 2);
    hold on;
    fGHz_p = dP.fGHz;
    h_meas2 = plot(fGHz_p, dP.ym_raw, 'o', 'Color', c_meas, ...
        'MarkerSize', 4, 'LineWidth', 1.0, 'MarkerFaceColor', 'none');
    h_fit2  = plot(fGHz_p, dP.fit_mag, '-', 'Color', c_para, 'LineWidth', 1.8);

    set(gca, 'FontName', fn, 'FontSize', fs_axes, 'Box', 'on', 'LineWidth', 1.0);
    xlabel('Frequency (GHz)', 'FontName', fn, 'FontSize', fs_label);
    ylabel('|S_{21}|', 'FontName', fn, 'FontSize', fs_label);
    xlim([18, 44]); grid on;

    title(sprintf('%s - 平行纹理  (RMSE=%.3f dB)', pair_labels{p}, dP.rmse), ...
        'FontName', fn, 'FontSize', fs_label-1, 'FontWeight', 'bold');
    legend([h_meas2, h_fit2], {'Measured', 'Fit (DE-CMA-ES)'}, ...
        'FontName', fn, 'FontSize', fs_legend, 'Location', 'southwest');

    txt = sprintf('\\epsilon_r''=%.2f  tan\\delta=%.3f  d=%.0fmm', ...
        dP.params(1), 10^dP.params(2), dP.params(4)*1000);
    text(0.98, 0.12, txt, 'Units', 'normalized', 'FontName', fn, ...
        'FontSize', fs_axes-1, 'HorizontalAlignment', 'right', ...
        'BackgroundColor', [1 1 1 0.7], 'EdgeColor', [0.5 0.5 0.5]);
end

sgtitle('木材 S_{21} 透射拟合: 垂直纹理 vs 平行纹理', 'FontName', fn, 'FontSize', 18, 'FontWeight', 'bold');
exportgraphics(fig1, fullfile(resDir, 'Wood_Fitting_Comparison.png'), 'Resolution', 600);

%% ========================================================================
%% Figure 2: 纹理方向对比 (相同尺寸, 叠加图)
%% ========================================================================
fig2 = figure('Color','w', 'Position',[50 100 1100 500]);

for p = 1:2
    snV = pairs{p}{1}; snP = pairs{p}{2};
    snVFld = matlab.lang.makeValidName(snV);
    snPFld = matlab.lang.makeValidName(snP);

    dV = allData.(snVFld).(algos{1});
    dP = allData.(snPFld).(algos{1});

    subplot(1, 2, p);
    hold on;

    % 各用自己的频率轴
    h_m_V = plot(dV.fGHz, dV.ym_raw, 'o', 'Color', c_vert, ...
        'MarkerSize', 3, 'LineWidth', 0.8, 'MarkerFaceColor', 'none');
    h_f_V = plot(dV.fGHz, dV.fit_mag, '-', 'Color', c_vert, 'LineWidth', 2.0);

    h_m_P = plot(dP.fGHz, dP.ym_raw, '^', 'Color', c_para, ...
        'MarkerSize', 3, 'LineWidth', 0.8, 'MarkerFaceColor', 'none');
    h_f_P = plot(dP.fGHz, dP.fit_mag, '--', 'Color', c_para, 'LineWidth', 2.0);

    set(gca, 'FontName', fn, 'FontSize', fs_axes, 'Box', 'on', 'LineWidth', 1.0);
    xlabel('Frequency (GHz)', 'FontName', fn, 'FontSize', fs_label);
    ylabel('|S_{21}|', 'FontName', fn, 'FontSize', fs_label);
    xlim([18, 44]); grid on;

    title(sprintf('%s: 垂直 vs 平行纹理', pair_labels{p}), ...
        'FontName', fn, 'FontSize', fs_label, 'FontWeight', 'bold');

    % 图例: 用颜色区分纹理, 用线型区分实测/拟合
    lgd = legend([h_f_V, h_f_P, h_m_V, h_m_P], ...
        {'垂直纹理 Fit', '平行纹理 Fit', '垂直纹理 Meas', '平行纹理 Meas'}, ...
        'FontName', fn, 'FontSize', fs_legend, 'Location', 'southwest', 'NumColumns', 2);

    % 标注参数差异
    d_eps = dV.params(1) - dP.params(1);
    d_tand = 10^dV.params(2) - 10^dP.params(2);
    txt = sprintf(['\\Delta\\epsilon_r'' = %+.2f  ' ...
                   '\\Deltatan\\delta = %+.3f'], d_eps, d_tand);
    text(0.98, 0.90, txt, 'Units', 'normalized', 'FontName', fn, ...
        'FontSize', fs_axes, 'HorizontalAlignment', 'right', ...
        'BackgroundColor', [1 1 1 0.7], 'EdgeColor', [0.3 0.3 0.3]);
end

sgtitle('纹理方向对微波透射的影响', 'FontName', fn, 'FontSize', 18, 'FontWeight', 'bold');
exportgraphics(fig2, fullfile(resDir, 'Wood_Texture_Comparison.png'), 'Resolution', 600);

%% ========================================================================
%% Figure 3: 5算法收敛曲线对比 (4样品)
%% ========================================================================
fig3 = figure('Color','w', 'Position',[50 100 1100 800]);

all_samples = [pairs{1}, pairs{2}];
sample_labels = {'LargeWood-Vertical', 'LargeWood-Parallel', ...
                 'SmallWood-Vertical', 'SmallWood-Parallel'};
colors = {c_fit, [74 125 186]/255, [143 199 222]/255, c_para, [250 135 79]/255};

for i = 1:4
    sn = all_samples{i};
    snFld = matlab.lang.makeValidName(sn);

    subplot(2, 2, i);
    hold on;

    for a = 1:nAlgo
        if isfield(allData.(snFld), algos{a})
            loss = allData.(snFld).(algos{a}).loss_curve;
            semilogy(loss, '-', 'Color', colors{a}, 'LineWidth', 1.5);
        end
    end

    set(gca, 'FontName', fn, 'FontSize', fs_axes-1, 'Box', 'on', ...
        'XGrid', 'on', 'YGrid', 'on', 'LineWidth', 1.0);
    xlabel('Iteration', 'FontName', fn, 'FontSize', fs_label);
    ylabel('RMSE', 'FontName', fn, 'FontSize', fs_label);
    title(sample_labels{i}, 'FontName', fn, 'FontSize', fs_label-1, 'FontWeight', 'bold');
    grid on; grid minor;

    % 标注最终 RMSE
    final_vals = zeros(1, nAlgo);
    for a = 1:nAlgo
        if isfield(allData.(snFld), algos{a})
            lc = allData.(snFld).(algos{a}).loss_curve;
            final_vals(a) = lc(end);
        end
    end
end

lgd = legend(algo_labels, 'FontName', fn, 'FontSize', fs_legend, ...
    'Orientation', 'horizontal', 'Position', [0.32 0.01 0.40 0.03]);
sgtitle('5种算法收敛曲线对比', 'FontName', fn, 'FontSize', 18, 'FontWeight', 'bold');
exportgraphics(fig3, fullfile(resDir, 'Wood_Convergence.png'), 'Resolution', 600);

%% ========================================================================
%% Figure 4: 介电参数条形对比图 (纹理效应一目了然)
%% ========================================================================
fig4 = figure('Color','w', 'Position',[50 100 1000 550]);

param_names = {'\epsilon_r''', 'tan\delta @30GHz', 'n (频散指数)', 'd (mm)'};
param_vert = zeros(4, 4);  % [sample x param]
param_para = zeros(4, 4);

for p = 1:2
    snV = pairs{p}{1}; snP = pairs{p}{2};
    snVFld = matlab.lang.makeValidName(snV);
    snPFld = matlab.lang.makeValidName(snP);

    if isfield(allData, snVFld) && isfield(allData.(snVFld), algos{1})
        dv = allData.(snVFld).(algos{1}).params;
        param_vert(p, :) = [dv(1), 10^dv(2), dv(3), dv(4)*1000];
    end
    if isfield(allData, snPFld) && isfield(allData.(snPFld), algos{1})
        dp = allData.(snPFld).(algos{1}).params;
        param_para(p, :) = [dp(1), 10^dp(2), dp(3), dp(4)*1000];
    end
end

for i = 1:4
    subplot(2, 2, i);
    data = [param_vert(1,i), param_para(1,i); param_vert(2,i), param_para(2,i)];
    hb = bar(data);
    hb(1).FaceColor = c_vert; hb(2).FaceColor = c_para;

    set(gca, 'XTickLabel', pair_labels, 'FontName', fn, 'FontSize', fs_axes, 'Box', 'on');
    ylabel(param_names{i}, 'FontName', fn, 'FontSize', fs_label);
    grid on;

    if i == 1
        legend('Vertical (垂直)', 'Parallel (平行)', ...
            'FontName', fn, 'FontSize', fs_legend, 'Location', 'best');
    end
end

sgtitle('木材介电参数: 纹理方向对比', 'FontName', fn, 'FontSize', 18, 'FontWeight', 'bold');
exportgraphics(fig4, fullfile(resDir, 'Wood_Parameter_BarChart.png'), 'Resolution', 600);

%% ========================================================================
%% 控制台汇总输出
%% ========================================================================
fprintf('\n============ 木材纹理对比汇总 ============\n');
fprintf('%-22s %8s %8s %8s %8s\n', '样品', 'eps_r', 'tand', 'n', 'd(mm)');
fprintf('%-22s %8s %8s %8s %8s\n', '--------------------', '------', '------', '------', '------');
for p = 1:2
    snV = pairs{p}{1}; snP = pairs{p}{2};
    snVFld = matlab.lang.makeValidName(snV);
    snPFld = matlab.lang.makeValidName(snP);
    if isfield(allData, snVFld) && isfield(allData.(snVFld), algos{1})
        dv = allData.(snVFld).(algos{1}).params;
        fprintf('%-22s %8.2f %8.3f %8.2f %8.0f\n', [pair_labels{p} '-垂直'], ...
            dv(1), 10^dv(2), dv(3), dv(4)*1000);
    end
    if isfield(allData, snPFld) && isfield(allData.(snPFld), algos{1})
        dp = allData.(snPFld).(algos{1}).params;
        fprintf('%-22s %8.2f %8.3f %8.2f %8.0f\n', [pair_labels{p} '-平行'], ...
            dp(1), 10^dp(2), dp(3), dp(4)*1000);
    end
    fprintf('\n');
end
fprintf('图片已保存至: %s\n', resDir);

end
