%% 从已保存的 .mat 文件加载数据并绘图
clc; clear; close all;

% 1. 选择要加载的文件 (手动指定文件名)
[file, path] = uigetfile('*.mat', '请选择要调出的仿真数据文件');
if isequal(file, 0), disp('用户取消了操作'); return; end
load(fullfile(path, file));

% 2. 绘图配置 
fig = figure('Color', 'w', 'Position', [100, 100, 900, 700]);
% 总标题包含中文，保留雅黑防止乱码
sgtitle(['已调出的数据: ', file], 'Interpreter', 'none', 'FontWeight', 'bold', 'FontName', 'Microsoft YaHei');

% 更新子图标题为图 2 样式
labels = {'\epsilon_r', 'c', 'd', 'h'};
y_labels = {'RMSE / CRLB', 'RMSE / CRLB', 'RMSE / CRLB', 'RMSE / CRLB'};

% 自动判断当前加载的是“单点文件”还是“完整的曲线文件”
is_single_point = (length(d_scan_mm) == 1);

for i = 1:4
    subplot(2, 2, i);
    hold on; grid on;
    
    if is_single_point
        % --- 【单点绘制模式】 ---
        p_crlb = plot(d_scan_mm, CRLB_records(i, :), 'ks', 'MarkerSize', 10, 'LineWidth', 2.0);
        p_rmse = plot(d_scan_mm, RMSE_records(i, :), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8, 'LineWidth', 1.5);
        xlim([d_scan_mm - 5, d_scan_mm + 5]); 
    else
        % --- 【多点/原版绘制模式】 ---
        p_crlb = plot(d_scan_mm, CRLB_records(i, :), 'k-', 'LineWidth', 2.0);
        p_rmse = plot(d_scan_mm, RMSE_records(i, :), 'k--.', 'MarkerSize', 15, 'LineWidth', 1.5);
        xlim([min(d_scan_mm) max(d_scan_mm)]);
    end
    
    % 设置坐标轴的西文字体和格式
    set(gca, 'YScale', 'log', 'FontName', 'Times New Roman', 'FontSize', 11, 'Box', 'on');
    
    % 强制为标题和坐标轴标签绑定 Times New Roman
    title(labels{i}, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
    xlabel('Sample Thickness (mm)', 'FontName', 'Times New Roman');
    ylabel(y_labels{i}, 'FontName', 'Times New Roman');
    
    % 动态获取当前 Y 轴范围，并在顶部留出空白 (乘以 3 倍)
    current_ylim = ylim; 
    ylim([current_ylim(1), current_ylim(2) * 3]); 
end

% 提取图例并强制绑定 Times New Roman
lgd = legend(subplot(2,2,1), [p_rmse, p_crlb], {'RMSE', 'CRLB'}, 'Location', 'best');
set(lgd, 'FontName', 'Times New Roman');

disp('>>> 数据调出并绘图完成！');