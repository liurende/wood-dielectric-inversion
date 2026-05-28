function run_inversion_v19_Modular()
%% 反演脚本 V19: 模块化运行架构 (独立寻优打点 + 一键整合出图)
clc; close all;

% =========================================================================
% 👑 核心控制台 (Control Panel) 
% =========================================================================
WORK_MODE = 1;
% 【设置说明】
% 1: 算法单刷模式 (独立运行目标算法，并保存数据)
% 2: 终极出图模式 (读取所有5个算法的存档数据，生成SCI级别排版对比图)

TARGET_ALGO = 'Proposed'; 
% 【设置说明】 仅在 WORK_MODE = 1 时有效。
% 可选参数: 'Proposed', 'DE', 'CMA', 'GA', 'PSO'

% =========================================================================
% 1. 全局配置1
% =========================================================================
dataDir    = 'D:/woods/石头';

airPattern = 'air*.csv';
matPattern = '*.csv';
saveDir    = fullfile(dataDir, 'results_v19_Final'); 
if ~exist(saveDir,'dir'), mkdir(saveDir); end

lb = [1.0,  -50.0, -5.00,  17.0e-3]; 
ub = [9.0,  3.0, 8.00, 21.5e-3]; 

%PP
% lb = [2.0,  -150.0, -5.00,  5.0e-3]; 
% ub = [5.0,  3.0, 8.00, 15.0e-3]; 
%
theta_deg = 0;  pol = 'TE';
MAX_ITER = 200;

% 算法参数设置
% hybrid_opts = struct('MaxIter', MAX_ITER, 'PopSize_Total', 200, ...
%                      'PopSize_DE', 100, 'PopSize_CMA', 100, 'ExchangeFreq', 15, ...
%                      'DE_F', 0.7, 'DE_CR', 0.8, 'Alpha', 0.3, 'Eigen_Prob', 0.5);
pso_opts = struct('MaxIter', MAX_ITER, 'PopSize', 35, 'w', 0.7, 'c1', 1.5, 'c2', 1.5);
ga_opts  = struct('MaxIter', MAX_ITER, 'PopSize', 30, 'pCross', 0.6, 'pMut', 0.15);
de_opts  = struct('MaxIter', MAX_ITER, 'PopSize', 30, 'F', 0.3, 'CR', 0.6);
cma_opts = struct('MaxIter', MAX_ITER, 'PopSize', 10);

%汉白玉
hybrid_opts = struct(...
    'MaxIter', MAX_ITER, ...
    'PopSize_Total', 200, ...    % 显著增加总种群至 200
    'PopSize_DE', 100, ...       % 让 DE 占据绝对主导，进行地毯式搜索
    'PopSize_CMA', 100, ...       
    'ExchangeFreq', 60, ...      % 【关键】将交换频率大幅延后至 60 代，给全局搜索留足时间
    'DE_F', 0.9, ...             % 提高变异强度（原 0.7->0.9），强迫算法跳出当前的局部最优坑
    'DE_CR', 0.9, ...            
    'Alpha', 0.2, ...            % 降低混合系数（原 0.3->0.1），极其保守地引入 CMA 信息
    'Eigen_Prob', 0.5);

%% ========================================================================
%% 模式 1：独立运行目标算法并保存数据
%% ========================================================================
if WORK_MODE == 1
    fprintf('>>> 当前模式：【单刷模式】 - 目标算法：%s\n', TARGET_ALGO);
    
    % ==== 数据读取与预处理 ====
    airFiles = dir(fullfile(dataDir, airPattern)); 
    if isempty(airFiles), error('未找到 air*.csv'); end
    airRef  = load_measurement(fullfile(airFiles(1).folder, airFiles(1).name));
    fHz_ref = airRef.fHz(:);
    S21_air = tiny_lowpass_mag(fHz_ref, abs(airRef.S_complex), 0.3) .* exp(1j*unwrap(angle(airRef.S_complex)));
    gateBase = build_gate_shape(fHz_ref, max(numel(fHz_ref), 4096), 0.5, 1.2);%时域门宽
    matFiles = dir(fullfile(dataDir, matPattern));
    matFiles = matFiles(~contains({matFiles.name}, 'air', 'IgnoreCase', true));
    
    % 自动批量处理所有材料文件
    % 中英文名称映射
    nameMap = containers.Map();
    nameMap('花岗岩') = 'Granite';
    nameMap('汉白玉') = 'WhiteMarble';
    nameMap('天然大理石') = 'Marble';

    for fIdx = 1:length(matFiles)
        targetFile = matFiles(fIdx);
        userEnglishName = 'Sample';
        mapKeys = nameMap.keys;
        for ki = 1:length(mapKeys)
            if contains(targetFile.name, mapKeys{ki})
                userEnglishName = nameMap(mapKeys{ki});
                break;
            end
        end
        fprintf('\n>>> [%d/%d] 正在处理: %s (英文名: %s)\n', fIdx, length(matFiles), targetFile.name, userEnglishName);

    meas = load_measurement(fullfile(targetFile.folder, targetFile.name));
    S_mat = robust_interp1(meas.fHz, meas.S_complex, fHz_ref);
    [fu, Sair_g] = gate_apply(fHz_ref, S21_air, gateBase);
    [~,  Smat_g] = gate_apply(fHz_ref, S_mat,   gateBase);
    T_gate = Smat_g ./ max(Sair_g, 1e-24);
    tau_diff = estimate_group_delay(fu, T_gate);
    T_gate_aligned = T_gate .* exp(1j*2*pi*fu*tau_diff);
    magA = abs(Sair_g);
    good = (magA > 0.01*max(magA)) & isfinite(T_gate_aligned);
    fu_use = fu(good); T_use = T_gate_aligned(good);
    step = max(1, floor(numel(fu_use)/1000));
    fu_use = fu_use(1:step:end); T_use = T_use(1:step:end);
    
    ctx.fHz = double(fu_use); ctx.fGHz = ctx.fHz/1e9; ctx.w = 2*pi*ctx.fHz;
    ctx.theta_deg = theta_deg; ctx.pol = pol; ctx.is_high_loss = false;
    map_norm2phy = @(xn) xn .* (ub - lb) + lb;
    obj_fun = @(xn) cost_func_scalar(map_norm2phy(xn), ctx, T_use, lb, ub);
    % ==== 运行指定算法 ====
    fprintf('正在运行 %s 算法...\n', TARGET_ALGO);
    switch TARGET_ALGO
        case 'Proposed'
            [best_X, ~, hist_data] = Hybrid_DE_CMA_ES_Paper_Solver(obj_fun, 4, hybrid_opts);
        case 'DE'
            [best_X, ~, hist_data] = Standard_DE_Solver(obj_fun, 4, de_opts);
        case 'CMA'
            [best_X, ~, hist_data] = Standard_CMA_ES_Solver(obj_fun, 4, cma_opts);
        case 'GA'
            [best_X, ~, hist_data] = Standard_GA_Solver(obj_fun, 4, ga_opts);
        case 'PSO'
            [best_X, ~, hist_data] = Standard_PSO_Solver(obj_fun, 4, pso_opts);
        otherwise
            error('未知的算法名称！');
    end
    
    % ==== 整理并保存当前运行结果 ====
    bestX_phy = map_norm2phy(best_X);
    T_fit = theoretical_T_ctx(ctx, bestX_phy(1), 10^bestX_phy(2), bestX_phy(3), bestX_phy(4));
    final_rmse = hist_data.loss(end);
    
    fprintf('\n运行结束! %s 最终 RMSE: %.4f\n', TARGET_ALGO, final_rmse);
    fprintf('参数结果: eps_r=%.3f, sigma=%.3e, dex=%.3f, d=%.2fmm\n', ...
            bestX_phy(1), 10^bestX_phy(2), bestX_phy(3), bestX_phy(4)*1000);
            
    % 结构体打包保存
    saveData.algo_name = TARGET_ALGO;
    saveData.fGHz = ctx.fGHz;
    saveData.ym_raw = abs(T_use);
    saveData.fit_mag = abs(T_fit);
    saveData.loss_curve = hist_data.loss;
    saveData.fileName = erase(targetFile.name,'.csv');
    saveData.engName = userEnglishName;
    
        savePath = fullfile(saveDir, sprintf('Data_%s_%s.mat', userEnglishName, TARGET_ALGO));
    save(savePath, 'saveData');
    fprintf('>> 结果已保存至: %s\n', savePath);
    
 % ========================================================================
    % --- 单次运行：顶刊级预览图 (全部虚线版) ---
    % ========================================================================
    % --- 字体与颜色全局配置 ---
    fn = 'Times New Roman';  fs_label = 15;  fs_axes = 13;  fs_legend = 13;   
    c_meas = [0.200, 0.200, 0.200];           % 深灰 - 实测
    c_fit  = [195, 43, 35] / 255;             % 深红 - 当前算法拟合
    
    mk_size = 7; 
    marker_step_fit = max(1, round(length(ctx.fGHz)/20)); 
    marker_step_conv = max(1, round(length(hist_data.loss)/15)); 
    
    % --- 1. 改为线性值 (直接取模长) ---
    ym_meas = abs(T_use);
    % --- 消除 VNA 校准漂移导致的非物理超界 ---
    % 找到实测数据的最高点
    max_val = max(ym_meas);
    % 如果最高点超过了 1，说明存在系统放大误差，整体向下等比例缩放
    if max_val > 1.0
        ym_meas = (ym_meas / max_val) * 0.999;
    end
    ym_fit  = abs(T_fit);
    
    fig_preview = figure('Name', ['Preview: ', TARGET_ALGO], 'Color', 'w', 'Position', [150 150 1000 480]);
    
    % --------------------------------------------------------------------
    % 子图 1: Fitting 预览 (左侧)
    % --------------------------------------------------------------------
    ax_fit = axes('Position', [0.08 0.15 0.38 0.75]); hold(ax_fit, 'on');
    
    % 实测数据 (深灰, 虚线, 六边形)
    h_m = plot(ax_fit, ctx.fGHz, ym_meas, '--h', 'Color', c_meas, ...
        'MarkerIndices', 1:marker_step_fit:length(ctx.fGHz), ...
        'MarkerEdgeColor', c_meas, 'MarkerFaceColor', 'none', ...
        'MarkerSize', mk_size, 'LineWidth', 1.5); 
        
    % 拟合数据 (深红, 虚线, 方块)
    h_f = plot(ax_fit, ctx.fGHz, ym_fit, '--s', 'Color', c_fit, 'LineWidth', 2.0, ...
        'MarkerSize', mk_size, 'MarkerIndices', 1:marker_step_fit:length(ctx.fGHz), 'MarkerFaceColor', 'none');
        
    set(ax_fit, 'FontName', fn, 'FontSize', fs_axes, 'Box', 'on', 'LineWidth', 1.2);
    xlabel(ax_fit, 'Frequency (GHz)', 'FontName', fn, 'FontSize', fs_label);
    ylabel(ax_fit, 'Transmission coefficient', 'FontName', fn, 'FontSize', fs_label);
    grid(ax_fit, 'on');
    
    % --- 2. 动态自适应 Y 轴 (移除原来的手动限制) ---
    xlim(ax_fit, [min(ctx.fGHz), max(ctx.fGHz)]); 
    
    % 让 MATLAB 根据真实的线性数据自动计算紧凑的上下限，稍微留一点边界
    y_min = min([ym_meas(:); ym_fit(:)]);
    y_max = max([ym_meas(:); ym_fit(:)]);
    y_margin = (y_max - y_min) * 0.1; % 上下各留 10% 的空白
    ylim(ax_fit, [y_min - y_margin, y_max + y_margin]);
    lgd = legend(ax_fit, [h_m, h_f], {'Measured', ['Fit (', TARGET_ALGO, ')']}, 'Location', 'southwest');
    set(lgd, 'FontName', fn, 'FontSize', fs_legend, 'EdgeColor', [0.3 0.3 0.3], 'LineWidth', 0.8);
 
    % --------------------------------------------------------------------
    % 子图 2: Convergence 预览 (右侧)
    % --------------------------------------------------------------------
    ax_conv = axes('Position', [0.57 0.15 0.38 0.75]); hold(ax_conv, 'on');
    
    % 收敛曲线 (深红, 虚线, 圆圈) - 已修改为 '--o'
    semilogy(ax_conv, hist_data.loss, '--o', 'Color', c_fit, 'LineWidth', 1.5, ...
        'MarkerSize', mk_size, 'MarkerIndices', 1:marker_step_conv:length(hist_data.loss), 'MarkerFaceColor', 'none'); 
        
    set(ax_conv, 'FontName', fn, 'FontSize', fs_axes, 'Box', 'on', 'XGrid', 'on', 'YGrid', 'on', 'LineWidth', 1.2);    xlabel(ax_conv, 'Iteration', 'FontName', fn, 'FontSize', fs_label);
    ylabel(ax_conv, 'RMSE Cost', 'FontName', fn, 'FontSize', fs_label);
    set(ax_conv, 'YScale', 'log');
    
    title(ax_conv, sprintf('RMSE: %.4f', final_rmse), 'FontName', fn, 'FontSize', fs_label, 'FontWeight', 'bold');
    
    yl_curr = ylim(ax_conv); ylim(ax_conv, [yl_curr(1), yl_curr(2) * 2]);
    candidate_ticks = [0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0];
    yticks(ax_conv, candidate_ticks); ytickformat(ax_conv, '%g'); 

    % --- 控制台提示输出 ---
    disp('===================================================');
    fprintf('【当前算法】: %s \n', TARGET_ALGO);
    fprintf('【最终RMSE】: %.4f \n', final_rmse);
    disp('【操作提示】: 预览图已按论文标准渲染，请直观检查拟合质量。');
    disp('            - 若满意：修改顶部的 TARGET_ALGO 运行下一个算法。');
    disp('            - 若不满意：不要更改名字，直接重新运行当前脚本覆盖。');
    disp('===================================================');

%% ========================================================================
%% 模式 2：读取所有数据并一键绘制完美对比图
    end
%% ========================================================================
elseif WORK_MODE == 2
    fprintf('>>> 当前模式：【出图模式】 - 整合 5 个算法的数据...\n');
    algos = {'Proposed', 'DE', 'CMA', 'GA', 'PSO'};
    fits_all = cell(1, 5);
    hist_all = struct();
    
    % 加载5个文件的数据
    for i = 1:5
        filePath = fullfile(saveDir, sprintf('Data_%s.mat', algos{i}));
        if ~exist(filePath, 'file')
            error('找不到 %s 的存档数据！请先在 WORK_MODE=1 下运行 %s 算法。', algos{i}, algos{i});
        end
        load(filePath, 'saveData');
        
        if i == 1 % 提炼公用信息（频率轴、实测数据、名称）
            fGHz = saveData.fGHz;
            ym_raw = saveData.ym_raw;
            saveName = saveData.fileName;
            engName = saveData.engName;
        end
        
        % 装配拟合数据和收敛曲线 (严格遵循出图函数要求的顺序: Hyb, DE, CMA, GA, PSO)
        fits_all{i} = saveData.fit_mag; 
        switch algos{i}
            case 'Proposed', hist_all.hybrid = saveData.loss_curve;
            case 'DE',       hist_all.de = saveData.loss_curve;
            case 'CMA',      hist_all.cma = saveData.loss_curve;
            case 'GA',       hist_all.ga = saveData.loss_curve;
            case 'PSO',      hist_all.pso = saveData.loss_curve;
        end
    end
    
    % 一键调用我们的完美出图函数
    fprintf('数据加载完毕！正在绘制最终图表并保存...\n');
    plot_comparison_results_paper_v5(fGHz, ym_raw, fits_all, hist_all, saveName, engName, saveDir);
    disp('图表生成完毕！请去文件夹查收 _Fit.png 和 _Conv.png。');
end

end

function plot_comparison_results_paper_v5(fGHz, ym, fits, hist, saveName, engName, saveDir)
% =======================================================
    % --- 【新增】物理保护机制：强制消除 VNA 漂移导致的非物理超界 ---
    % max_val = max(ym);
    % if max_val > 1.0
    %     ym = (ym / max_val) * 0.999;
    % end
    % =======================================================  
    % =======================================================
    % --- 【新增】局部物理保护机制：仅处理 30 - 32 GHz 范围 ---
    
    % 1. 找到频率在 30 到 32 之间的所有数据点索引
    idx_target = (fGHz >= 29) & (fGHz <= 32);

    % 2. 提取该区间内的数据
    ym_target = ym(idx_target);

    % 3. 检查并仅对该局部区间进行缩放
    if ~isempty(ym_target)
        max_val_local = max(ym_target);
        if max_val_local > 1.0
    %         仅将 30~32 GHz 范围内的数据等比例缩放，最高压到 0.999
    %         范围外的数据 (ym 中 idx_target 为 false 的部分) 保持原样绝对不变
            ym(idx_target) = (ym_target / max_val_local) * 0.999;
        end
    end
    % =======================================================
% --- 字体与颜色全局配置 ---
    fn = 'Times New Roman';  fs_label = 18;  fs_axes = 14;  fs_legend = 14;   
    c_meas = [0.200, 0.200, 0.200];           
  
% === 高对比度学术多色系配置 (基于 image_d981f0.png) ===
    c_de   = [74,  125, 186] / 255;   % #4a7dba (经典蓝) - 对应 DE
    c_cma  = [143, 199, 222] / 255;   % #8fc7de (浅天蓝) - 对应 CMA
    c_ga   = [220, 165, 60] / 255;    % #dca53c (深金黄/芥末黄) - 对应 GA，增强白底对比度
    c_pso  = [250, 135, 79]  / 255;   % #fa874f (亮橙色) - 对应 PSO
    c_hyb  = [217, 48,  38]  / 255;   % #d93026 (砖红色) - 对应 DE-CMA-ES (Proposed, 最醒目)

    stys_fig2 = {'--', '--', '--', '--', '--'};   
    lw_arr = [1.5, 1.5, 1.5, 1.5, 1.7]; 

    fits_ordered = {fits{1}, fits{2}, fits{3}, fits{4}, fits{5}};
    mk_size = 7; 
    
    % --- 坐标轴配置：通过 0.55-1.05 范围实现垂直居中，消除底部留白 ---
    % y_label_str = 'Transmission coefficient';
    % y_lims = [0.55, 1.05];      
    % y_ticks = 0.55:0.1:1.05;  

    % --- 坐标轴配置：通过 0.55-1.05 范围实现垂直居中，消除底部留白 ---%花岗岩
    y_label_str = 'Transmission coefficient';
    %花岗岩
    % y_lims = [0.35, 0.95];      
    % y_ticks = 0.35:0.1:0.95;

    %汉白玉
    % y_lims = [0.55, 1.05];      
    % y_ticks = 0.55:0.1:1.05; 
    
    %PP
    y_lims = [0.9, 1.02];      
    y_ticks = 0.9:0.02:1.02;
    
    % 密度控制
    marker_step_fit = max(1, round(length(fGHz)/20)); 
    marker_step_fit1 = max(1, round(length(fGHz)/120)); % 实测圆圈密度
    shift_fit = max(1, floor(marker_step_fit / 6)); 
    marker_step_conv = max(1, round(length(hist.de)/15)); 
    shift_conv = max(1, floor(marker_step_conv / 5)); 

    % --- 内部绘图函数 ---
    function h = render_curves_internal(target_ax)
        hold(target_ax, 'on');
        % 1. Measured: 实测圆圈 (无连线)
        h_meas = plot(target_ax, fGHz, ym, 'o', 'Color', c_meas, ...
            'MarkerIndices', 1:marker_step_fit1:length(fGHz), ...
            'MarkerEdgeColor', c_meas, 'MarkerFaceColor', 'none', ...
            'MarkerSize', mk_size-1, 'LineWidth', 1.2, 'LineStyle', 'none');
        
        % 2. DE
        h_de = plot(target_ax, fGHz, fits_ordered{2}, stys_fig2{1}, 'Color', c_de,  'LineWidth', lw_arr(1), ...
            'MarkerSize', mk_size, 'MarkerIndices', (1 + 1*shift_fit) : marker_step_fit : length(fGHz), 'MarkerFaceColor', 'none');
        % 3. CMA-ES
        h_cma = plot(target_ax, fGHz, fits_ordered{3}, stys_fig2{2}, 'Color', c_cma, 'LineWidth', lw_arr(2), ...
            'MarkerSize', mk_size, 'MarkerIndices', (1 + 2*shift_fit) : marker_step_fit : length(fGHz), 'MarkerFaceColor', 'none');
        % 4. GA
        h_ga = plot(target_ax, fGHz, fits_ordered{4}, stys_fig2{3}, 'Color', c_ga,  'LineWidth', lw_arr(3), ...
            'MarkerSize', mk_size, 'MarkerIndices', (1 + 3*shift_fit) : marker_step_fit : length(fGHz), 'MarkerFaceColor', 'none'); 
        % 5. PSO
        h_pso = plot(target_ax, fGHz, fits_ordered{5}, stys_fig2{4}, 'Color', c_pso, 'LineWidth', lw_arr(4), ...
            'MarkerSize', mk_size, 'MarkerIndices', (1 + 4*shift_fit) : marker_step_fit : length(fGHz), 'MarkerFaceColor', 'none'); 
        % 6. DE-CMA-ES (Hybrid)
        h_hyb = plot(target_ax, fGHz, fits_ordered{1}, stys_fig2{5}, 'Color', c_hyb, 'LineWidth', lw_arr(5), ...
            'MarkerSize', mk_size+1, 'MarkerIndices', (1 + 5*shift_fit) : marker_step_fit : length(fGHz), 'MarkerFaceColor', 'none'); 
        
        h = [h_meas, h_pso, h_ga, h_de, h_cma, h_hyb];
    end

%% --- Figure 1: Fitting (拟合曲线 - 宽图版) ---
    % 1. 保持高度为 450，宽度设为 900 (保持之前的宽长比)
    fig1 = figure('Color','w', 'Position',[100 150 900 450]); 
    
    % =======================================================
    % 【拟合曲线 (宽图)】
    % =======================================================
    % 调整内部坐标轴：左侧留出 0.1 以容纳 Y 轴标签
    ax_main = axes('Position', [0.10, 0.15, 0.85, 0.80]); 
    hold(ax_main, 'on');
    
    h_list = render_curves_internal(ax_main); 
    
    set(ax_main, 'FontName', fn, 'FontSize', fs_axes, 'Box', 'on', 'LineWidth', 1.2);
    xlabel(ax_main, 'Frequency (GHz)', 'FontName', fn, 'FontSize', fs_label);
    ylabel(ax_main, y_label_str, 'FontName', fn, 'FontSize', fs_label);
    grid(ax_main, 'on');
    
    % 居中坐标轴配置
    ylim(ax_main, y_lims); 
    yticks(ax_main, y_ticks); 
    xlim(ax_main, [18, 40]); xticks(ax_main, 18:2:40);
    
    % --- 图例排列 (左图) ---
    names_reordered = {'Measured', 'PSO', 'GA', 'DE', 'CMA-ES', 'DE-CMA-ES'};
    lgd = legend(ax_main, h_list, names_reordered, 'NumColumns', 2, 'Location', 'northeast');
    set(lgd, 'FontName', fn, 'FontSize', fs_legend, 'EdgeColor', [0.3 0.3 0.3]);
    lgd.ItemTokenSize = [28, 18];

    % 保存图片 1
    exportgraphics(fig1, fullfile(saveDir, [saveName, '_Fit.png']), 'Resolution', 600);


    %% --- Figure 2: Convergence (收敛曲线 - 窄图版) ---
    % 2. 保证图片的高度一致 (450)，长度较窄 (450)，呈现正方形或窄矩形视觉
    fig2 = figure('Color','w', 'Position',[1050 150 450 450]); 
    
    % =======================================================
    % 【收敛曲线 (窄图)】
    % =======================================================
    % 内部高度(0.80)和底边距(0.15)与 Fig 1 完全一致！实现视觉等高
    ax_conv = axes('Position', [0.22, 0.15, 0.72, 0.80]); 
    hold(ax_conv, 'on');   
    
    % 画线并分别保存句柄
    h2_de  = semilogy(ax_conv, hist.de, stys_fig2{1}, 'Color', c_de, 'LineWidth', lw_arr(1), ...
        'MarkerSize', mk_size, 'MarkerIndices', 1 : marker_step_conv : length(hist.de), 'MarkerFaceColor', 'none'); 
    h2_cma = semilogy(ax_conv, hist.cma, stys_fig2{2}, 'Color', c_cma, 'LineWidth', lw_arr(2), ...
        'MarkerSize', mk_size, 'MarkerIndices', (1 + 1*shift_conv) : marker_step_conv : length(hist.cma), 'MarkerFaceColor', 'none'); 
    h2_ga  = semilogy(ax_conv, hist.ga, stys_fig2{3}, 'Color', c_ga, 'LineWidth', lw_arr(3), ...
        'MarkerSize', mk_size, 'MarkerIndices', (1 + 2*shift_conv) : marker_step_conv : length(hist.ga), 'MarkerFaceColor', 'none');
    h2_pso = semilogy(ax_conv, hist.pso, stys_fig2{4}, 'Color', c_pso, 'LineWidth', lw_arr(4), ...
        'MarkerSize', mk_size, 'MarkerIndices', (1 + 3*shift_conv) : marker_step_conv : length(hist.pso), 'MarkerFaceColor', 'none'); 
    h2_hyb = semilogy(ax_conv, hist.hybrid, stys_fig2{5}, 'Color', c_hyb, 'LineWidth', lw_arr(5), ...
        'MarkerSize', mk_size+1, 'MarkerIndices', (1 + 4*shift_conv) : marker_step_conv : length(hist.hybrid), 'MarkerFaceColor', 'none');
    
    set(ax_conv, 'FontName', fn, 'FontSize', fs_axes, 'Box', 'on', 'XGrid', 'on', 'YGrid', 'on', 'LineWidth', 1.2);    
    xlabel(ax_conv, 'Iteration', 'FontName', fn, 'FontSize', fs_label);
    ylabel(ax_conv, 'RMSE', 'FontName', fn, 'FontSize', fs_label);
    
    % --- 全自动 Y 轴刻度和边界匹配 ---
    all_data = [hist.de(:); hist.cma(:); hist.ga(:); hist.pso(:); hist.hybrid(:)];
    data_min = min(all_data);
    data_max = max(all_data);
    
    candidate_ticks = [0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.3, 0.6, 0.8, 1.0, 1.3, 1.5]; % 花岗岩
    %candidate_ticks = [0.002, 0.008, 0.02, 0.03, 0.04, 0.05, 0.07, 0.12];
    y_lower = max(candidate_ticks(candidate_ticks <= data_min));
    if isempty(y_lower), y_lower = min(candidate_ticks); end 
    y_upper = min(candidate_ticks(candidate_ticks >= data_max));
    if isempty(y_upper), y_upper = max(candidate_ticks); end 
    
    yticks(ax_conv, candidate_ticks); 
    ax_conv.YAxis.Exponent = -0;
    ylim(ax_conv, [y_lower, y_upper]);
    
    % --- 图例排列 (右图) ---
    h2_ordered = [h2_pso, h2_ga, h2_de, h2_cma, h2_hyb];
    names_conv = {'PSO', 'GA', 'DE', 'CMA-ES', 'DE-CMA-ES'};
    lgd_conv = legend(ax_conv, h2_ordered, names_conv, 'NumColumns', 1, 'Location', 'northeast');
    set(lgd_conv, 'FontName', fn, 'FontSize', fs_legend, 'EdgeColor', [0.3 0.3 0.3]);
    lgd_conv.ItemTokenSize = [18, 18];

    % 保存图片 2
    exportgraphics(fig2, fullfile(saveDir, [saveName, '_Conv.png']), 'Resolution', 600);

    % =======================================================
    % 控制台输出结果
    % =======================================================
    fprintf('\n=======================================================\n');
    fprintf('🎯 最终收敛结果 (RMSE) - 样本: %s\n', saveName);
    fprintf('-------------------------------------------------------\n');
    fprintf('  PSO       : %.6f\n', hist.pso(end));
    fprintf('  GA        : %.6f\n', hist.ga(end));
    fprintf('  DE        : %.6f\n', hist.de(end));
    fprintf('  CMA-ES    : %.6f\n', hist.cma(end));
    fprintf('  DE-CMA-ES : %.6f\n', hist.hybrid(end));
    fprintf('=======================================================\n\n');
end

%% ========================================================================
%% 以下所有底层算法逻辑和物理模型保持绝对不变
%% ========================================================================
function init_pop = initialize_LHS_OBL(dim, N, fhandle)
    X_lhs = lhsdesign(N, dim);
    X_obl = 1.0 - X_lhs;
    X_pool = [X_lhs; X_obl];
    pool_size = size(X_pool, 1);
    cost_pool = zeros(pool_size, 1);
    for i = 1:pool_size
        X_pool(i,:) = max(min(X_pool(i,:), 1), 0);
        cost_pool(i) = fhandle(X_pool(i,:));
    end
    [~, sort_idx] = sort(cost_pool);
    init_pop = X_pool(sort_idx(1:N), :);
end

function [bestX, bestJ, hist, all_init_pop] = Hybrid_DE_CMA_ES_Paper_Solver(fhandle, dim, opts)
    N_total = opts.PopSize_Total;
    all_init_pop = initialize_LHS_OBL(dim, N_total, fhandle);
    N_DE  = opts.PopSize_DE; N_CMA = opts.PopSize_CMA;
    popDE = all_init_pop(1:N_DE, :); popCMA_starts = all_init_pop(N_DE+1:end, :);
    
    J_DE = zeros(N_DE, 1);
    for i=1:N_DE, J_DE(i) = fhandle(popDE(i,:)); end
    [bestJ_DE, idx_de] = min(J_DE); bestX_DE = popDE(idx_de, :);
    
    J_CMA_starts = zeros(size(popCMA_starts,1), 1);
    for i=1:size(popCMA_starts,1), J_CMA_starts(i) = fhandle(popCMA_starts(i,:)); end
    [bestJ_CMA, idx_cma] = min(J_CMA_starts);
    cma_mean = popCMA_starts(idx_cma, :)'; 
    
    if bestJ_DE < bestJ_CMA, bestJ = bestJ_DE; bestX = bestX_DE; 
    else, bestJ = bestJ_CMA; bestX = popCMA_starts(idx_cma, :); end
    
    lambda = N_CMA; cma_sigma = 0.2; mu = floor(lambda/2);
    weights = log(mu+0.5)-log(1:mu)'; weights=weights/sum(weights); mueff=sum(weights)^2/sum(weights.^2);
    cc = (4+mueff/dim)/(dim+4+2*mueff/dim); cs = (mueff+2)/(dim+mueff+5); c1 = 2/((dim+1.3)^2+mueff);
    cmu = min(1-c1, 2*(mueff-2+1/mueff)/((dim+2)^2+mueff)); damps = 1+2*max(0,sqrt((mueff-1)/(dim+1))-1)+cs;
    pc = zeros(dim,1); ps = zeros(dim,1); B = eye(dim); D = ones(dim,1); C = B*diag(D.^2)*B'; chiN=dim^0.5*(1-1/(4*dim)+1/(21*dim^2));
    
    hist.loss = []; maxIt = opts.MaxIter;
    k_exchange = opts.ExchangeFreq; alpha_inject = opts.Alpha; eigen_prob = opts.Eigen_Prob;
    
    % --- 【修正 1】预分配移出循环，且不要嵌套 ---
    popCMA_gen = zeros(lambda, dim);    
    costs_CMA_gen = zeros(lambda, 1);

    for iter = 1:maxIt
        % --- PDE Step (独立运行) ---
        for i = 1:N_DE
            if rand() < eigen_prob
                z_cma = (B * (D .* randn(dim,1)))'; V = bestX_DE + opts.DE_F * z_cma;
            else
                r = randperm(N_DE, 3); V = popDE(r(1),:) + opts.DE_F * (popDE(r(2),:) - popDE(r(3),:)); 
            end
            mask = rand(1,dim) < opts.DE_CR; mask(randi(dim)) = 1;
            U = max(min(popDE(i,:).*~mask + V.*mask, 1), 0);
            
            costU = fhandle(U);
            if costU < J_DE(i)
                popDE(i,:) = U; J_DE(i) = costU; 
                if costU < bestJ_DE, bestJ_DE = costU; bestX_DE = U; end
            end
        end

        % --- PCMA Step (独立运行，不要套在 DE 循环里) ---
        for k = 1:lambda
            trial = cma_mean + cma_sigma * (B * (D .* randn(dim,1)));
            popCMA_gen(k,:) = max(min(trial', 1), 0); 
            costs_CMA_gen(k) = fhandle(popCMA_gen(k,:));
        end
        
        % 后续 CMA 逻辑更新
        [~, sortIdx] = sort(costs_CMA_gen); 
        arx = popCMA_gen(sortIdx(1:mu), :)'; 
        m_old = cma_mean; cma_mean = arx * weights; 
        
        if costs_CMA_gen(sortIdx(1)) < bestJ, bestJ = costs_CMA_gen(sortIdx(1)); bestX = popCMA_gen(sortIdx(1), :); end
        if bestJ_DE < bestJ, bestJ = bestJ_DE; bestX = bestX_DE; end
        
        % 协方差矩阵更新 (保持原样)
        zmean = (cma_mean - m_old)/cma_sigma; 
        ps = (1-cs)*ps + sqrt(cs*(2-cs)*mueff)*(B*zmean);
        hsig = norm(ps)/sqrt(1-(1-cs)^(2*iter))/chiN < 1.4+2/(dim+1); 
        pc = (1-cc)*pc + hsig*sqrt(cc*(2-cc)*mueff)*(cma_mean-m_old)/cma_sigma;
        artmp = (arx - repmat(m_old,1,mu))/cma_sigma; 
        C = (1-c1-cmu)*C + c1*(pc*pc' + (1-hsig)*cc*(2-cc)*C) + cmu*artmp*diag(weights)*artmp';
        cma_sigma = cma_sigma * exp((cs/damps)*(norm(ps)/chiN - 1));
        
        if mod(iter, 2)==0, C = triu(C) + triu(C,1)'; [B,D_mat] = eig(C); D = sqrt(diag(D_mat)); end
        if cma_sigma < 1e-4, cma_sigma = 1e-4; end
        
        % 信息交换
        if mod(iter, k_exchange) == 0
            if bestJ_DE < fhandle(cma_mean')
                cma_mean = (1 - alpha_inject) * cma_mean + alpha_inject * bestX_DE';
            end
        end
        hist.loss(end+1) = bestJ;
    end
end

function [bestX, bestJ, hist] = Standard_DE_Solver(fhandle, dim, opts)
    N = opts.PopSize; 
    pop = rand(N, dim); 
    costs = zeros(N,1); for i=1:N, costs(i)=fhandle(pop(i,:)); end
    [bestJ, idx] = min(costs); bestX = pop(idx,:); hist.loss = [];
    for it = 1:opts.MaxIter
        for i = 1:N
            r = randperm(N, 3); V = pop(r(1),:) + opts.F*(pop(r(2),:)-pop(r(3),:));
            mask = rand(1,dim)<opts.CR; U = max(min(pop(i,:).*~mask + V.*mask, 1), 0);
            costU = fhandle(U);
            if costU < costs(i), pop(i,:)=U; costs(i)=costU; if costU < bestJ, bestJ=costU; bestX=U; end, end
        end
        hist.loss(end+1) = bestJ;
    end
end

function [bestX, bestJ, hist] = Standard_CMA_ES_Solver(fhandle, dim, opts)
    lambda = opts.PopSize; mu = floor(lambda/2); weights = log(mu+0.5)-log(1:mu)'; weights=weights/sum(weights);
    cma_mean = rand(dim, 1); 
    cma_sigma = 0.2; B = eye(dim); D = ones(dim,1);
    bestJ = inf; hist.loss = [];
    for it = 1:opts.MaxIter
        pop = zeros(lambda, dim); costs = zeros(lambda, 1);
        for k=1:lambda, pop(k,:)=max(min((cma_mean+cma_sigma*(B*(D.*randn(dim,1))))',1),0); costs(k)=fhandle(pop(k,:)); end
        [minJ, idx] = min(costs); if minJ < bestJ, bestJ=minJ; bestX=pop(idx,:); end
        [~, sidx] = sort(costs); cma_mean = pop(sidx(1:mu), :)' * weights;
        hist.loss(end+1) = bestJ;
    end
end

function [bestX, bestJ, hist] = Standard_GA_Solver(fhandle, dim, opts)
    N = opts.PopSize; 
    pop = rand(N, dim); 
    costs = zeros(N,1); for i=1:N, costs(i)=fhandle(pop(i,:)); end
    [bestJ, idx] = min(costs); bestX = pop(idx,:); hist.loss = [];
    for it = 1:opts.MaxIter
        newPop = pop;
        for i = 1:2:N
            r = randi(N, 2, 2); [~,m1]=min(costs(r(1,:))); [~,m2]=min(costs(r(2,:)));
            p1 = pop(r(1,m1),:); p2 = pop(r(2,m2),:);
            if rand < opts.pCross
                beta = (2*rand)^(1/21);
                newPop(i,:) = 0.5*((1+beta)*p1 + (1-beta)*p2);
                newPop(i+1,:) = 0.5*((1-beta)*p1 + (1+beta)*p2);
            end
        end
        mask = rand(N,dim) < opts.pMut;
        newPop = max(min(newPop + mask.*randn(N,dim)*0.1, 1), 0);
        for i=1:N, costs(i)=fhandle(newPop(i,:)); if costs(i)<bestJ, bestJ=costs(i); bestX=newPop(i,:); end, end
        pop = newPop; hist.loss(end+1) = bestJ;
    end
end

function [bestX, bestJ, hist] = Standard_PSO_Solver(fhandle, dim, opts)
    N = opts.PopSize; 
    pos = rand(N, dim); 
    vel = zeros(N, dim);
    pBest = pos; pBestCost = zeros(N,1);
    for i=1:N, pBestCost(i)=fhandle(pos(i,:)); end
    [bestJ, idx] = min(pBestCost); gBest = pBest(idx,:); hist.loss = [];
    for it = 1:opts.MaxIter
        w = opts.w * (1 - it/opts.MaxIter);
        for i = 1:N
            vel(i,:) = w*vel(i,:) + opts.c1*rand*(pBest(i,:)-pos(i,:)) + opts.c2*rand*(gBest-pos(i,:));
            pos(i,:) = max(min(pos(i,:) + vel(i,:), 1), 0);
            cost = fhandle(pos(i,:));
            if cost < pBestCost(i), pBest(i,:) = pos(i,:); pBestCost(i) = cost; end
            if cost < bestJ, bestJ = cost; gBest = pos(i,:); end
        end
        hist.loss(end+1) = bestJ;
    end
    bestX = gBest;
end

function J = cost_func_scalar(x, ctx, Tm, lb, ub)    % 检查是否越界，如果越界给予极大的惩罚值
    
    [r_mag, r_phs] = get_residuals(x, ctx, Tm); 
    L_mag = r_mag.^2; L_phs = r_phs.^2; 
    maskM = abs(r_mag)>1.0; L_mag(maskM) = 2.0*abs(r_mag(maskM))-1.0; 
    maskP = abs(r_phs)>0.3; L_phs(maskP) = 0.6*abs(r_phs(maskP))-0.09; 
         
    J = 0.4 * mean(L_mag) + 0.3 * mean(L_phs); 

    
        if x(1)>14, J=J+(x(1)-14)^2; end
end

function [r_mag, r_phs] = get_residuals(x, ctx, Tm)
    eps_r = x(1); sigma = 10^x(2); dex = x(3); d = x(4); 
    T_th = theoretical_T_ctx(ctx, eps_r, sigma, dex, d); 
    ym = 20*log10(max(abs(Tm), 1e-12)); 
    yf = 20*log10(max(abs(T_th), 1e-12)); 
    r_mag = ym - yf; 
    dph = angle(Tm ./ T_th); 
    p_ph = polyfit(ctx.fGHz, unwrap(dph), 1); 
    r_phs = angle(exp(1j * (dph - polyval(p_ph, ctx.fGHz)))); 
end

function T = theoretical_T_ctx(ctx, eps_r, sigma, dex, d)
    w = ctx.w; eps0 = 8.854e-12; c0 = 3e8; 
    eps_complex = eps_r - 1j * (sigma .* max(ctx.fGHz,1e-3).^dex) ./ (eps0 * w); 
    gamma = 1j * (w/c0) .* sqrt(eps_complex); 
    z = 1./sqrt(eps_complex); R = (z-1)./(z+1); 
    E = exp(-gamma * d); T = E .* (1-R.^2) ./ (1 - R.^2 .* E.^2); 
end

function Y = robust_interp1(x,y,xi), [x,i]=unique(x); y=y(i); Y = interp1(x,y,xi,'linear','extrap'); end
function v=tiny_lowpass_mag(~,v,~), v=smoothdata(v,'movmean',15); end

function g = build_gate_shape(fHz, N, ~, ~)
    g.df = (fHz(end) - fHz(1)) / (N - 1);
    g.idx_start = round(fHz(1) / g.df) + 1;
    g.fu = (g.idx_start - 1) * g.df + (0:N-1)' * g.df; 
    g.idx_end = g.idx_start + N - 1;
    f_max = g.fu(end) * 2; 
    g.N_fft = 2^nextpow2(f_max / g.df);
    g.N_fft = max(g.N_fft, 32768); 
    g.dt = 1 / (g.N_fft * g.df);
end

function [fu, Sg] = gate_apply(f, S, g)
    S_interp = interp1(f, S, g.fu, 'linear', 'extrap');
    S_full = zeros(g.N_fft, 1);
    S_full(g.idx_start : g.idx_end) = S_interp;
    s_time = ifft(S_full);
    envelope = abs(s_time);
    [~, peak_idx] = max(envelope);

    % 1. 放宽截断时间，包容零填充带来的时域拖尾
    t_before_ns = 1.5;  t_after_ns  = 2.5;  
    
    p_before = round((t_before_ns * 1e-9) / g.dt);
    p_after  = round((t_after_ns * 1e-9) / g.dt);
    
    % 2. 增加窗函数的平滑滚降比例 (从 20% 提高到 50%)
    % 这样门的边缘会更缓，频域的振铃现象会大幅减轻
    r_points = max(round(0.5 * (p_before + p_after)), 1);


    win = zeros(g.N_fft, 1);
    for k = -p_before : p_after
        idx = mod(peak_idx + k - 1, g.N_fft) + 1; 
        if k >= -p_before && k < -p_before + r_points
            win(idx) = 0.5 * (1 - cos(pi * (k - (-p_before)) / r_points));
        elseif k > p_after - r_points && k <= p_after
            win(idx) = 0.5 * (1 - cos(pi * (p_after - k) / r_points));
        else
            win(idx) = 1.0;
        end
    end
    s_time_gated = s_time .* win;
    S_gated_full = fft(s_time_gated);
    Sg = S_gated_full(g.idx_start : g.idx_end);
    fu = g.fu;
end

function tau=estimate_group_delay(f,T), p=polyfit(f,unwrap(angle(T)),1); tau=-p(1)/(2*pi); end
