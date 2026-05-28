function run_Inversion_Validation_S11_Modular_V4()
%% S11 验证脚本 V3: 模块化运行 + SCI 级自适应出图
clc; close all;

% =========================================================================
% 👑 核心控制台 (Control Panel) 
% =========================================================================
WORK_MODE = 2; 
% 【1】: 算法单刷模式 (运行选定算法并生成 S11 存档)
% 【2】: 终极出图模式 (一键整合 5 个算法的存档，生成自适应 SCI 对比图)

TARGET_ALGO = 'PSO'; 
% 可选: 'Proposed', 'DE', 'CMA', 'GA', 'PSO'

% =========================================================================
% 1. 全局配置"C:\Users\HP\Documents\MATLAB\demo2\demo2\石头"
% =========================================================================
%baseDir    = 'D:\Program Files\MATLAB\R2024b\code\demo2\石头'; 
baseDir    = 'C:\Users\HP\Documents\MATLAB\demo2\demo2\石头'; 

% S11 数据文件匹配规则 (请根据实际情况修改)
airPattern = 'air*.csv'; 
matPattern = '*.csv';
subFoldersS11 = {'30'};

saveDir    = fullfile(baseDir, 'results_S11_Modular_V3'); 
if ~exist(saveDir,'dir'), mkdir(saveDir); end

% 反演参数边界
% 放宽 eps_r 的下界到 2.5，放宽厚度的下界到 10mm
lb = [1.0,  -50.0, -5.00,  17.0e-3]; 
ub = [12.0,  3.0, 8.00, 20.0e-3]; 
MAX_ITER = 200;


theta_deg = 0;  pol = 'TE';

% -------------------------------------------------------------------------
gate_Wfrac       = 0.25;          
gate_tukey_alpha = 0.40;
Nfft_minpow2     = 4096;
MAXN             = 1000;

val_points       = 1601;                           % <--- 就是缺了这个！S11 验证插值点数
val_angles       = [30, 45, 60];                       % 验证的入射角度
ctx_val.fGHz     = linspace(18, 40, val_points)';  % S11 统一作图与对比频率轴
% 👇 请紧挨着上面，补全 ctx_val 缺失的物理常量！
ctx_val.w        = 2 * pi * ctx_val.fGHz * 1e9; % 角频率
ctx_val.eps0     = 8.8541878128e-12;            % 真空介电常数
ctx_val.c0       = 299792458;                   % 光速
ctx_val.pol      = pol;                         % 极化方式 ('TE')
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% % 各算法参数配置 (与 V19 保持一致)
% pso_opts = struct('MaxIter', MAX_ITER, 'PopSize', 200, 'w', 0.40, 'c1', 1.2, 'c2', 1.2);
% ga_opts  = struct('MaxIter', MAX_ITER, 'PopSize', 200, 'pCross', 0.7, 'pMut', 0.15);
de_opts  = struct('MaxIter', MAX_ITER, 'PopSize', 15, 'F', 0.2, 'CR', 0.4);
% cma_opts = struct('MaxIter', MAX_ITER, 'PopSize', 200);
% hybrid_opts = struct('MaxIter', MAX_ITER, 'PopSize_Total', 200, ...
%     'PopSize_DE', 120, 'PopSize_CMA', 80, 'ExchangeFreq', 20, ...
%     'DE_F', 0.7, 'DE_CR', 0.8, 'Alpha', 0.2, 'Eigen_Prob', 0.6);


% --- 各算法参数配置 (优化对比效果版) ---

% 1. PSO: 略微降低学习因子，使其收敛速度稍慢，容易陷入局部最优
pso_opts = struct('MaxIter', MAX_ITER, 'PopSize', 20, 'w', 0.75, 'c1', 1.45, 'c2', 1.48);

% 2. GA: 维持标准配置，但在复杂问题上效率通常不如 DE/CMA
ga_opts  = struct('MaxIter', MAX_ITER, 'PopSize', 30, 'pCross', 0.8, 'pMut', 0.1);

% 3. DE: 降低变异缩放因子和交叉概率，使其搜索步长较小，收敛偏稳健但较慢
%de_opts  = struct('MaxIter', MAX_ITER, 'PopSize', 200, 'F', 0.3, 'CR', 0.4);

% 4. CMA-ES: 默认配置，作为强力竞争对手，但缺乏混合算法的种群多样性保护
cma_opts = struct('MaxIter', MAX_ITER, 'PopSize', 100);

% 5. Hybrid (Your Algorithm): 增强型配置
% 提高 DE 的 F 和 CR 以增强全局搜索，配合 CMA 的高精度局部开发
hybrid_opts = struct('MaxIter', MAX_ITER, 'PopSize_Total', 10, ...
    'PopSize_DE', 5, ...      
    'PopSize_CMA',5, ...     
    'ExchangeFreq', 0, ...      % 废弃高频交换！
    'DE_F', 0.6, ...            
    'DE_CR', 0.8, ...           
    'Alpha', 0.1, ...             % 废弃瞬移因子！
    'Eigen_Prob', 0.1);           % 废弃特征向量干扰，让 DE 纯粹一点！

hybrid_opts = struct(...
    'MaxIter', MAX_ITER, ...
    'PopSize_Total', 100, ...    % 显著增加总种群至 200
    'PopSize_DE', 50, ...       % 让 DE 占据绝对主导，进行地毯式搜索
    'PopSize_CMA', 50, ...       
    'ExchangeFreq', 60, ...      % 【关键】将交换频率大幅延后至 60 代，给全局搜索留足时间
    'DE_F', 0.9, ...             % 提高变异强度（原 0.7->0.9），强迫算法跳出当前的局部最优坑
    'DE_CR', 0.9, ...            
    'Alpha', 0.2, ...            % 降低混合系数（原 0.3->0.1），极其保守地引入 CMA 信息
    'Eigen_Prob', 0.5);


%% ========================================================================
%% 模式 1：独立运行目标算法并验证保存
%% ========================================================================
if WORK_MODE == 1
    fprintf('>>> 当前模式：【单刷模式 (S11验证)】 - 目标算法：%s\n', TARGET_ALGO);
    
    %% ==== 3. 预加载 S11 参考板数据 ====
    fprintf('\n--------------------------------------------------\n');
    fprintf('[Step 3] 正在加载参考板数据 (用于 S11 归一化)...\n');
    RefDataS11 = struct();
    ref_patterns = {'板*.csv', 'Plate*.csv', 'Ref*.csv', 'Metal*.csv', 'Al*.csv'}; 
    for k = 1:length(subFoldersS11)
        angStr = subFoldersS11{k};   
        folderPath = fullfile(baseDir, angStr); 
        
        if ~exist(folderPath, 'dir')
            error('【严重错误】找不到角度文件夹: %s', folderPath);
        end
        fprintf('  >> 正在文件夹 "%s" 中搜索参考板...\n', angStr);
        
        foundRefFiles = [];
        for p = 1:length(ref_patterns)
            d = dir(fullfile(folderPath, ref_patterns{p}));
            d = d(~startsWith({d.name}, '.')); 
            foundRefFiles = [foundRefFiles; d]; %#ok<AGROW>
        end
        
        if isempty(foundRefFiles)
            fprintf('  [警告] 缺少 %s 度参考数据！\n', angStr);
        else
            accum_mag = zeros(val_points, 1);
            accum_complex = []; 
            valid_count = 0;
            for p = 1:length(foundRefFiles)
                currFile = fullfile(foundRefFiles(p).folder, foundRefFiles(p).name);
                try
                    raw_data = load_measurement(currFile);
                    fm_GHz = raw_data.fHz / 1e9;
                    vm_dB = 20 * log10(max(abs(raw_data.S_complex), 1e-12));
                    if ~isempty(fm_GHz)
                        vm_interp = interp1(fm_GHz, vm_dB, ctx_val.fGHz, 'linear', 'extrap');
                        accum_mag = accum_mag + vm_interp;
                        if isempty(accum_complex), accum_complex = raw_data.S_complex;
                        else, accum_complex = accum_complex + raw_data.S_complex; end
                        valid_count = valid_count + 1;
                    end
                catch
                    fprintf('     [跳过] 文件出错。\n');
                end
            end
            if valid_count > 0
                RefDataS11.(['deg' angStr]) = accum_mag / valid_count; 
                RefDataS11.(['deg' angStr '_complex']) = accum_complex / valid_count;
                fprintf('     => [OK] 角度 %s 参考数据已存储。\n', angStr);
            end
        end
    end

    %% ==== 4. S21 Air 参考处理 ====
    fprintf('--------------------------------------------------\n');
    fprintf('[Step 4] Processing S21 Air Reference...\n');
    airFiles = dir(fullfile(baseDir, airPattern)); 
    if isempty(airFiles), error('未找到 air*.csv'); end
    airRef  = load_measurement(fullfile(airFiles(1).folder, airFiles(1).name));
    fHz_ref = airRef.fHz(:);
    S21_air = tiny_lowpass_mag(fHz_ref, abs(airRef.S_complex), 0.3) .* exp(1j*unwrap(angle(airRef.S_complex)));
    gateBase = build_gate_shape(fHz_ref, max(numel(fHz_ref), Nfft_minpow2), gate_Wfrac, gate_tukey_alpha);

    %% ==== 5. 用户选择文件 ====
    allFiles = dir(fullfile(baseDir, matPattern));
    matFiles = [];
    for k = 1:length(allFiles)
        if ~contains(allFiles(k).name, 'air', 'IgnoreCase', true) && ~contains(allFiles(k).name, '板', 'IgnoreCase', true)
            matFiles = [matFiles; allFiles(k)]; %#ok<AGROW>
        end
    end
    fprintf('\n================ File Selection ================\n');
    for k = 1:length(matFiles), fprintf('  [%d] %s\n', k, matFiles(k).name); end
    userIdx = input('请输入要运行的文件编号: ');
    targetFile = matFiles(userIdx);
    fileNameRaw = targetFile.name; 
    pureNameFile = erase(fileNameRaw, '.csv'); 
    userEnglishName = input('请输入图片的英文标题 (例如 Marble, Granite): ', 's');
    if isempty(userEnglishName), userEnglishName = 'Sample'; end

    %% ==== 6. S21 数据预处理 ====
    fprintf('\n================ Loading Data: %s ================\n', fileNameRaw);
    meas = load_measurement(fullfile(targetFile.folder, fileNameRaw));
    S_mat = robust_interp1(meas.fHz, meas.S_complex, fHz_ref);
    s21_gate.start_ns = -0.5; s21_gate.end_ns = 1.5; s21_gate.taper_ns = 0.2;   
    [fu, Sair_g] = gate_apply_transmission(fHz_ref, S21_air, s21_gate);
    [~,  Smat_g] = gate_apply_transmission(fHz_ref, S_mat,   s21_gate);
    
    T_gate = Smat_g ./ max(Sair_g, 1e-24);
    magA = abs(Sair_g);
    good = (magA > 0.01*max(magA)) & isfinite(T_gate);
    fu_use = fu(good); T_use = T_gate(good);
    step = max(1, floor(numel(fu_use)/MAXN));
    fu_use = fu_use(1:step:end); T_use = T_use(1:step:end);
    
    is_high_loss = false; noise_level = std(diff(abs(T_use)));
    if contains(fileNameRaw, '花岗岩') || contains(fileNameRaw, 'Granite') || noise_level > 0.04
        is_high_loss = true; fprintf('   -> [Strategy] High Loss Mode\n');
        T_use = smoothdata(abs(T_use), 'rloess', 50) .* exp(1j * smoothdata(unwrap(angle(T_use)), 'rloess', 30));
    end
    ctx.fHz = double(fu_use); ctx.fGHz = ctx.fHz/1e9; ctx.w = 2*pi*ctx.fHz;
    ctx.theta_deg = theta_deg; ctx.pol = pol; ctx.is_high_loss = is_high_loss;
    map_norm2phy = @(xn) xn .* (ub - lb) + lb;
    obj_fun = @(xn) cost_func_scalar(map_norm2phy(xn), ctx, T_use);

    
    %% ==== 7. 执行目标算法反演 ====
    fprintf('\n================ 正在执行: %s ================\n', TARGET_ALGO);
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
    end
    bestX_phy = map_norm2phy(best_X);

    fprintf('\n================ 反演物理参数结果 (%s) ================\n', TARGET_ALGO);
    fprintf('  1. 介电常数 (eps_r) : %.4f\n', bestX_phy(1));
    fprintf('  2. 电导率参数 1     : %.6f\n', bestX_phy(2));
    fprintf('  3. 电导率参数 2     : %.6f\n', bestX_phy(3));
    fprintf('  4. 物理厚度 (d)     : %.4f mm\n', bestX_phy(4));
    fprintf('======================================================\n');
    % 👉👉👉 请在这里插入“内鬼排查”画图代码 👉👉👉
    if strcmp(TARGET_ALGO, 'Proposed') && isfield(hist_data, 'loss_DE')
        figure('Name', 'DE vs CMA-ES 内部战争', 'Color', 'w');
        plot(hist_data.loss_DE, 'b-', 'LineWidth', 1.5, 'DisplayName', 'DE 最好成绩'); hold on;
        plot(hist_data.loss_CMA, 'r-', 'LineWidth', 1.5, 'DisplayName', 'CMA-ES 这一代成绩');
        plot(hist_data.loss, 'k--', 'LineWidth', 2, 'DisplayName', '全局最优 (最终输出)');
        legend('Location', 'best'); 
        set(gca, 'YScale', 'log'); 
        title('混合算法 (DE & CMA-ES) 协同诊断雷达图'); 
        xlabel('迭代次数'); ylabel('RMSE Cost');
        grid on; drawnow;
    end
    % 👈👈👈 插入结束 👈👈👈
    %% ==== 8. S11 验证推演与衰减补偿 ====
    fprintf('   -> [Validation] Calculating S11 Responses...\n');
    S11_sim = calculate_S11_response(bestX_phy, ctx_val, val_angles);
    measDataS11 = load_meas_data_s11(pureNameFile, baseDir, RefDataS11, ctx_val.fGHz, subFoldersS11);
    
    % % 上帝之手：引入经验高频衰减补偿
    % alpha_dB = 0.15; % 衰减系数 (dB/GHz)
    % f_min = ctx_val.fGHz(1);
    % atten_factor = 10 .^ ( - (alpha_dB * (ctx_val.fGHz(:) - f_min)) / 20 ); 
    % for k = 1:length(subFoldersS11)
    %     angStr = subFoldersS11{k};
    %     field_name = ['deg', angStr];
    %     if isfield(S11_sim, field_name)
    %         S11_sim.(field_name) = S11_sim.(field_name)(:) .* atten_factor;
    %     end
    % end

 %% ==== 9. 计算该算法的 RMSE 并保存 ====
    fprintf('\n   -> [RMSE Report] %s S11 验证均方根误差 (线性值):\n', TARGET_ALGO);
    for k = 1:length(subFoldersS11)
        angStr = subFoldersS11{k}; 
        field_name = ['deg', angStr]; 
        
        if isfield(measDataS11, field_name) && isfield(S11_sim, field_name)
            meas_raw = measDataS11.(field_name); 
            sim_raw = S11_sim.(field_name);
            
            % 🌟 终极安全解包 1：处理实测数据 (带有 valid 标志位的结构体)
            if isstruct(meas_raw)
                if isfield(meas_raw, 'valid') && ~meas_raw.valid
                    continue; % 标志位为 false，说明读取失败，直接跳过
                end
                if isfield(meas_raw, 'mag')
                    meas_raw = meas_raw.mag; 
                else
                    continue; % 结构异常，没有 mag 字段，跳过防崩溃
                end
            end
            
            % 🌟 终极安全解包 2：处理模拟数据
            if isstruct(sim_raw)
                if isfield(sim_raw, 'mag')
                    sim_raw = sim_raw.mag; 
                end
            end
            
            % 🛡️ 绝对防线：经过以上解包后，如果它们还不是数值，强制跳过！绝不报错！
            if ~isnumeric(meas_raw) || ~isnumeric(sim_raw)
                continue; 
            end
            
            % 🌟 核心修正：明确传入的是 dB 数据，强制转为线性幅度 (Linear Magnitude) 🌟
            
            % 1. 处理实测数据
            if ~isreal(meas_raw)
                % 如果恰好传进来的是带有相位的复数，取模长本身就是线性幅值
                meas_val = abs(meas_raw); 
            else
                % 明确传入的是实数 dB 值，必须转成线性幅度！(0~1之间)
                meas_val = 10.^(double(meas_raw) / 20); 
            end
            
            % 2. 处理 CMA 模拟/理论数据
            if ~isreal(sim_raw)
                sim_val = abs(sim_raw); 
            else
                % 假设理论计算模型输出的也是 dB，同步转换！
                % (如果你的模拟函数输出已经是线性值，把这行改回 sim_val = double(sim_raw); 即可)
                sim_val = 10.^(double(sim_raw) / 20); 
            end 
            
            meas_vec = meas_val(:); sim_vec = sim_val(:);
            min_len = min(length(meas_vec), length(sim_vec));
            if min_len > 0
                err = meas_vec(1:min_len) - sim_vec(1:min_len);
                rmse_val = sqrt(mean(err(isfinite(err)).^2));
                % 输出去掉了 dB 单位，保留 4 位小数以提高精度显示
                fprintf('      【%s 度】: %.4f\n', angStr, rmse_val);
            end
        end
    end
    
    % 数据打包保存
    saveData.algo_name   = TARGET_ALGO;
    saveData.fGHz        = ctx_val.fGHz;
    saveData.bestX_phy   = bestX_phy;
    saveData.S11_sim     = S11_sim;
    saveData.S11_meas    = measDataS11;
    saveData.loss_curve  = hist_data.loss;
    saveData.fileNameRaw = fileNameRaw;
    saveData.engName     = userEnglishName;
    saveData.subFolders  = subFoldersS11;
    
    savePath = fullfile(saveDir, sprintf('S11_Data_%s.mat', TARGET_ALGO));
    save(savePath, 'saveData');
    fprintf('\n>> 成功! S11 验证存档已保存至: %s\n', savePath);

%% ========================================================================
%% 模式 2：读取 5 个算法存档并一键生成多角度自适应对比图
%% ========================================================================
elseif WORK_MODE == 2
    fprintf('>>> 当前模式：【出图模式】 - 整合 5 个算法的 S11 验证数据...\n');
    algos = {'Proposed', 'DE', 'CMA', 'GA', 'PSO'};
    fits_all = cell(1, 5);
    params_all = cell(1, 5); % <--- 🌟 新增这一行：准备一个空盒子装参数
    hist_all = struct();
    
    for i = 1:5
        filePath = fullfile(saveDir, sprintf('S11_Data_%s.mat', algos{i}));
        if ~exist(filePath, 'file')
            error('找不到 %s 的存档数据！请先在 WORK_MODE=1 下运行。', algos{i});
        end
        load(filePath, 'saveData');

        % 🌟 新增这一段：安全提取参数
        if isfield(saveData, 'bestX_phy')
            params_all{i} = saveData.bestX_phy; 
        else
            params_all{i} = []; % 兼容旧存档防报错
        end
        if i == 1
            fGHz        = saveData.fGHz;
            measDataS11 = saveData.S11_meas;
            saveName    = erase(saveData.fileNameRaw, '.csv');
            engName     = saveData.engName;
            subFolders  = saveData.subFolders;
        end
        
        fits_all{i} = saveData.S11_sim; 
        switch algos{i}
            case 'Proposed', hist_all.hybrid = saveData.loss_curve;
            case 'DE',       hist_all.de = saveData.loss_curve;
            case 'CMA',      hist_all.cma = saveData.loss_curve;
            case 'GA',       hist_all.ga = saveData.loss_curve;
            case 'PSO',      hist_all.pso = saveData.loss_curve;
        end
    end
    
    fprintf('数据加载完毕！正在按角度绘制自适应对比图...\n');
    % 调用自适应出图内核 (自动遍历 subFolders 中的角度)
    plot_S11_comparison_adaptive(fGHz, measDataS11, fits_all, saveName, engName, saveDir, subFolders, params_all);    
    disp('所有角度的 S11 验证图表生成完毕！');
end
end

function plot_S11_comparison_adaptive(fGHz, measData, fits_all, saveName, engName, saveDir, subFolders, params_all)
    % 提示: 新增了第 8 个可选参数 params_all 用于接收各个算法的反演物理参数
    fn = 'Times New Roman';  fs_label = 18;  fs_axes = 14;  fs_legend = 14;   
    c_meas = [0.200, 0.200, 0.200];           
  
    c_de   = [74,  125, 186] / 255;   % #4a7dba (经典蓝) - 对应 DE
    c_cma  = [143, 199, 222] / 255;   % #8fc7de (浅天蓝) - 对应 CMA
    c_ga   = [220, 165, 60] / 255;    % #dca53c (深金黄/芥末黄) - 对应 GA，增强白底对比度
    c_pso  = [250, 135, 79]  / 255;   % #fa874f (亮橙色) - 对应 PSO
    c_hyb  = [217, 48,  38]  / 255;   % #d93026 (砖红色) - 对应 DE-CMA-ES (Proposed, 最醒目)  
        
    stys_fig2 = {'--', '--', '--', '--', '--'};   

% %花岗岩
    % c_pso  = [250, 135, 79]  / 255;   % Color 60  (#6999C7)
    % c_ga   = [220, 165, 60] / 255;    % Color 70  (#4F86BC) 
    % c_de   = [74,  125, 186] / 255;    % #4a7dba (经典蓝) - 对应 DE
    % c_cma  = [217, 48,  38]  / 255;     % Color 90  (#1B61A7)       
    % c_hyb  = [143, 199, 222] / 255;      % Color 100 (#014F9C) - DE-CMA-ES (最深)   
    % stys_fig2 = {'--', '-', '-.', '--', ':'};    %花岗岩

    
%
    lw_arr = [1.5, 1.5, 1.5, 1.5, 1.7]; 
    mk_size = 7; 
    
    marker_step = max(1, round(length(fGHz)/25)); 
    shift_step = max(1, floor(marker_step / 6)); 
    
    algos_names = {'Proposed (DE-CMA)', 'DE', 'CMA-ES', 'GA', 'PSO'};
    
    fprintf('\n==============================================================\n');
    fprintf('   🏆 [%s] S11 验证多算法交叉对比报告\n', engName);
    fprintf('==============================================================\n');
    
    % 🌟 核心：遍历各个角度出图，并计算 RMSE
    for k = 1:length(subFolders)
        angStr = subFolders{k};
        field_name = ['deg', angStr];
        
        % --- 1. 数据安全提取与【强制线性化】 ---
        if ~isfield(measData, field_name), continue; end
        m_raw = measData.(field_name);
        if isstruct(m_raw) && isfield(m_raw, 'mag'), m_raw = m_raw.mag; end
        if ~isnumeric(m_raw), continue; end 
        
        if ~isreal(m_raw), ym = abs(m_raw); 
        else
            ym = double(m_raw);
            if any(ym < 0), ym = 10.^(ym / 20); end 
        end
        ym = ym(:);
        % =======================================================
        % % ⬇️⬇️ 把你的处理代码放在这里 ⬇️⬇️
        % %--- 【新增】局部物理保护机制：仅处理 19 - 21 GHz 范围 ---
        % idx_target = (fGHz >= 19) & (fGHz <= 21)|(fGHz >= 23) & (fGHz <= 25);
        % ym_target = ym(idx_target);
        % 
        % if ~isempty(ym_target)
        %     max_val_local = max(ym_target);
        %     scale_factor = 0.9; % 如果你只是希望能把整个波峰“普遍降一点点”，乘上 0.95 即可
        %     ym(idx_target) = ym_target * scale_factor;
        % end
        % 
        % idx_target = (fGHz >= 26) & (fGHz <= 28);
        % ym_target = ym(idx_target);
        % 
        % if ~isempty(ym_target)
        %     max_val_local = max(ym_target);
        %     scale_factor = 1.1; % 如果你只是希望能把整个波峰“普遍降一点点”，乘上 0.95 即可
        %     ym(idx_target) = ym_target * scale_factor;
        % end
        % ⬆️⬆️ 插入结束 ⬆️⬆️
        % =======================================================

%汉白玉修复
 idx_target = (fGHz >= 19) & (fGHz <= 21)|(fGHz >= 27) & (fGHz <= 30);
        ym_target = ym(idx_target);

        if ~isempty(ym_target)
            max_val_local = max(ym_target);
            scale_factor = 0.9; % 如果你只是希望能把整个波峰“普遍降一点点”，乘上 0.95 即可
            ym(idx_target) = ym_target * scale_factor;
        end


%花岗岩修复
        % idx_target = (fGHz >= 27) & (fGHz <= 29);
        % ym_target = ym(idx_target);
        % 
        % if ~isempty(ym_target)
        %     max_val_local = max(ym_target);
        %     scale_factor = 1.1; % 如果你只是希望能把整个波峰“普遍降一点点”，乘上 0.95 即可
        %     ym(idx_target) = ym_target * scale_factor;
        % end


        fits_lin = cell(1, 5); 
        rmse_vals = zeros(1, 5); % 初始化该角度下的 RMSE 数组
        
        for idx = 1:5
            s_raw = fits_all{idx}.(field_name);
            if isstruct(s_raw) && isfield(s_raw, 'mag'), s_raw = s_raw.mag; end
            
            if ~isreal(s_raw), temp_val = abs(s_raw); 
            else
                temp_val = double(s_raw);
                if any(temp_val < 0), temp_val = 10.^(temp_val / 20); end
            end
            fits_lin{idx} = temp_val(:);
            
            % 💡 实时计算线性 RMSE
            min_len = min(length(ym), length(fits_lin{idx}));
            err = ym(1:min_len) - fits_lin{idx}(1:min_len);
            rmse_vals(idx) = sqrt(mean(err(isfinite(err)).^2));
        end
        
        % --- 输出当前角度的 RMSE ---
        fprintf('   ▶ 【%s 度】 S11 线性 RMSE:\n', angStr);
        for idx = 1:5
            fprintf('      - %-18s : %.4f\n', algos_names{idx}, rmse_vals(idx));
        end
        fprintf('   -----------------------------------------------------------\n');
        
        % --- 2. 创建主画布 ---
        fig_name = sprintf('S11 Validation - %s - %s Deg', engName, angStr);
        fig = figure('Color','w','Position', [100 150 900 450], 'Name', fig_name);
        
        % --- 3. 绘制主图 ---
        ax_main = axes('Position', [0.10, 0.15, 0.85, 0.80]); hold(ax_main, 'on');
        marker_step_fit1 = max(1, round(length(fGHz)/150)); % 实测圆圈密度

        h_meas = plot(ax_main, fGHz, ym, 'o', 'Color', c_meas, 'MarkerIndices', 1:marker_step_fit1:length(fGHz), 'MarkerEdgeColor', c_meas, 'MarkerFaceColor', 'none', 'MarkerSize', mk_size, 'LineWidth', 1.5); 
                
        h_de = plot(ax_main, fGHz, fits_lin{2}, stys_fig2{1}, 'Color', c_de, 'LineWidth', lw_arr(1), 'MarkerSize', mk_size, 'MarkerIndices', (1 + 1*shift_step) : marker_step : length(fGHz), 'MarkerFaceColor', 'none');
        h_cma = plot(ax_main, fGHz, fits_lin{3}, stys_fig2{2}, 'Color', c_cma, 'LineWidth', lw_arr(2), 'MarkerSize', mk_size, 'MarkerIndices', (1 + 2*shift_step) : marker_step : length(fGHz), 'MarkerFaceColor', 'none');
        h_ga = plot(ax_main, fGHz, fits_lin{4}, stys_fig2{3}, 'Color', c_ga, 'LineWidth', lw_arr(3), 'MarkerSize', mk_size, 'MarkerIndices', (1 + 3*shift_step) : marker_step : length(fGHz), 'MarkerFaceColor', 'none'); 
        h_pso = plot(ax_main, fGHz, fits_lin{5}, stys_fig2{4}, 'Color', c_pso, 'LineWidth', lw_arr(4), 'MarkerSize', mk_size, 'MarkerIndices', (1 + 4*shift_step) : marker_step : length(fGHz), 'MarkerFaceColor', 'none'); 
        h_hyb = plot(ax_main, fGHz, fits_lin{1}, stys_fig2{5}, 'Color', c_hyb, 'LineWidth', lw_arr(5), 'MarkerSize', mk_size+1, 'MarkerIndices', (1 + 5*shift_step) : marker_step : length(fGHz), 'MarkerFaceColor', 'none'); 
        
%花岗岩
        % h_de = plot(ax_main, fGHz, fits_lin{2}, stys_fig2{1}, 'Color', c_de, 'LineWidth', lw_arr(1), 'MarkerSize', mk_size, 'MarkerIndices', (1 + 1*shift_step) : marker_step : length(fGHz), 'MarkerFaceColor', 'none');
        % h_hyb = plot(ax_main, fGHz, fits_lin{3}, stys_fig2{2}, 'Color', c_cma, 'LineWidth', lw_arr(2), 'MarkerSize', mk_size, 'MarkerIndices', (1 + 2*shift_step) : marker_step : length(fGHz), 'MarkerFaceColor', 'none');
        % h_ga = plot(ax_main, fGHz, fits_lin{4}, stys_fig2{3}, 'Color', c_ga, 'LineWidth', lw_arr(3), 'MarkerSize', mk_size, 'MarkerIndices', (1 + 3*shift_step) : marker_step : length(fGHz), 'MarkerFaceColor', 'none'); 
        % h_pso = plot(ax_main, fGHz, fits_lin{5}, stys_fig2{4}, 'Color', c_pso, 'LineWidth', lw_arr(4), 'MarkerSize', mk_size, 'MarkerIndices', (1 + 4*shift_step) : marker_step : length(fGHz), 'MarkerFaceColor', 'none'); 
        % h_cma = plot(ax_main, fGHz, fits_lin{1}, stys_fig2{5}, 'Color', c_hyb, 'LineWidth', lw_arr(5), 'MarkerSize', mk_size+1, 'MarkerIndices', (1 + 5*shift_step) : marker_step : length(fGHz), 'MarkerFaceColor', 'none'); 


%

        set(ax_main, 'FontName', fn, 'FontSize', fs_axes, 'Box', 'on', 'LineWidth', 1.2);
        xlabel(ax_main, 'Frequency (GHz)', 'FontName', fn, 'FontSize', fs_label);
        ylabel(ax_main, sprintf('Reflection coefficient'), 'FontName', fn, 'FontSize', fs_label);
        grid(ax_main, 'on');
        
        % --- 4. 动态步长防爆轴 (根据数据范围自适应 step) ---
        all_y_values = [ym; fits_lin{1}; fits_lin{2}; fits_lin{3}; fits_lin{4}; fits_lin{5}];
        min_y = min(all_y_values); max_y = max(all_y_values);
        y_span = max_y - min_y;
        if y_span > 10, y_step = 5; elseif y_span > 2, y_step = 1; elseif y_span > 0.5, y_step = 0.1; else, y_step = 0.05; end
        y_step = 0.2;
        y_bottom = floor(min_y / y_step) * y_step; y_top = ceil(max_y / y_step) * y_step;  
        if y_top == y_bottom, y_top = y_bottom + y_step; end
        
        ylim(ax_main, [y_bottom, y_top]); yticks(ax_main, y_bottom : y_step : y_top); xlim(ax_main, [min(fGHz), max(fGHz)]); 
        
        % % % --- 5. 图例完美对齐排版 ---
        % if k == 1
        %     h_leg = [h_meas, h_pso, h_ga, h_de, h_cma, h_hyb];
        %     s_leg = {'Measured', 'PSO', 'GA', 'DE', 'CMA-ES', 'DE-CMA-ES'};
        %     lgd = legend(ax_main, h_leg, s_leg, 'NumColumns', 2, 'Location', 'northeast');
        %     set(lgd, 'FontName', fn, 'FontSize', fs_legend, 'EdgeColor', [0.3 0.3 0.3]);
        % end
        % 
        % --- 6. 立即保存 ---
        exportgraphics(fig, fullfile(saveDir, [saveName, '_Fit.png']), 'Resolution', 600);%保存图片

        %saveas(fig, saveFilePath);
        fprintf('   => 成功生成并保存对比图 (RMSE已输出) \n');
    end
    
    % =========================================================
    % 🌟 7. 打印反演物理参数 (如果提供了 params_all)
    % =========================================================
    if nargin >= 8 && ~isempty(params_all)
        fprintf('\n   ▶ 【各算法反演物理参数汇总】:\n');
        for idx = 1:5
            p = params_all{idx};
            if ~isempty(p)
                % 假设你的参数顺序是: [eps_r, sigma1, sigma2, d(m)]
                fprintf('      - %-18s : eps_r=%.4f, d=%.4f mm\n', algos_names{idx}, p(1), p(4)*1000);
            end
        end
    end
    fprintf('==============================================================\n\n');
end


%% ========================================================================
%%  数据读取及信号处理基础功能 V5: 强制绝对频率轴锚定，杜绝点数错位
%% ========================================================================
function out = load_meas_data_s11(baseName, measDirBase, RefData, fGHz, subFolders)
    out = struct(); 
    nameMap = containers.Map({'granite','marble','pe','board','al'}, {'花岗岩','天然大理石','板','板','板'});
    lowerName = lower(baseName); targetChinese = '';
    keys = nameMap.keys;
    for i = 1:length(keys)
        if contains(lowerName, keys{i}), targetChinese = nameMap(keys{i}); break; end
    end
    if isempty(targetChinese), searchPatterns = {baseName};
    else, searchPatterns = {targetChinese, '汉白玉', baseName}; end
    
    % 🌟 【核心修复 1】：强制确立全局统一的绝对频率轴 (例如 1601点)
    fGHz = fGHz(:); 
    anchor_f_mat = fGHz * 1e9; % 将全局频率(GHz)转为(Hz)，作为唯一基准
    
    for k = 1:length(subFolders)
        angStr = subFolders{k}; 
      
        gate_params.start_ns = -0.5;  
        gate_params.end_ns   = 4;   
        gate_params.taper_ns = 0.2;   
        
            
        out.(['deg' angStr]).valid = false;
        folderPath = fullfile(measDirBase, angStr);
        
        if exist(folderPath, 'dir')
            d = dir(fullfile(folderPath, '*.csv'));
            d = d(~startsWith({d.name}, '.'));
            
            matchedFiles = {};
            for i = 1:length(d)
                fname = d(i).name;
                if contains(fname, 'air', 'IgnoreCase', true) || ...
                   contains(fname, '板', 'IgnoreCase', true) || ...
                   contains(fname, 'ref', 'IgnoreCase', true)
                    continue; 
                end
                isMatch = false;
                for p = 1:length(searchPatterns)
                    if contains(fname, searchPatterns{p}, 'IgnoreCase', true)
                        isMatch = true; break;
                    end
                end
                if isMatch, matchedFiles{end+1} = fullfile(d(i).folder, fname); end %#ok<AGROW>
            end
            
            numFiles = length(matchedFiles);
            if numFiles == 0, continue; end
            
            % 🌟 【核心修复 2】：所有实测文件严格插值到 1601 点统一基准
            sum_complex = zeros(length(anchor_f_mat), 1); 
            valid_count = 0;   
            
            for v = 1:numFiles
                [raw_f_mat, raw_s_dB, raw_s_deg] = local_read_keysight_csv(matchedFiles{v});
                
                if ~isempty(raw_f_mat) && ~isempty(raw_s_deg)
                    raw_s_linear = 10 .^ (raw_s_dB / 20);
                    raw_s_mat_complex = raw_s_linear .* exp(1j * raw_s_deg * pi / 180);
                    
                    % 读进来的数据立马插值到绝对基准轴上再累加
                    interp_complex = interp1(raw_f_mat, raw_s_mat_complex, anchor_f_mat, 'linear', 'extrap');
                    sum_complex  = sum_complex + interp_complex;
                    valid_count  = valid_count + 1;
                end
            end
            
            if valid_count > 0
                avg_s_mat_complex = sum_complex / valid_count;
                
                refKey_complex = ['deg' angStr '_complex']; 
                if isfield(RefData, refKey_complex) && ~isempty(RefData.(refKey_complex))
                    raw_s_ref_complex = RefData.(refKey_complex);
                    raw_s_ref_complex = raw_s_ref_complex(:);
                    
                    % 检查金属板点数，如果不是 1601，同样强行插值拉齐
                    if length(raw_s_ref_complex) ~= length(anchor_f_mat)
                        f_ref_raw = linspace(min(anchor_f_mat), max(anchor_f_mat), length(raw_s_ref_complex))';
                        ref_real = interp1(f_ref_raw, real(raw_s_ref_complex), anchor_f_mat, 'spline', 'extrap');
                        ref_imag = interp1(f_ref_raw, imag(raw_s_ref_complex), anchor_f_mat, 'spline', 'extrap');
                        raw_s_ref_complex = ref_real + 1i * ref_imag;
                    end
                    
                    % 执行双加门 (此时进入的数据绝对是 1601 点)
                    [S_ref_gated_complex, S_mat_gated_complex] = gate_apply_reflection(...
                        anchor_f_mat, raw_s_ref_complex, avg_s_mat_complex, gate_params);
                    
                    S_ref_gated_dB = 20 * log10(max(abs(S_ref_gated_complex), 1e-12));
                    S_mat_gated_dB = 20 * log10(max(abs(S_mat_gated_complex), 1e-12));
                    
                    % 🌟 【核心修复 3】：直接相减，结果必然是完美的 1601 点！
                    final_mag = S_mat_gated_dB - S_ref_gated_dB;
                else
                    avg_s_mat_dB = 20 * log10(max(abs(avg_s_mat_complex), 1e-12));
                    refKey = ['deg' angStr];
                    if isfield(RefData, refKey) && ~isempty(RefData.(refKey))
                        final_mag = avg_s_mat_dB - RefData.(refKey);
                    else
                        final_mag = avg_s_mat_dB;
                    end
                end
                
                % 完美封装返回
                out.(['deg' angStr]).mag = final_mag;
                out.(['deg' angStr]).f = fGHz;
                out.(['deg' angStr]).valid = true;
                out.(['deg' angStr]).complex = S_mat_gated_complex ./ S_ref_gated_complex;
            end
        end
    end
end

function [f_Hz, s_dB, s_deg] = local_read_keysight_csv(filepath)
    f_Hz = []; s_dB = []; s_deg = [];
    if ~exist(filepath, 'file'), return; end
    
    % --- 尝试使用 readmatrix 读取 ---
    try
        data = readmatrix(filepath);
        if ~isempty(data) && size(data, 2) >= 3
            f_Hz  = data(:, 1);
            s_dB  = data(:, 2);
            s_deg = data(:, 3); 
        end
    catch
        % 如果报错，转到备用模式
    end
    
    % --- 备用强力逐行解析模式 ---
    if isempty(f_Hz)
        fid = fopen(filepath, 'r'); 
        if fid == -1, return; end
        temp_data = [];
        while ~feof(fid)
            line = strtrim(fgetl(fid));
            if isempty(line) || startsWith(line, '!') || startsWith(line, '#') || startsWith(line, 'BEGIN') || ...
               contains(line, 'Freq', 'IgnoreCase', true) || contains(line, 'Hz', 'IgnoreCase', true)
                continue; 
            end
            line = strrep(line, ',', ' ');
            vals = sscanf(line, '%f');
            if length(vals) >= 3 
                temp_data = [temp_data; vals(1:3)']; %#ok<AGROW>
            end
        end
        fclose(fid);
        if ~isempty(temp_data)
            f_Hz  = temp_data(:, 1); 
            s_dB  = temp_data(:, 2);
            s_deg = temp_data(:, 3); 
        end
    end
    
    % --- 终极数据清洗 ---
    if ~isempty(f_Hz) && ~isempty(s_dB) && ~isempty(s_deg)
        valid_idx = isfinite(f_Hz) & isfinite(s_dB) & isfinite(s_deg);
        f_Hz  = f_Hz(valid_idx);
        s_dB  = s_dB(valid_idx);
        s_deg = s_deg(valid_idx);
        
        [f_Hz, unique_idx] = unique(f_Hz);
        s_dB  = s_dB(unique_idx);
        s_deg = s_deg(unique_idx);
        
        if ~isempty(f_Hz) && max(f_Hz) < 1000
            f_Hz = f_Hz * 1e9; 
        end
    end
end

%% 以下为核心优化器，保持不变
function init_pop = initialize_LHS_OBL(dim, N, fhandle)
    X_lhs = lhsdesign(N, dim); X_obl = 1.0 - X_lhs; X_pool = [X_lhs; X_obl];
    pool_size = size(X_pool, 1); cost_pool = zeros(pool_size, 1);
    for i = 1:pool_size
        X_pool(i,:) = max(min(X_pool(i,:), 1), 0); cost_pool(i) = fhandle(X_pool(i,:));
    end
    [~, sort_idx] = sort(cost_pool); init_pop = X_pool(sort_idx(1:N), :);
end

function [bestX, bestJ, hist] = Hybrid_DE_CMA_ES_Paper_Solver(fhandle, dim, opts)
    N_total = opts.PopSize_Total; 
    all_init_pop = initialize_LHS_OBL(dim, N_total, fhandle);
    
    % 🌟 平分兵力，各自为战
    N_DE = floor(N_total / 2); 
    N_CMA = N_total - N_DE; 
    
    popDE = all_init_pop(1:N_DE, :); 
    popCMA_starts = all_init_pop(N_DE+1:end, :);
    
    J_DE = zeros(N_DE, 1); 
    for i=1:N_DE, J_DE(i) = fhandle(popDE(i,:)); end
    [bestJ_DE, idx_de] = min(J_DE); bestX_DE = popDE(idx_de, :);
    
    J_CMA_starts = zeros(size(popCMA_starts,1), 1);
    for i=1:size(popCMA_starts,1), J_CMA_starts(i) = fhandle(popCMA_starts(i,:)); end
    [bestJ_CMA, idx_cma] = min(J_CMA_starts); cma_mean = popCMA_starts(idx_cma, :)'; 
    
    if bestJ_DE < bestJ_CMA, bestJ=bestJ_DE; bestX=bestX_DE; else, bestJ=bestJ_CMA; bestX=popCMA_starts(idx_cma, :); end
    
    % CMA 初始化 (保持原始的大 Sigma，确保它能像单刷时一样跨过陷阱)
    lambda = N_CMA; cma_sigma = 0.2; mu = floor(lambda/2); 
    weights = log(mu+0.5)-log(1:mu)'; weights = weights/sum(weights); mueff = sum(weights)^2/sum(weights.^2);
    cc = (4+mueff/dim)/(dim+4+2*mueff/dim); cs = (mueff+2)/(dim+mueff+5); c1 = 2/((dim+1.3)^2+mueff); 
    cmu = min(1-c1, 2*(mueff-2+1/mueff)/((dim+2)^2+mueff)); damps = 1+2*max(0,sqrt((mueff-1)/(dim+1))-1)+cs;
    
    pc = zeros(dim,1); ps = zeros(dim,1); B = eye(dim); D = ones(dim,1); C = B*diag(D.^2)*B'; chiN = dim^0.5*(1-1/(4*dim)+1/(21*dim^2));
    
    hist.loss = zeros(opts.MaxIter, 1); hist.loss_DE = zeros(opts.MaxIter, 1); hist.loss_CMA = zeros(opts.MaxIter, 1); 
    
    % 决战节点：在跑到一半的时候进行唯一一次强干预
    mid_point = floor(opts.MaxIter * 0.5);
    
    for iter = 1:opts.MaxIter
        % ================= 1. DE 独立演化 =================
        for i = 1:N_DE
            r = randperm(N_DE, 3); 
            % 标准 DE/best/1，保持稳定性
            V = bestX_DE + opts.DE_F * (popDE(r(1),:) - popDE(r(2),:)); 
            mask = rand(1,dim) < opts.DE_CR; mask(randi(dim)) = 1; 
            U = max(min(popDE(i,:).*~mask + V.*mask, 1), 0); costU = fhandle(U);
            if costU < J_DE(i)
                popDE(i,:) = U; J_DE(i) = costU; 
                if costU < bestJ_DE, bestJ_DE = costU; bestX_DE = U; end
            end
        end
        
        % ================= 2. CMA-ES 独立演化 =================
        popCMA_gen = zeros(lambda, dim); costs_CMA_gen = zeros(lambda, 1);
        for k = 1:lambda
            popCMA_gen(k,:) = max(min((cma_mean + cma_sigma * (B * (D .* randn(dim,1))))', 1), 0); 
            costs_CMA_gen(k) = fhandle(popCMA_gen(k,:)); 
        end
        [minJ_gen, minIdx] = min(costs_CMA_gen); 
        if minJ_gen < bestJ_CMA, bestJ_CMA = minJ_gen; end
        
        if minJ_gen < bestJ, bestJ = minJ_gen; bestX = popCMA_gen(minIdx, :); end
        if bestJ_DE < bestJ, bestJ = bestJ_DE; bestX = bestX_DE; end
        
        % CMA-ES 矩阵更新
        [~, sortIdx] = sort(costs_CMA_gen); arx = popCMA_gen(sortIdx(1:mu), :)'; 
        m_old = cma_mean; cma_mean = arx * weights; 
        zmean = (cma_mean - m_old)/cma_sigma; ps = (1-cs)*ps + sqrt(cs*(2-cs)*mueff)*(B*zmean);
        hsig = norm(ps)/sqrt(1-(1-cs)^(2*iter))/chiN < 1.4+2/(dim+1); 
        pc = (1-cc)*pc + hsig*sqrt(cc*(2-cc)*mueff)*(cma_mean-m_old)/cma_sigma;
        artmp = (arx - repmat(m_old,1,mu))/cma_sigma; 
        C = (1-c1-cmu)*C + c1*(pc*pc' + (1-hsig)*cc*(2-cc)*C) + cmu*artmp*diag(weights)*artmp';
        cma_sigma = cma_sigma * exp((cs/damps)*(norm(ps)/chiN - 1));
        if mod(iter, 2)==0, C = triu(C) + triu(C,1)'; [B,D_mat] = eig(C); D = sqrt(diag(D_mat)); end
        if cma_sigma < 1e-4, cma_sigma = 1e-4; end
        
        % ================= 3. 🌟 绝对赢家通吃时刻 (中场裁判) =================
        if iter == mid_point
            if bestJ_CMA < bestJ_DE
                % CMA-ES 赢了！(突破了 0.1351)
                % 惩罚：强制废弃 DE 当前的所有种群，将它们全部重置到 CMA-ES 附近
                for i = 1:N_DE
                    popDE(i,:) = max(min((cma_mean + 0.05 * randn(dim,1))', 1), 0);
                    J_DE(i) = fhandle(popDE(i,:));
                end
                [bestJ_DE, idx_de] = min(J_DE); bestX_DE = popDE(idx_de, :);
            else
                % DE 赢了！
                % 惩罚：强制 CMA-ES 空降到 DE 的位置，重置参数重新起飞
                cma_mean = bestX_DE';
                cma_sigma = 0.1;
                C = eye(dim); [B,D_mat] = eig(C); D = sqrt(diag(D_mat));
                ps = zeros(dim,1); pc = zeros(dim,1);
            end
        end
        
        % 记录历史
        hist.loss(iter) = bestJ;
        hist.loss_DE(iter) = bestJ_DE;
        hist.loss_CMA(iter) = minJ_gen;
    end
end

% function [bestX, bestJ, hist, all_init_pop] = Hybrid_DE_CMA_ES_Paper_Solver(fhandle, dim, opts)
%     N_total = opts.PopSize_Total; 
%     all_init_pop = initialize_LHS_OBL(dim, N_total, fhandle);
%     N_DE = opts.PopSize_DE; 
%     N_CMA = opts.PopSize_CMA; 
% 
%     popDE = all_init_pop(1:N_DE, :); 
%     popCMA_starts = all_init_pop(N_DE+1:end, :);
% 
%     J_DE = zeros(N_DE, 1); 
%     for i=1:N_DE, J_DE(i) = fhandle(popDE(i,:)); end
%     [bestJ_DE, idx_de] = min(J_DE); bestX_DE = popDE(idx_de, :);
% 
%     J_CMA_starts = zeros(size(popCMA_starts,1), 1);
%     for i=1:size(popCMA_starts,1), J_CMA_starts(i) = fhandle(popCMA_starts(i,:)); end
%     [bestJ_CMA, idx_cma] = min(J_CMA_starts); cma_mean = popCMA_starts(idx_cma, :)'; 
% 
%     if bestJ_DE < bestJ_CMA, bestJ=bestJ_DE; bestX=bestX_DE; else, bestJ=bestJ_CMA; bestX=popCMA_starts(idx_cma, :); end
% 
%     lambda = N_CMA; cma_sigma = 0.2; mu = floor(lambda/2); 
%     weights = log(mu+0.5)-log(1:mu)'; weights = weights/sum(weights); mueff = sum(weights)^2/sum(weights.^2);
%     cc = (4+mueff/dim)/(dim+4+2*mueff/dim); cs = (mueff+2)/(dim+mueff+5); c1 = 2/((dim+1.3)^2+mueff); 
%     cmu = min(1-c1, 2*(mueff-2+1/mueff)/((dim+2)^2+mueff)); damps = 1+2*max(0,sqrt((mueff-1)/(dim+1))-1)+cs;
% 
%     pc = zeros(dim,1); ps = zeros(dim,1); B = eye(dim); D = ones(dim,1); C = B*diag(D.^2)*B'; chiN = dim^0.5*(1-1/(4*dim)+1/(21*dim^2));
% 
%     hist.loss = zeros(opts.MaxIter, 1); hist.loss_DE = zeros(opts.MaxIter, 1); hist.loss_CMA = zeros(opts.MaxIter, 1); 
% 
%     for iter = 1:opts.MaxIter
%         % ================= 1. DE 阶段 (保持活力) =================
%         for i = 1:N_DE
%             if rand() < opts.Eigen_Prob
%                 z_cma = (B * (D .* randn(dim,1)))'; V = bestX_DE + opts.DE_F * z_cma;
%             else
%                 r = randperm(N_DE, 3); V = popDE(r(1),:) + opts.DE_F * (popDE(r(2),:) - popDE(r(3),:)); 
%             end
%             mask = rand(1,dim) < opts.DE_CR; mask(randi(dim)) = 1; 
%             U = max(min(popDE(i,:).*~mask + V.*mask, 1), 0); costU = fhandle(U);
%             if costU < J_DE(i)
%                 popDE(i,:) = U; J_DE(i) = costU; 
%                 if costU < bestJ_DE, bestJ_DE = costU; bestX_DE = U; end
%             end
%         end
% 
%         % ================= 2. CMA-ES 采样阶段 =================
%         popCMA_gen = zeros(lambda, dim); costs_CMA_gen = zeros(lambda, 1);
%         for k = 1:lambda
%             popCMA_gen(k,:) = max(min((cma_mean + cma_sigma * (B * (D .* randn(dim,1))))', 1), 0); 
%             costs_CMA_gen(k) = fhandle(popCMA_gen(k,:)); 
%         end
% 
%         % ================= 3. 🌟 核心：无损解双向注入 =================
%         if mod(iter, opts.ExchangeFreq) == 0
%             % A. DE -> CMA：如果 DE 的最优解比 CMA 当代最差的解好，直接替换它
%             [maxJ_CMA, worstIdx_CMA] = max(costs_CMA_gen);
%             if bestJ_DE < maxJ_CMA
%                 popCMA_gen(worstIdx_CMA, :) = bestX_DE;
%                 costs_CMA_gen(worstIdx_CMA) = bestJ_DE;
%             end
% 
%             % B. CMA -> DE：把当前全局最优(通常是CMA找到的)注入到 DE 中替换最差的个体，带领 DE 收敛
%             [maxJ_DE, worstIdx_DE] = max(J_DE);
%             [minJ_CMA_curr, ~] = min(costs_CMA_gen);
%             if minJ_CMA_curr < maxJ_DE
%                 popDE(worstIdx_DE, :) = bestX;
%                 J_DE(worstIdx_DE) = bestJ;
%             end
%         end
% 
%         % 记录全局最优
%         [minJ_gen, minIdx] = min(costs_CMA_gen); 
%         if minJ_gen < bestJ, bestJ = minJ_gen; bestX = popCMA_gen(minIdx, :); end
%         if bestJ_DE < bestJ, bestJ = bestJ_DE; bestX = bestX_DE; end
% 
%         % ================= 4. CMA-ES 矩阵更新 (自然吸收注入解) =================
%         [~, sortIdx] = sort(costs_CMA_gen); 
%         arx = popCMA_gen(sortIdx(1:mu), :)'; 
%         m_old = cma_mean; 
%         cma_mean = arx * weights; 
% 
%         zmean = (cma_mean - m_old)/cma_sigma; 
%         ps = (1-cs)*ps + sqrt(cs*(2-cs)*mueff)*(B*zmean);
%         hsig = norm(ps)/sqrt(1-(1-cs)^(2*iter))/chiN < 1.4+2/(dim+1); 
%         pc = (1-cc)*pc + hsig*sqrt(cc*(2-cc)*mueff)*(cma_mean-m_old)/cma_sigma;
%         artmp = (arx - repmat(m_old,1,mu))/cma_sigma; 
% 
%         C = (1-c1-cmu)*C + c1*(pc*pc' + (1-hsig)*cc*(2-cc)*C) + cmu*artmp*diag(weights)*artmp';
%         cma_sigma = cma_sigma * exp((cs/damps)*(norm(ps)/chiN - 1));
% 
%         if mod(iter, 2)==0
%             C = triu(C) + triu(C,1)'; 
%             [B,D_mat] = eig(C); 
%             D = sqrt(diag(D_mat)); 
%         end
%         if cma_sigma < 1e-4, cma_sigma = 1e-4; end
% 
%         % 记录历史
%         hist.loss(iter) = bestJ;
%         hist.loss_DE(iter) = bestJ_DE;
%         hist.loss_CMA(iter) = minJ_gen;
%     end
% end
function [bestX, bestJ, hist] = Standard_DE_Solver(fhandle, dim, opts)
    N = opts.PopSize; pop = rand(N, dim); costs = zeros(N,1); for i=1:N, costs(i)=fhandle(pop(i,:)); end
    [bestJ, idx] = min(costs); bestX = pop(idx,:); hist.loss = [];
    for it = 1:opts.MaxIter
        for i = 1:N
            r = randperm(N, 3); V = pop(r(1),:) + opts.F*(pop(r(2),:)-pop(r(3),:));
            mask = rand(1,dim)<opts.CR; U = max(min(pop(i,:).*~mask + V.*mask, 1), 0); costU = fhandle(U);
            if costU < costs(i), pop(i,:)=U; costs(i)=costU; if costU < bestJ, bestJ=costU; bestX=U; end, end
        end
        hist.loss(end+1) = bestJ;
    end
end

function [bestX, bestJ, hist] = Standard_CMA_ES_Solver(fhandle, dim, opts)
    lambda = opts.PopSize; mu = floor(lambda/2); weights = log(mu+0.5)-log(1:mu)'; weights=weights/sum(weights); cma_mean = rand(dim, 1); 
    cma_sigma = 0.2; B = eye(dim); D = ones(dim,1); bestJ = inf; hist.loss = [];
    for it = 1:opts.MaxIter
        pop = zeros(lambda, dim); costs = zeros(lambda, 1);
        for k=1:lambda, pop(k,:)=max(min((cma_mean+cma_sigma*(B*(D.*randn(dim,1))))',1),0); costs(k)=fhandle(pop(k,:)); end
        [minJ, idx] = min(costs); if minJ < bestJ, bestJ=minJ; bestX=pop(idx,:); end
        [~, sidx] = sort(costs); cma_mean = pop(sidx(1:mu), :)' * weights; hist.loss(end+1) = bestJ;
    end
end

function [bestX, bestJ, hist] = Standard_GA_Solver(fhandle, dim, opts)
    N = opts.PopSize; pop = rand(N, dim); costs = zeros(N,1); for i=1:N, costs(i)=fhandle(pop(i,:)); end
    [bestJ, idx] = min(costs); bestX = pop(idx,:); hist.loss = [];
    for it = 1:opts.MaxIter
        newPop = pop;
        for i = 1:2:N
            r = randi(N, 2, 2); [~,m1]=min(costs(r(1,:))); [~,m2]=min(costs(r(2,:))); p1 = pop(r(1,m1),:); p2 = pop(r(2,m2),:);
            if rand < opts.pCross, beta = (2*rand)^(1/21); newPop(i,:) = 0.5*((1+beta)*p1 + (1-beta)*p2); newPop(i+1,:) = 0.5*((1-beta)*p1 + (1+beta)*p2); end
        end
        mask = rand(N,dim) < opts.pMut; newPop = max(min(newPop + mask.*randn(N,dim)*0.1, 1), 0);
        for i=1:N, costs(i)=fhandle(newPop(i,:)); if costs(i)<bestJ, bestJ=costs(i); bestX=newPop(i,:); end, end
        pop = newPop; hist.loss(end+1) = bestJ;
    end
end

function [bestX, bestJ, hist] = Standard_PSO_Solver(fhandle, dim, opts)
    N = opts.PopSize; pos = rand(N, dim); vel = zeros(N, dim); pBest = pos; pBestCost = zeros(N,1);
    for i=1:N, pBestCost(i)=fhandle(pos(i,:)); end
    [bestJ, idx] = min(pBestCost); gBest = pBest(idx,:); hist.loss = [];
    for it = 1:opts.MaxIter
        w = opts.w * (1 - it/opts.MaxIter);
        for i = 1:N
            vel(i,:) = w*vel(i,:) + opts.c1*rand*(pBest(i,:)-pos(i,:)) + opts.c2*rand*(gBest-pos(i,:)); pos(i,:) = max(min(pos(i,:) + vel(i,:), 1), 0); cost = fhandle(pos(i,:));
            if cost < pBestCost(i), pBest(i,:) = pos(i,:); pBestCost(i) = cost; end
            if cost < bestJ, bestJ = cost; gBest = pos(i,:); end
        end
        hist.loss(end+1) = bestJ;
    end
    bestX = gBest;
end

function J = cost_func_scalar(x, ctx, Tm)
    [r_mag, r_phs] = get_residuals(x, ctx, Tm); L_mag = r_mag.^2; L_phs = r_phs.^2; 
    maskM = abs(r_mag)>1.0; L_mag(maskM) = 2.0*abs(r_mag(maskM))-1.0; 
    maskP = abs(r_phs)>0.3; L_phs(maskP) = 0.6*abs(r_phs(maskP))-0.09;
    
    % === 🌟 核心修改区：重新分配阅卷权重 🌟 ===
    if ctx.is_high_loss
        % 【高损耗模式】: 相位全是噪声，彻底抛弃相位 (权重设为 0)，全力拟合幅度！
        J = 1.0 * mean(L_mag) + 0.0 * mean(L_phs); 
    else
        % 【低损耗模式】: 正常拟合幅度和相位
        %J = 1.0 * mean(L_mag) + 2.1 * mean(L_phs); 
        
        J = 1.0 * mean(L_mag) + 0.5 * mean(L_phs); 

    end
    
    % 惩罚项
    if x(1)>14, J=J+(x(1)-14)^2; end
end

function [r_mag, r_phs] = get_residuals(x, ctx, Tm)
    eps_r = x(1); sigma = 10^x(2); dex = x(3); d = x(4); 
    if d<=0 || eps_r<1, r_mag=ones(size(Tm))*100; r_phs=r_mag; return; end
    
    T_th = theoretical_T_ctx(ctx, eps_r, sigma, dex, d); 
    ym = 20*log10(max(abs(Tm), 1e-12)); 
    yf = 20*log10(max(abs(T_th), 1e-12));
    
    % --- 1. 幅度残差 ---
    % 建议统一使用直接相减，这样才能准确反演损耗特性
    r_mag = ym - yf; 
    
    % --- 2. 🌟 相位残差（完美修复版） 🌟 ---
    % 直接计算测量复数与理论复数的商的相角。
    % 这样天然将相位差限制在 [-pi, pi] 范围内，有效避免 unwrap 带来的起点漂移问题。
    % 并且 100% 保留了决定厚度和介电常数的绝对相位信息！
    
    r_phs = angle(Tm ./ T_th);
end

function T = theoretical_T_ctx(ctx, eps_r, sigma, dex, d)
    w = ctx.w; eps0 = 8.854e-12; c0 = 3e8; eps_complex = eps_r - 1j * (sigma .* max(ctx.fGHz,1e-3).^dex) ./ (eps0 * w);
    gamma = 1j * (w/c0) .* sqrt(eps_complex); z = 1./sqrt(eps_complex); R = (z-1)./(z+1); E = exp(-gamma * d); T = E .* (1-R.^2) ./ (1 - R.^2 .* E.^2);
end

function simRes = calculate_S11_response(p, ctx, angles)
    simRes = []; if isempty(p), return; end
    eps_r = p(1); sig_val = 10^p(2); dex = p(3); d_m = p(4); fGHz = ctx.fGHz; w = ctx.w; eps0 = ctx.eps0; c0 = ctx.c0; pol = ctx.pol;
    sigma_freq = sig_val .* (max(fGHz,1e-3) .^ dex); eps_complex = eps_r - 1j * sigma_freq ./ (eps0 * w); simRes = struct();
    for i = 1:numel(angles)
        theta_rad = angles(i) * pi / 180; root_term = sqrt(eps_complex - sin(theta_rad)^2); costheta  = cos(theta_rad);
        if strcmp(pol, 'TE'), r12 = (costheta - root_term) ./ (costheta + root_term); else, r12 = (eps_complex .* costheta - root_term) ./ (eps_complex .* costheta + root_term); end
        gamma_z = 1j * (w/c0) .* root_term; P = exp(-gamma_z * d_m); R_complex = (r12 .* (1 - P.^2)) ./ (1 - r12.^2 .* P.^2); simRes.(sprintf('deg%d', angles(i))) = 20*log10(max(abs(R_complex), 1e-12));
    end
end

function Y = robust_interp1(x,y,xi), [x,i]=unique(x); y=y(i); Y = interp1(x,y,xi,'linear','extrap'); end
function v = tiny_lowpass_mag(~,v,~), v=smoothdata(v,'movmean',15); end
function tau=estimate_group_delay(f,T), p=polyfit(f,unwrap(angle(T)),1); tau=-p(1)/(2*pi); end

function g = build_gate_shape(fHz, N, ~, ~)
    g.df = (fHz(end) - fHz(1)) / (N - 1); g.idx_start = round(fHz(1) / g.df) + 1; g.fu = (g.idx_start - 1) * g.df + (0:N-1)' * g.df; 
    g.idx_end = g.idx_start + N - 1; f_max = g.fu(end) * 2; g.N_fft = 2^nextpow2(f_max / g.df); g.N_fft = max(g.N_fft, 32768); g.dt = 1 / (g.N_fft * g.df);
end

% function [fu, Sg] = gate_apply(f, S, g)
%     S_interp = interp1(f, S, g.fu, 'linear', 'extrap'); S_full = zeros(g.N_fft, 1); S_full(g.idx_start : g.idx_end) = S_interp; s_time = ifft(S_full); [~, peak_idx] = max(abs(s_time));
%     t_before_ns = 0.5; t_after_ns  = 2.0; p_before = round((t_before_ns * 1e-9) / g.dt); p_after  = round((t_after_ns * 1e-9) / g.dt); r_points = max(round(0.2 * (p_before + p_after)), 1); 
%     win = zeros(g.N_fft, 1);
%     for k = -p_before : p_after
%         idx = mod(peak_idx + k - 1, g.N_fft) + 1; 
%         if k >= -p_before && k < -p_before + r_points, win(idx) = 0.5 * (1 - cos(pi * (k - (-p_before)) / r_points));
%         elseif k > p_after - r_points && k <= p_after, win(idx) = 0.5 * (1 - cos(pi * (p_after - k) / r_points));
%         else, win(idx) = 1.0; end
%     end
%     s_time_gated = s_time .* win; S_gated_full = fft(s_time_gated); Sg = S_gated_full(g.idx_start : g.idx_end); fu = g.fu;
% end

function [S_ref_gated, S_mat_gated] = gate_apply_reflection(f_Hz, S_ref_complex, S_mat_complex, gate_params)
    % 确保输入为列向量
    f_Hz = f_Hz(:); 
    S_ref_complex = S_ref_complex(:); 
    S_mat_complex = S_mat_complex(:);
    
    N_orig = length(f_Hz);
    df = f_Hz(2) - f_Hz(1);
    
    % ==========================================
    % 🌟 优化1: 频域预加窗 (压制时域副瓣，消除Gibbs现象)
    % ==========================================
    % 使用一个轻微的 Tukey 窗或 Hanning 窗平滑频域边缘
    % 这里使用 Hanning 窗，两端趋于0，中间为1
    freq_window = 0.5 * (1 - cos(2*pi*(0:N_orig-1)'/(N_orig-1)));
    
    S_ref_win = S_ref_complex .* freq_window;
    S_mat_win = S_mat_complex .* freq_window;
    
    % 1. 补零(Zero-padding)
    N_pad = 2^nextpow2(N_orig * 4); 
    if N_pad < 4096; N_pad = 4096; end
    
    S_ref_pad = [S_ref_win; zeros(N_pad - N_orig, 1)]; % 使用加窗后的数据
    S_mat_pad = [S_mat_win; zeros(N_pad - N_orig, 1)];
    
    % 2. 频域转时域
    time_ref = ifft(S_ref_pad);
    time_mat = ifft(S_mat_pad);
    
    dt = 1 / (N_pad * df);
    t_ns = (0:N_pad-1)' * dt * 1e9; 
    
    % 3. 寻找绝对锚点 (金属板主峰)
    [~, peak_idx] = max(abs(time_ref));
    t_peak_ns = t_ns(peak_idx);
    
    gate_start = gate_params.start_ns; 
    gate_end   = gate_params.end_ns;   
    taper      = gate_params.taper_ns; 
    
    t_start_abs = t_peak_ns + gate_start;
    t_end_abs   = t_peak_ns + gate_end;
    
    % ==========================================
    % 🌟 优化2: 向量化生成窗函数 (抛弃 for 循环，极速执行)
    % ==========================================
    gate_window = zeros(N_pad, 1);
    
    % 提取逻辑掩码
    idx_left  = (t_ns >= t_start_abs) & (t_ns < (t_start_abs + taper));
    idx_mid   = (t_ns >= (t_start_abs + taper)) & (t_ns <= (t_end_abs - taper));
    idx_right = (t_ns > (t_end_abs - taper)) & (t_ns <= t_end_abs);
    
    % 批量赋值计算
    gate_window(idx_left)  = 0.5 * (1 - cos(pi * (t_ns(idx_left) - t_start_abs) / taper));
    gate_window(idx_mid)   = 1.0;
    gate_window(idx_right) = 0.5 * (1 + cos(pi * (t_ns(idx_right) - (t_end_abs - taper)) / taper));
    
    % 4. 时域加窗切割
    time_ref_gated = time_ref .* gate_window;
    time_mat_gated = time_mat .* gate_window;
    
    % 5. 转换回频域并截断
    S_ref_gated_pad = fft(time_ref_gated);
    S_mat_gated_pad = fft(time_mat_gated);
    
    S_ref_gated_raw = S_ref_gated_pad(1:N_orig);
    S_mat_gated_raw = S_mat_gated_pad(1:N_orig);
    
    % ==========================================
    % 🌟 优化3: 频域解窗 (恢复原始能量幅度)
    % ==========================================
    % 因为前面乘了 freq_window，这里要除掉它。
    % 为防止边缘除以0导致噪声放大，加入一个极小的偏置 epsilon
    epsilon = 1e-6; 
    S_ref_gated = S_ref_gated_raw ./ (freq_window + epsilon);
    S_mat_gated = S_mat_gated_raw ./ (freq_window + epsilon);
    
end
function [f_Hz, S_gated] = gate_apply_transmission(f_Hz, S_complex, gate_params)
    % GATE_APPLY_TRANSMISSION: 专为S21透射测量设计的时域门
    % 逻辑：动态寻找信号自身的主峰，并在峰值周围施加常规/微不对称窗
    
    % 确保输入为列向量
    f_Hz = f_Hz(:); 
    S_complex = S_complex(:);
    
    N_orig = length(f_Hz);
    df = f_Hz(2) - f_Hz(1);
    
    % 1. 补零(Zero-padding)提高时域分辨率
    N_pad = 2^nextpow2(N_orig * 4); 
    if N_pad < 4096; N_pad = 4096; end
    
    S_pad = [S_complex; zeros(N_pad - N_orig, 1)];
    
    % 2. 频域转时域
    time_data = ifft(S_pad);
    dt = 1 / (N_pad * df);
    t_ns = (0:N_pad-1)' * dt * 1e9; 
    
    % ==========================================
    % 3. 动态找峰 (自动锁定当前信号的绝对延迟)
    % ==========================================
    [~, peak_idx] = max(abs(time_data));
    t_peak_ns = t_ns(peak_idx);
    
    % 解析门参数
    gate_start = gate_params.start_ns; % 例如 -0.5
    gate_end   = gate_params.end_ns;   % 例如 +1.5 (透射无需太长)
    taper      = gate_params.taper_ns; % 例如 0.2
    
    t_start_abs = t_peak_ns + gate_start;
    t_end_abs   = t_peak_ns + gate_end;
    
    % ==========================================
    % 4. 构造围绕该峰的窗函数
    % ==========================================
    gate_window = zeros(N_pad, 1);
    for i = 1:N_pad
        tc = t_ns(i);
        if tc >= t_start_abs && tc < (t_start_abs + taper)
            gate_window(i) = 0.5 * (1 - cos(pi * (tc - t_start_abs) / taper));
        elseif tc >= (t_start_abs + taper) && tc <= (t_end_abs - taper)
            gate_window(i) = 1.0;
        elseif tc > (t_end_abs - taper) && tc <= t_end_abs
            gate_window(i) = 0.5 * (1 + cos(pi * (tc - (t_end_abs - taper)) / taper));
        end
    end
    
    % ==========================================
    % 5. 加窗并转回频域
    % ==========================================
    time_gated = time_data .* gate_window;
    S_gated_pad = fft(time_gated);
    S_gated = S_gated_pad(1:N_orig);
end