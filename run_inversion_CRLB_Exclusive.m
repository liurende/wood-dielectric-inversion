function run_MonteCarlo_CRLB_Verification()
%% 蒙特卡洛 CRLB 验证脚本 - 专门用于生成 RMSE vs CRLB 2x2 顶级验证图
clc; close all;

% =========================================================================
% 0. 启动并行计算池 (检查与初始化)
% =========================================================================
fprintf('>>> 正在检查并分配并行计算资源 (Worker)...\n');
poolobj = gcp('nocreate'); % 检查是否已有并行池
if isempty(poolobj)
    poolobj = parpool; % 如果没有，则启动默认并行池
end

% 兼容不同 MATLAB 版本的提示信息获取
try
    num_workers = poolobj.NumWorkers;
    fprintf('>>> 成功激活并行池！当前可用核心数 (Workers): %d 个 🚀\n\n', num_workers);
catch
    % 如果当前版本的 MATLAB 无法直接读取 NumWorkers，则跳过详细打印
    fprintf('>>> 成功激活并行多核计算环境！🚀\n\n');
end

% =========================================================================
% 1. 全局配置与“快速测试”开关
% =========================================================================
% ⚠️ 强烈建议：首次运行保持 FAST_MODE = true 验证逻辑。出正式图时设为 false！
FAST_MODE = true; 
FAST_MODE = false;

if FAST_MODE
    fprintf('>>> 当前处于 [快速测试模式]：极少量蒙特卡洛和扫描点...\n');
    d_scan_mm = linspace(10, 50, 5); % 仅扫描 5 个厚度点
    N_MC = 3;                        % 每个点跑 3 次蒙特卡洛
    MAX_ITER = 100;                   % 算法极速迭代
    PopSize = 50;
else
    % fprintf('>>> 当前处于 [完整论文模式]：计算量巨大，请耐心等待...\n');
    % d_scan_mm = linspace(10, 50, 10); % 扫描 15 个厚度点 (如需更平滑可改回15或更高)
    % N_MC = 500;                      % 每个点跑 500 次蒙特卡洛
    % MAX_ITER = 300;                  % 算法满血迭代
    % PopSize = 200;
    
    
    fprintf('>>> 当前处于 [异常点重跑模式]：计算量巨大，请耐心等待...\n');
    % 假设你查到那个异常点是 23.33 mm
    d_scan_mm = [32.22, 36.67, 41.11, 45.56, 50.0]; 

    % 为了攻克这个难点，可以针对性地把这个点的参数调高
    N_MC = 400;      % 蒙特卡洛次数拉高
    MAX_ITER = 400;  % 迭代次数拉高
    PopSize = 300;   % 种群放大
end

% 噪声设置 (严格对应图标题)
sigma_noise = 2.0e-03; 
noise_variance = sigma_noise^2;

% 仿真系统参数
ctx.fGHz = linspace(18, 40, 1601)';  % 频率 18GHz ~ 40GHz (已对齐高频要求)
ctx.w = 2 * pi * ctx.fGHz * 1e9;
Nf = length(ctx.fGHz);
% 
% %--- 真实参数设定：Marble（大理石）---
% GT_eps_r = 7.749;
% GT_c_log = log10(2.37e-3);   
% GT_dex   = 1;
% 参数搜索边界
% lb = [ 4.0,  -5.0,  0.0,  1.0e-3]; 
% ub = [10.0,  -1.0,  3.0, 60.0e-3];

% --- 真实参数设定：Granite (花岗岩) ---
% GT_eps_r = 4.365;         % 介电常数 \epsilon_r
% GT_c_log     = log10(3.817e-5);      % 线性系数 c
% GT_dex     = 2.522;         % 指数/厚度参数 d
% % 针对 Granite (花岗岩) 优化的边界设定
% % 真实值参考：[4.365, -4.42(log10), 2.522, 厚度(通常10~50mm)]
% lb = [ 1.0,  -6.0,  0.0,  1.0e-3]; 
% ub = [10.0,  -2.0,  4.0, 60.0e-3];

% --- 真实参数设定：White marble (汉白玉) ---
GT_eps_r = 8.214;         % 介电常数 \epsilon_r
GT_c_log     = log10(2.895e-7);      % 线性系数 c
GT_dex     = 2.885;         % 指数/厚度参数 d

% 针对 White marble (汉白玉) 优化的边界设定
lb = [ 2.0,  -8.0,  1.0,  1.0e-3]; 
ub = [15.0,  -4.0,  5.0, 60.0e-3];

% DE-CMA-ES 核心参数
hybrid_opts = struct(...
    'MaxIter', MAX_ITER, 'PopSize_Total', PopSize*2, ...    
    'PopSize_DE', PopSize, 'PopSize_CMA', PopSize, ...       
    'ExchangeFreq', 20, 'DE_F', 0.8, 'DE_CR', 0.9, ...            
    'Alpha', 0.3, 'Eigen_Prob', 0.5);

% 初始化结果记录矩阵
num_pts = length(d_scan_mm);
RMSE_records = zeros(4, num_pts);
CRLB_records = zeros(4, num_pts);

% =========================================================================
% 2. 核心双层循环：厚度扫描 + 蒙特卡洛仿真 (带进度条与心跳)
% =========================================================================
fprintf('开始厚度扫描与蒙特卡洛验证...\n');

% --- 进度条与数据队列初始化 ---
total_tasks = num_pts * N_MC;       % 总需执行的蒙特卡洛次数
completed_tasks = 0;                % 已完成任务计数器
start_time = tic;                   % 记录开始时间
h_waitbar = waitbar(0, '正在初始化多核计算环境...', 'Name', '仿真进度');
D = parallel.pool.DataQueue;        % 创建并行数据队列
afterEach(D, @update_progress);     % 将队列与更新进度条函数绑定

% 计算控制台打印频率：如果次数很少就每次都打，如果很多就分批打
print_step = max(1, floor(N_MC / 10)); 

map_norm2phy = @(xn) xn .* (ub - lb) + lb;

for i = 1:num_pts
    GT_d = d_scan_mm(i) * 1e-3; % 真实厚度 (米)
    fprintf('\n[%d/%d] 正在分配 Sample Thickness = %.1f mm 的计算任务到并行池...\n', i, num_pts, GT_d * 1000);
    
    % 生成理想无噪声的 S21
    T_true = theoretical_T_ctx(ctx, GT_eps_r, 10^GT_c_log, GT_dex, GT_d);
    
    % 计算当前真实参数下的理论 CRLB (红线数据)
    GT_params_linear = [GT_eps_r, 10^GT_c_log, GT_dex, GT_d];
    CRLB_records(:, i) = compute_CRLB_theoretical(GT_params_linear, ctx, noise_variance);
    
    % 蒙特卡洛内层循环 (蓝线数据) - 开启 parfor 多核计算
    mc_results = zeros(4, N_MC);
    
    parfor mc = 1:N_MC
        % 注入复数高斯白噪声 (AWGN)
        noise = (randn(Nf, 1) + 1j * randn(Nf, 1)) * (sigma_noise / sqrt(2));
        T_noisy = T_true + noise;
        
        obj_fun = @(xn) compute_loss(map_norm2phy(xn), T_noisy, ctx);
        
        % 运行算法
        [best_X_norm, ~, ~] = Hybrid_DE_CMA_ES_Paper_Solver(obj_fun, 4, hybrid_opts);
        best_X_norm = real(best_X_norm);
        bestX_phy = map_norm2phy(best_X_norm);
        
        % 保存当前结果
        mc_results(:, mc) = [bestX_phy(1); 
                             10^bestX_phy(2); 
                             bestX_phy(3); 
                             bestX_phy(4)];
                         
        % --- 【并行心跳：控制台打印 Worker 信息】 ---
        if mod(mc, print_step) == 0 || mc == 1
            task = getCurrentTask();
            if ~isempty(task)
                fprintf('  -> [Worker %d 汇报] 厚度 %.1fmm 的第 %d/%d 次仿真已完成！\n', task.ID, GT_d*1000, mc, N_MC);
            end
        end
        
        % --- 【进度条信号：向主线程发送完成信号】 ---
        send(D, 1);
    end
    
    % 统计当前厚度下的 RMSE
    RMSE_records(1, i) = sqrt(mean((mc_results(1, :) - GT_eps_r).^2));
    RMSE_records(2, i) = sqrt(mean((mc_results(2, :) - 10^GT_c_log).^2));
    RMSE_records(3, i) = sqrt(mean((mc_results(3, :) - GT_dex).^2));
    RMSE_records(4, i) = sqrt(mean((mc_results(4, :) - GT_d).^2));
    
    fprintf('   eps_r RMSE: %.4e | CRLB: %.4e\n', RMSE_records(1, i), CRLB_records(1, i));
end

% 关闭进度条
if isvalid(h_waitbar)
    close(h_waitbar);
end

% =========================================================================
% 3. 完美复刻：2x2 高精度验证图绘制
% =========================================================================
fig = figure('Name', 'Monte Carlo CRLB Verification', 'Color', 'w', 'Position', [100, 100, 900, 700]);
sgtitle('天然大理石\_best\_gate.mat (Noise \sigma=2.0e-03) - High Precision Verification', ...
    'Interpreter', 'tex', 'FontWeight', 'bold', 'FontSize', 14, 'FontName', 'Microsoft YaHei');

labels = {'\epsilon_r', 'c', 'd', 'h'};
y_labels = {'RMSE / CRLB', 'RMSE / CRLB', 'RMSE / CRLB', 'RMSE / CRLB'};

for i = 1:4
    subplot(2, 2, i);
    hold on; grid on;
    
    p_crlb = plot(d_scan_mm, CRLB_records(i, :), 'k-', 'LineWidth', 2.0, 'DisplayName', 'CRLB');
    p_rmse = plot(d_scan_mm, RMSE_records(i, :), 'k--.', 'MarkerSize', 15, 'LineWidth', 1.5, 'DisplayName', 'RMSE');
    
    set(gca, 'YScale', 'log', 'FontName', 'Times New Roman', 'FontSize', 11, 'Box', 'on');
    title(labels{i}, 'FontWeight', 'bold', 'FontSize', 12);
    xlabel('Sample Thickness (mm)');
    ylabel(y_labels{i});
    xlim([10 50]);
    
    ylim_min = min([CRLB_records(i,:), RMSE_records(i,:)]) * 0.1;
    ylim_max = max([CRLB_records(i,:), RMSE_records(i,:)]) * 10;
    if ylim_min < ylim_max && isfinite(ylim_min) && isfinite(ylim_max)
        ylim([ylim_min, ylim_max]);
    end
end
legend(subplot(2,2,1), [p_rmse, p_crlb], 'Location', 'best');
disp('>>> 仿真验证与绘图完成！');

% =========================================================================
% 4. 保存数据结果 
% =========================================================================
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
filename = ['Verification_Result_', timestamp, '.mat'];

save_data = struct();
save_data.d_scan_mm    = d_scan_mm;
save_data.RMSE_records = RMSE_records;
save_data.CRLB_records = CRLB_records;
save_data.GT_params    = GT_params_linear;
save_data.sigma_noise  = sigma_noise;
save_data.ctx          = ctx;
save_data.date         = timestamp;

%save(filename, '-struct', 'save_data');

save('Patch_3222366741114556500.mat', '-struct', 'save_data');

fprintf('\n>>> 运行数据已成功保存至: %s\n', filename);

% =========================================================================
% 🎯 嵌套函数：专门用于更新进度条 (可访问主函数变量)
% =========================================================================
    function update_progress(~)
        completed_tasks = completed_tasks + 1;
        elapsed = toc(start_time); % 已用时间（秒）
        
        % 计算剩余时间并格式化
        rem_time = (elapsed / completed_tasks) * (total_tasks - completed_tasks);
        hours = floor(rem_time / 3600);
        mins = floor(mod(rem_time, 3600) / 60);
        secs = floor(mod(rem_time, 60));
        
        prog_pct = completed_tasks / total_tasks;
        
        if hours > 0
            msg = sprintf('总进度: %.1f%% | 预计剩余: %d小时 %d分 %d秒', prog_pct*100, hours, mins, secs);
        else
            msg = sprintf('总进度: %.1f%% | 预计剩余: %d分 %d秒', prog_pct*100, mins, secs);
        end
        
        % 安全更新进度条
        if isvalid(h_waitbar)
            waitbar(prog_pct, h_waitbar, msg);
        end
    end

end % 主函数结束

%% ========================================================================
%% 核心函数：已知真值与噪声下计算理论 CRLB (终极防线：使用伪逆 pinv)
%% ========================================================================
function std_bounds = compute_CRLB_theoretical(GT_linear, ctx, sigma2_noise)
    num_params = 4;
    Nf = length(ctx.fGHz); 
    J_real = zeros(Nf, num_params);
    J_imag = zeros(Nf, num_params);
    
    for k = 1:num_params
        delta = max(1e-6 * abs(GT_linear(k)), 1e-12); 
        x_plus = GT_linear;  x_plus(k) = x_plus(k) + delta;
        x_minus = GT_linear; x_minus(k) = x_minus(k) - delta;
        
        T_plus  = theoretical_T_ctx(ctx, x_plus(1), x_plus(2), x_plus(3), x_plus(4));
        T_minus = theoretical_T_ctx(ctx, x_minus(1), x_minus(2), x_minus(3), x_minus(4));
        
        dT = (T_plus(:) - T_minus(:)) / (2 * delta);
        J_real(:, k) = real(dT);
        J_imag(:, k) = imag(dT);
    end
    
    J_stacked = [J_real; J_imag];
    scale_factors = max(abs(GT_linear(:)), 1e-15);
    H = diag(scale_factors);
    J_scaled = J_stacked * H; 
    FIM_scaled = (2 / sigma2_noise) * (J_scaled' * J_scaled);
    Cov_scaled = pinv(FIM_scaled);
    Cov_Matrix = H * Cov_scaled * H;
    std_bounds = sqrt(abs(diag(Cov_Matrix))); 
end

%% ========================================================================
%% 物理模型
%% ========================================================================
function T = theoretical_T_ctx(ctx, eps_r, sigma, dex, d)
    w = ctx.w; eps0 = 8.854e-12; c0 = 3e8; 
    eps_complex = eps_r - 1j * (sigma .* max(ctx.fGHz,1e-3).^dex) ./ (eps0 * w); 
    gamma = 1j * (w/c0) .* sqrt(eps_complex); 
    z = 1./sqrt(eps_complex); R = (z-1)./(z+1); 
    E = exp(-gamma * d); T = E .* (1-R.^2) ./ (1 - R.^2 .* E.^2); 
end

%% ========================================================================
%% DE-CMA-ES 核心算法 (极致纯净版)
%% ========================================================================
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
    N_DE = opts.PopSize_DE; N_CMA = opts.PopSize_CMA;
    popDE = all_init_pop(1:N_DE, :); popCMA_starts = all_init_pop(N_DE+1:end, :);
    
    J_DE = zeros(N_DE, 1); for i=1:N_DE, J_DE(i) = fhandle(popDE(i,:)); end
    [bestJ_DE, idx_de] = min(J_DE); bestX_DE = popDE(idx_de, :);
    
    J_CMA_starts = zeros(size(popCMA_starts,1), 1);
    for i=1:size(popCMA_starts,1), J_CMA_starts(i) = fhandle(popCMA_starts(i,:)); end
    [bestJ_CMA, idx_cma] = min(J_CMA_starts); cma_mean = popCMA_starts(idx_cma, :)'; 
    
    if bestJ_DE < bestJ_CMA, bestJ = bestJ_DE; bestX = bestX_DE; 
    else, bestJ = bestJ_CMA; bestX = popCMA_starts(idx_cma, :); end
    
    lambda = N_CMA; cma_sigma = 0.2; mu = floor(lambda/2);
    weights = log(mu+0.5)-log(1:mu)'; weights=weights/sum(weights); mueff=sum(weights)^2/sum(weights.^2);
    cc = (4+mueff/dim)/(dim+4+2*mueff/dim); cs = (mueff+2)/(dim+mueff+5); c1 = 2/((dim+1.3)^2+mueff);
    cmu = min(1-c1, 2*(mueff-2+1/mueff)/((dim+2)^2+mueff)); damps = 1+2*max(0,sqrt((mueff-1)/(dim+1))-1)+cs;
    pc = zeros(dim,1); ps = zeros(dim,1); B = eye(dim); D = ones(dim,1); C = B*diag(D.^2)*B'; chiN=dim^0.5*(1-1/(4*dim)+1/(21*dim^2));
    
    hist.loss = []; maxIt = opts.MaxIter;
    k_exchange = opts.ExchangeFreq; alpha_inject = opts.Alpha; eigen_prob = opts.Eigen_Prob;
    popCMA_gen = zeros(lambda, dim); costs_CMA_gen = zeros(lambda, 1);
    
    for iter = 1:maxIt
        for i = 1:N_DE
            if rand() < eigen_prob, z_cma = (B * (D .* randn(dim,1)))'; V = bestX_DE + opts.DE_F * z_cma;
            else, r = randperm(N_DE, 3); V = popDE(r(1),:) + opts.DE_F * (popDE(r(2),:) - popDE(r(3),:)); end
            mask = rand(1,dim) < opts.DE_CR; mask(randi(dim)) = 1;
            U = max(min(popDE(i,:).*~mask + V.*mask, 1), 0);
            costU = fhandle(U);
            if costU < J_DE(i)
                popDE(i,:) = U; J_DE(i) = costU; 
                if costU < bestJ_DE, bestJ_DE = costU; bestX_DE = U; end
            end
        end
        for k = 1:lambda
            trial = cma_mean + cma_sigma * (B * (D .* randn(dim,1)));
            popCMA_gen(k,:) = max(min(trial', 1), 0); costs_CMA_gen(k) = fhandle(popCMA_gen(k,:));
        end
        [~, sortIdx] = sort(costs_CMA_gen); arx = popCMA_gen(sortIdx(1:mu), :)'; 
        m_old = cma_mean; cma_mean = arx * weights; 
        
        if costs_CMA_gen(sortIdx(1)) < bestJ, bestJ = costs_CMA_gen(sortIdx(1)); bestX = popCMA_gen(sortIdx(1), :); end
        if bestJ_DE < bestJ, bestJ = bestJ_DE; bestX = bestX_DE; end
        
        zmean = (cma_mean - m_old)/cma_sigma; ps = (1-cs)*ps + sqrt(cs*(2-cs)*mueff)*(B*zmean);
        hsig = norm(ps)/sqrt(1-(1-cs)^(2*iter))/chiN < 1.4+2/(dim+1); 
        pc = (1-cc)*pc + hsig*sqrt(cc*(2-cc)*mueff)*(cma_mean-m_old)/cma_sigma;
        artmp = (arx - repmat(m_old,1,mu))/cma_sigma; 
        C = (1-c1-cmu)*C + c1*(pc*pc' + (1-hsig)*cc*(2-cc)*C) + cmu*artmp*diag(weights)*artmp';
        cma_sigma = cma_sigma * exp((cs/damps)*(norm(ps)/chiN - 1));
        if mod(iter, 2)==0, C = triu(C) + triu(C,1)'; [B,D_mat] = eig(C); D = sqrt(diag(D_mat)); end
        if cma_sigma < 1e-9, cma_sigma = 1e-9; end        
        if mod(iter, k_exchange) == 0 && bestJ_DE < fhandle(cma_mean')
            cma_mean = (1 - alpha_inject) * cma_mean + alpha_inject * bestX_DE';
        end
        hist.loss(end+1) = bestJ;
    end
end

function loss = compute_loss(phy, T_noisy, ctx)
    T_est = theoretical_T_ctx(ctx, phy(1), 10^phy(2), phy(3), phy(4));
    loss = sum(abs(T_noisy - T_est).^2) / length(T_noisy);
end