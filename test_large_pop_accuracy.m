function verify_algorithm_performance_4param_PureL2_HighN()
%% 算法性能验证脚本 (4参数 + N=200 + 纯净L2 + 强制对数谱 + ETA)
% 
% 修复核心：
% 1. [关键] 绘图时对 RMSE 数据加 1e-20。防止 RMSE=0 时 semilogy 失效变回线性图。
% 2. [视觉] 强制打开 'YScale' = 'log' 和 'MinorGrid'，确保出现“下宽上窄”的网格。
clc; clear; 

%% ==== 1. 并行环境启动 ====
fprintf('正在初始化并行计算环境 (N_mc=200)...\n');
delete(gcp('nocreate')); 
try
    parpool('Threads'); 
catch
    parpool(); 
end

%% ==== 2. 全局参数配置 ====
dataDir    = 'C:\Users\LENOVO\Documents\MATLAB\demo2\石头'; 
resultDir  = fullfile(dataDir, 'results');
airPattern = 'air*.csv';

d_sweep_mm = linspace(10, 100, 20); 
n_points   = length(d_sweep_mm);
sigma      = 0.01;      
N_mc       = 5;        % 保持 200
rng(42);

popSize         = 150;   
maxIter         = 100;  
stagnationThres = 20;
hybrid_opts     = struct('enable',true,'useDE',true,'useGA',true);

lb = [1,  -18.0, 0.20, 1e-3];   
ub = [20,  -2.0, 10.20, 50e-3]; 
theta_deg = 0; pol = 'TE';

fmin_opts = optimoptions('fmincon', 'Display', 'off', ...
    'Algorithm', 'sqp', 'StepTolerance', 1e-12, 'OptimalityTolerance', 1e-12, ... 
    'MaxFunctionEvaluations', 1000);

%% ==== 3. 准备工作 ====
fprintf('正在读取基础数据...\n');
airFiles = dir(fullfile(dataDir, airPattern));
if isempty(airFiles), error('未找到 air*.csv。'); end
airBase = load_measurement(fullfile(airFiles(1).folder, airFiles(1).name));

fHz = airBase.fHz(:);
ctx.fHz = double(fHz); ctx.fGHz = double(fHz/1e9); 
ctx.w = double(2*pi*fHz); ctx.theta_deg = double(theta_deg); ctx.pol = pol;

matFiles = dir(fullfile(resultDir, '*_best_gate.mat'));
if isempty(matFiles), error('results 文件夹为空！'); end

%% ==== 4. 多材料循环验证 ====
for m_idx = 1:length(matFiles)
    targetFile = matFiles(m_idx);
    matName = targetFile.name;
    baseName = erase(matName, '_best_gate.mat');
    
    tmp = load(fullfile(targetFile.folder, matName), 'bestX');
    base_params_full = double(tmp.bestX); 
    
    ref_eps = base_params_full(1); 
    ref_pC = base_params_full(2); 
    ref_dex = base_params_full(3);
    
    fprintf('\n>>> [%d/%d] 启动材料: %s (eps=%.2f)\n', m_idx, length(matFiles), baseName, ref_eps);
    
    results_rmse = zeros(n_points, 4); 
    results_crlb = zeros(n_points, 4); 
    
    t_mat_start = tic; 
    
    for i = 1:n_points
        d_curr_m = d_sweep_mm(i) / 1000;
        theta_true_4param = [ref_eps; ref_pC; ref_dex; d_curr_m];
        
        % === A. CRLB 计算 ===
        dT_dtheta = finite_diff_grad_4param(ctx, theta_true_4param);
        J_fisher  = zeros(4, 4);
        noise_var = sigma^2;
        for r = 1:4
            for c = 1:4
                J_fisher(r,c) = (2/noise_var) * real(sum(conj(dT_dtheta(r,:)) .* dT_dtheta(c,:)));
            end
        end
        J_reg = J_fisher + 1e-15*eye(4);
        crlb_raw_local = sqrt(abs(diag(inv(J_reg)))).';
        val_c_true = 10^theta_true_4param(2);
        crlb_raw_local(2) = crlb_raw_local(2) * (val_c_true * log(10));
        results_crlb(i, :) = crlb_raw_local;
        
        % === B. Monte Carlo ===
        T_clean = theoretical_T_ctx(ctx, theta_true_4param(1), 10^theta_true_4param(2), theta_true_4param(3), theta_true_4param(4));
        sq_err_acc = zeros(1, 4); 
        
        par_pop = popSize; par_max = maxIter; par_stag = stagnationThres; 
        par_lb = lb; par_ub = ub; par_hyb = hybrid_opts; par_fmin = fmin_opts;
        
        parfor k = 1:N_mc
            noise = (randn(size(T_clean)) + 1j*randn(size(T_clean))) * (sigma / sqrt(2));
            T_meas = T_clean + noise;
            currObj = @(x) objective_4param_pure_L2(x, ctx, T_meas);
            
            woa_internal = struct('useLocal',true,'localEvery',20,'pbestFraction',0.35, ...
                                  'eliteTopK',5,'restartOnStall',true,'restartFrac',0.20, ...
                                  'hybrid', par_hyb, 'useParallel', false);
            [best_run_X, ~] = WOA_Improved_fast(par_pop, par_max, par_lb, par_ub, 4, currObj, par_stag, woa_internal);
            best_run_X = best_run_X(:).';
            
            try
                [x_polished, ~] = fmincon(currObj, double(best_run_X), [],[],[],[], par_lb, par_ub, [], par_fmin);
                final_X = x_polished(:).'; 
            catch
                final_X = best_run_X;
            end
            
            vec_est  = final_X;
            vec_true = theta_true_4param(:).';
            err_vec = vec_est - vec_true; 
            c_est = 10^vec_est(2); c_true = 10^vec_true(2);
            err_vec(2) = c_est - c_true;
            
            sq_err_acc = sq_err_acc + err_vec.^2;
        end
        results_rmse(i, :) = sqrt(sq_err_acc / N_mc);
        
        % ETA
        t_elapsed = toc(t_mat_start);
        avg_time_per_point = t_elapsed / i;
        points_remain = n_points - i;
        eta_sec = avg_time_per_point * points_remain;
        if mod(i,1)==0
            eta_min = floor(eta_sec/60); eta_s = mod(eta_sec, 60);
            fprintf('   -> 进度: %d/%d (d=%.1f mm) | RMSE(eps): %.2e | 剩余: %02d分%02.0f秒\n', ...
                i, n_points, d_curr_m*1000, results_rmse(i,1), eta_min, eta_s);
        end
    end
    
    plot_results_safe_4param_LOG(baseName, d_sweep_mm, results_crlb, results_rmse, ref_eps, sigma);
end
fprintf('\n>> 验证完成。\n');
end

%% ================= 子函数 =================
function J = objective_4param_pure_L2(x, ctx, Tm)
    if any(~isfinite(x)), J = 1e8; return; end
    eps_r = x(1); pC = x(2); dex = x(3); dth = x(4);
    T_th = theoretical_T_ctx(ctx, eps_r, 10.^pC, dex, dth);
    residual = Tm - T_th;
    J = sum(abs(residual).^2);
end

function dT_dtheta = finite_diff_grad_4param(ctx, theta_4)
    d_step = [1e-6; 1e-6; 1e-6; 1e-8]; 
    nParams = 4; dT_dtheta = zeros(nParams, numel(ctx.fHz));
    for k = 1:nParams
        th_p = theta_4; th_p(k) = th_p(k) + d_step(k);
        th_m = theta_4; th_m(k) = th_m(k) - d_step(k);
        Tp = theoretical_T_ctx(ctx, th_p(1), 10^th_p(2), th_p(3), th_p(4));
        Tm = theoretical_T_ctx(ctx, th_m(1), 10^th_m(2), th_m(3), th_m(4));
        dT_dtheta(k, :) = (Tp - Tm) ./ (2 * d_step(k)); 
    end
end

% === [核心修改] 强制对数绘图函数 ===
function plot_results_safe_4param_LOG(matName, d_mm, res_crlb, res_rmse, ref_eps, sigma)
    paramNames = {'\epsilon_r', 'c (Linear)', 'd_{exp}', 'd_{thick} (m)'};
    fName = sprintf('Log Scale Validation: %s', matName);
    
    figure('Color','w','Name', fName, 'Position', [100, 100, 1000, 700]);
    
    for p = 1:4
        subplot(2, 2, p); hold on;
        
        % 1. 获取数据
        y_crlb = res_crlb(:, p);
        y_rmse = res_rmse(:, p);
        
        % 2. [关键修复]：防止 0 值导致对数轴失效
        % 如果值为0，semilogy 画不出来，导致变成线性轴
        % 我们给所有数据垫一个极小值 (1e-20)
        y_crlb = max(y_crlb, 1e-20);
        y_rmse = max(y_rmse, 1e-20);
        
        % 3. 绘图 (semilogy)
        h1 = semilogy(d_mm, y_crlb, 'r--', 'LineWidth', 2);
        h2 = semilogy(d_mm, y_rmse, 'bo-', 'LineWidth', 1.5, 'MarkerFaceColor','b', 'MarkerSize', 6);
        
        % 4. [视觉强制] 设置 Y 轴为 Log 并打开次级网格
        set(gca, 'YScale', 'log'); 
        grid on; 
        grid minor; % <--- 这行代码产生“下宽上窄”的细格线
        
        % 5. 调整 Y 轴范围，确保能看到格线
        % 自动寻找合适的范围，防止下限变成 -Inf
        maxY = max([max(y_crlb), max(y_rmse)]) * 10;
        minY = min([min(y_crlb), min(y_rmse)]) / 10;
        ylim([minY, maxY]);

        xlabel('Thickness (mm)'); 
        ylabel(['RMSE (Log) - ' paramNames{p}]);
        title(['Param: ' paramNames{p}]);
        box on;
        xlim([min(d_mm), max(d_mm)]);
        
        try
            legend([h1 h2], {'Theoretical Limit (CRLB)', 'WOA+fmincon RMSE'}, 'Location', 'best'); 
        catch
        end
    end
    sgtitle(sprintf('%s (Log Scale, N=200): (\\epsilon_r=%.2f)', matName, ref_eps));
    drawnow;
end

function T = theoretical_T_ctx(ctx, eps_r, c, d_exp, d_thick)
    fHz=ctx.fHz; fGHz=ctx.fGHz; w=ctx.w; 
    eps0=8.854187817e-12; c0=299792458.0; 
    sigma = c .* max(fGHz,1e-9).^max(d_exp,0);
    eps_c = eps_r - 1j*sigma./(eps0*w);
    th = ctx.theta_deg*pi/180;
    root = sqrt(eps_c - (sin(th))^2);
    if strcmpi(ctx.pol,'TE')
        Rp = (cos(th) - root) ./ (cos(th) + root);
    else
        Rp = (eps_c.*cos(th) - root) ./ (eps_c.*cos(th) + root);
    end
    q = 2*pi*(fHz/c0) * d_thick .* root;
    num = (1 - Rp.^2) .* exp(-1j*q);
    den = 1 - (Rp.^2).*exp(-1j*2*q); den(abs(den)<1e-15)=1e-15;
    T = num ./ den;
end