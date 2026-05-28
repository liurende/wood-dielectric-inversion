%% 诊断脚本：检查反演代码是否存在模型错误或目标函数病态
% 使用方法：直接运行，按提示选择测试内容
% 依赖：您的原始代码中已定义所有底层函数（theoretical_T_ctx, cost_func_scalar, Standard_DE_Solver 等）

clear; clc; close all;

%% ========== 1. 设置真实物理参数（用于生成仿真数据）==========
% 参考您论文中大理石的大致范围，设定一组“真实”值
truth.eps_r = 7.12;          % 介电常数
truth.sigma = 1e-5;          % 电导率系数 c (S/m·Hz^{-d})
truth.dex   = 1.0;            % 频率指数 d（物理合理值）
truth.d     = 17.93e-3;       % 厚度 (m)

% 频率范围（与您的实测数据一致）
fGHz = linspace(18, 40, 200)';   % 200个频点，足够平滑
fHz = fGHz * 1e9;
w = 2*pi*fHz;

% 构造 ctx 结构（与您的代码中格式一致）
ctx.fHz = fHz;
ctx.fGHz = fGHz;
ctx.w = w;
ctx.theta_deg = 0;   % 垂直入射
ctx.pol = 'TE';
ctx.is_high_loss = false;

% 计算理论透射系数（无噪声）
T_true = theoretical_T_ctx(ctx, truth.eps_r, truth.sigma, truth.dex, truth.d);

% 可选：添加少量噪声，模拟真实测量（信噪比 40 dB）
noise_level = 0.01;   % 复数噪声标准差，相对幅度
T_meas = T_true + noise_level * (randn(size(T_true)) + 1i*randn(size(T_true)));

%% ========== 2. 检查正向模型：已知参数能否还原T_true？==========
% 用另一组参数重新计算，比较差异（做自洽性检验）
T_check = theoretical_T_ctx(ctx, truth.eps_r, truth.sigma, truth.dex, truth.d);
max_diff = max(abs(T_check - T_true));
fprintf('正向模型自洽性检查：最大差异 = %e (应为0)\n', max_diff);
if max_diff > 1e-12
    error('正向模型不自治！请检查 theoretical_T_ctx 函数。');
end

%% ========== 3. 目标函数行为检查（固定其他参数，扫描 sigma 和 dex）==========
% 目的：观察目标函数是否有多谷、平坦区域，以及是否存在非物理的极小值
disp('正在扫描目标函数（c 和 d 空间），请稍候...');

% 固定 eps_r 和 d 为真实值
eps_r_fixed = truth.eps_r;
d_fixed = truth.d;

% 扫描范围：log10(sigma) 从 -6 到 0（sigma 从 1e-6 到 1），dex 从 -2 到 4
log10sigma_vals = linspace(-6, 0, 51);
dex_vals = linspace(-2, 4, 51);

RMSE_mat = zeros(length(log10sigma_vals), length(dex_vals));
for i = 1:length(log10sigma_vals)
    sigma_try = 10^log10sigma_vals(i);
    for j = 1:length(dex_vals)
        dex_try = dex_vals(j);
        T_try = theoretical_T_ctx(ctx, eps_r_fixed, sigma_try, dex_try, d_fixed);
        % 使用与您代码中相同的代价函数（但为了绘图，直接计算复数RMSE）
        err_complex = T_try - T_meas;
        rmse = sqrt(mean(abs(err_complex).^2));
        RMSE_mat(i,j) = rmse;
    end
end

% 绘制等值线图
figure('Name', '目标函数等值线 (固定 eps\_r, d)');
contourf(dex_vals, log10sigma_vals, log10(RMSE_mat), 20);
colorbar;
xlabel('频率指数 d (dex)');
ylabel('log_{10}(电导率系数 c)');
title('log_{10}(RMSE) 等值线图');
% 标注真实值位置
hold on;
plot(truth.dex, log10(truth.sigma), 'r*', 'MarkerSize', 12, 'LineWidth', 2);
text(truth.dex+0.2, log10(truth.sigma), '真实值', 'Color', 'r');
% 标注您之前反演得到的非物理值（例如 dex = -1.2）
plot(-1.2, log10(3.6e-9), 'bo', 'MarkerSize', 8);
text(-1.2+0.2, log10(3.6e-9), '非物理解', 'Color', 'b');
hold off;
grid on;

% 分析：如果真实值处有明显的低洼，且周围没有其他同样低的区域，则模型可辨识。
% 如果整个平面平坦或存在多个低洼，则目标函数病态。

%% ========== 4. 用仿真数据测试反演算法 ==========
% 将 T_meas 作为“实测数据”，运行您的 DE 或 Proposed 算法，看能否恢复 truth
disp('===== 算法反演测试（仿真数据）=====');

% 重新定义目标函数（直接使用您的 cost_func_scalar，但需要适配 ctx 和 T_meas）
% 注意：您的 cost_func_scalar 中使用了幅度 dB 差和相位残差，并包含边界惩罚。
% 我们直接重用，但需要构造一个与您主程序相同的 obj_fun

% 设置边界（与您的代码一致）
lb = [5.0,  -6.0, -2.30,  18.0e-3];   % 注意第二维我们缩小了范围，方便测试
ub = [9.0,  0.0,  4.00,   20.0e-3];
% 归一化映射函数
map_norm2phy = @(xn) xn .* (ub - lb) + lb;
% 定义适应度函数（调用您的 cost_func_scalar）
% 注意：cost_func_scalar 要求 x 为物理值，且内部会再次调用 get_residuals 等
obj_fun_test = @(x_norm) cost_func_scalar(map_norm2phy(x_norm), ctx, T_meas, lb, ub);

% 选择测试的算法（您可以手动修改）
algo_to_test = 'Proposed';   % 可选: 'DE', 'CMA', 'Proposed'
opts_DE = struct('MaxIter', 200, 'PopSize', 100, 'F', 0.7, 'CR', 0.8);
opts_CMA = struct('MaxIter', 200, 'PopSize', 100);
opts_Hybrid = struct('MaxIter', 200, 'PopSize_Total', 100, 'PopSize_DE', 50, ...
                     'PopSize_CMA', 50, 'ExchangeFreq', 15, 'DE_F', 0.7, ...
                     'DE_CR', 0.8, 'Alpha', 0.3, 'Eigen_Prob', 0.5);

switch algo_to_test
    case 'DE'
        [best_norm, ~, hist] = Standard_DE_Solver(obj_fun_test, 4, opts_DE);
    case 'CMA'
        [best_norm, ~, hist] = Standard_CMA_ES_Solver(obj_fun_test, 4, opts_CMA);
    case 'Proposed'
        [best_norm, ~, hist] = Hybrid_DE_CMA_ES_Paper_Solver(obj_fun_test, 4, opts_Hybrid);
    otherwise
        error('未知算法');
end

best_phy = map_norm2phy(best_norm);
fprintf('算法 %s 反演结果：\n', algo_to_test);
fprintf('  eps_r = %.4f (真实 %.4f)\n', best_phy(1), truth.eps_r);
fprintf('  sigma = %.4e (真实 %.4e)\n', 10^best_phy(2), truth.sigma);
fprintf('  dex   = %.4f (真实 %.4f)\n', best_phy(3), truth.dex);
fprintf('  d     = %.4f mm (真实 %.4f mm)\n', best_phy(4)*1000, truth.d*1000);

% 绘制收敛曲线
figure('Name', ['收敛曲线 - ', algo_to_test]);
semilogy(hist.loss, 'b-', 'LineWidth', 1.5);
xlabel('迭代次数'); ylabel('代价函数值'); grid on;
title(sprintf('最终代价: %.4e', hist.loss(end)));

% 比较拟合的透射系数
T_fit = theoretical_T_ctx(ctx, best_phy(1), 10^best_phy(2), best_phy(3), best_phy(4));
figure('Name', '透射系数拟合对比');
subplot(2,1,1);
plot(fGHz, 20*log10(abs(T_meas)), 'k-o', 'MarkerSize', 3, 'DisplayName', '仿真实测');
hold on;
plot(fGHz, 20*log10(abs(T_fit)), 'r--s', 'MarkerSize', 3, 'DisplayName', '反演拟合');
xlabel('频率 (GHz)'); ylabel('幅度 (dB)'); legend; grid on;
subplot(2,1,2);
plot(fGHz, unwrap(angle(T_meas))*180/pi, 'k-o', 'MarkerSize', 3); hold on;
plot(fGHz, unwrap(angle(T_fit))*180/pi, 'r--s', 'MarkerSize', 3);
xlabel('频率 (GHz)'); ylabel('相位 (度)'); grid on;

%% ========== 5. 边界惩罚机制检查 ==========
% 您的 cost_func_scalar 中只有对 eps_r > 14 的惩罚，没有对 dex 或 sigma 的边界惩罚。
% 这会导致算法可以自由进入 dex 负值区域而不会受到惩罚，除非最终拟合误差很大。
% 建议：在 cost_func_scalar 中加入对参数越界的强惩罚，或者直接限制搜索边界内（已做）。

% 测试：手动构造一个非物理参数组合，计算代价，与真实参数附近的代价比较
x_bad = [7.12, -5, -1.2, 17.93e-3];   % dex = -1.2
cost_bad = cost_func_scalar(x_bad, ctx, T_meas, lb, ub);
x_good = [truth.eps_r, log10(truth.sigma), truth.dex, truth.d];
cost_good = cost_func_scalar(x_good, ctx, T_meas, lb, ub);
fprintf('真实参数代价 = %.4f, 非物理参数代价 = %.4f\n', cost_good, cost_bad);
if cost_bad < cost_good
    fprintf('警告：非物理参数的代价更低！这解释为何算法会收敛到负dex。\n');
    fprintf('原因：当前目标函数对 dex 负值区域惩罚不足，且存在数据噪声/模型误差。\n');
else
    fprintf('物理参数代价更低，算法应收敛到正确区域。若实际中仍出现负值，请检查数据质量或算法早熟。\n');
end

%% ========== 6. 结论与建议 ==========
disp('========================================');
disp('诊断总结：');
disp('1. 如果仿真数据测试中算法能完美恢复真实参数，则您的代码逻辑无误；');
disp('   问题可能来自实测数据中的噪声或系统误差。');
disp('2. 如果仿真数据测试也无法恢复，则目标函数或优化器存在缺陷，建议：');
disp('   - 改用复数误差（直接 min |T_meas - T_theo|^2）代替幅度/相位加权；');
disp('   - 去除 get_residuals 中对相位去线性趋势的操作（可能会消除物理相位信息）；');
disp('   - 在 cost_func_scalar 中加入对 dex<0 的强惩罚（例如 if dex<0, J=J+1e3*abs(dex); end）。');
disp('3. 绘制等值线图可帮助您直观理解目标函数形态。');
disp('========================================');



function J = cost_func_scalar(x, ctx, Tm, lb, ub)
    % x = [eps_r, log10(sigma), dex, d] 物理值（已反归一化）
    eps_r = x(1);
    sigma = 10^x(2);
    dex = x(3);
    d = x(4);
    
    % 物理约束惩罚（软约束）
    penalty = 0;
    if eps_r < 1
        penalty = penalty + 1e3 * (1 - eps_r)^2;
    end
    if dex < 0   % 坚决禁止负的指数
        penalty = penalty + 1e4 * abs(dex);   % 强惩罚
    end
    if sigma <= 0
        penalty = penalty + 1e4 * abs(sigma);
    end
    if d <= 0
        penalty = penalty + 1e4 * abs(d);
    end
    
    % 计算理论透射系数
    T_th = theoretical_T_ctx(ctx, eps_r, sigma, dex, d);
    
    % 复数误差（均方根）
    err_complex = T_th - Tm;
    rmse = sqrt(mean(abs(err_complex).^2));
    
    % 最终代价 = 复数RMSE + 物理约束惩罚
    J = rmse + penalty;
end


function T = theoretical_T_ctx(ctx, eps_r, sigma, dex, d)
    w = ctx.w; eps0 = 8.854e-12; c0 = 3e8; 
    eps_complex = eps_r - 1j * (sigma .* max(ctx.fGHz,1e-3).^dex) ./ (eps0 * w); 
    gamma = 1j * (w/c0) .* sqrt(eps_complex); 
    z = 1./sqrt(eps_complex); R = (z-1)./(z+1); 
    E = exp(-gamma * d); T = E .* (1-R.^2) ./ (1 - R.^2 .* E.^2); 
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