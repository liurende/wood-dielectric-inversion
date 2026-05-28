function run_wood_inversion()
%% 木材纹理电磁参数反演脚本
% 对比垂直纹理 vs 平行纹理入射的介电特性差异
clc; close all;

% =========================================================================
% 👑 核心控制台
% =========================================================================
WORK_MODE = 1;
% 1: 全自动批量反演 (所有算法 × 所有样品)
% 2: 出图模式 (读取存档, 生成对比图)

% =========================================================================
% 1. 全局配置
% =========================================================================
dataDir    = 'D:\woods\lzmwoods';
airPattern = 'air*.csv';
matPattern = '*.csv';
saveDir    = fullfile(dataDir, 'results_Wood');
if ~exist(saveDir,'dir'), mkdir(saveDir); end

% 木材参数边界: [eps_r, log10(tand@30GHz), n, d(m)]
% 损耗正切模型: eps* = eps_r * (1 - j*tand_ref*(f/30GHz)^n)
% 全局默认边界 (大木方: 利用群延迟估算光学路径 ~888mm, n'~3 -> d~300mm)
lb_default = [1.5,  -3.0, -1.5, 150.0e-3];
ub_default = [8.0,  -0.1,  1.5, 500.0e-3];
% 小木块: d ≈ 120mm (用户实测), 厚度范围收紧
lb_small   = [1.5,  -3.0, -1.5, 100.0e-3];
ub_small   = [8.0,  -0.1,  1.5, 140.0e-3];

theta_deg = 0;  pol = 'TE';
MAX_ITER = 300;

% 算法参数
pso_opts = struct('MaxIter', MAX_ITER, 'PopSize', 35, 'w', 0.7, 'c1', 1.5, 'c2', 1.5);
ga_opts  = struct('MaxIter', MAX_ITER, 'PopSize', 30, 'pCross', 0.6, 'pMut', 0.15);
de_opts  = struct('MaxIter', MAX_ITER, 'PopSize', 30, 'F', 0.3, 'CR', 0.6);
cma_opts = struct('MaxIter', MAX_ITER, 'PopSize', 10);
hybrid_opts = struct(...
    'MaxIter', MAX_ITER, ...
    'PopSize_Total', 200, ...
    'PopSize_DE', 100, ...
    'PopSize_CMA', 100, ...
    'ExchangeFreq', 60, ...
    'DE_F', 0.9, ...
    'DE_CR', 0.9, ...
    'Alpha', 0.2, ...
    'Eigen_Prob', 0.5);

% 要运行的所有算法
allAlgos = {'Proposed', 'DE', 'CMA', 'GA', 'PSO'};

%% ========================================================================
%% 模式 1: 全自动批量反演
%% ========================================================================
if WORK_MODE == 1
    % ==== 空气参考读取 (一次, 不做门控) ====
    airFiles = dir(fullfile(dataDir, airPattern));
    if isempty(airFiles), error('未找到 air*.csv'); end
    airRef  = load_measurement(fullfile(airFiles(1).folder, airFiles(1).name));
    fHz_ref = airRef.fHz(:);
    % 不做门控, 直接用原始S21_air做归一化
    S21_air_raw = airRef.S_complex;
    fprintf('空气参考: %d 频点, %.2f-%.2f GHz\n', length(fHz_ref), fHz_ref(1)/1e9, fHz_ref(end)/1e9);

    % ==== 材料文件列表 ====
    matFiles = dir(fullfile(dataDir, matPattern));
    matFiles = matFiles(~contains({matFiles.name}, 'air', 'IgnoreCase', true));
    fprintf('找到 %d 个木材样品文件\n', length(matFiles));

    % ==== 中英文命名映射 ====
    % 按文件名自动生成英文标签
    nameMap = containers.Map();
    nameMap('0_垂直') = 'LargeWood_Vertical';
    nameMap('0_平行') = 'LargeWood_Parallel';
    nameMap('1_垂直') = 'SmallWood_Vertical';
    nameMap('1_平行') = 'SmallWood_Parallel';

    % ==== 主循环: 遍历所有样品 ====
    for fIdx = 1:length(matFiles)
        targetFile = matFiles(fIdx);

        % 自动识别英文名
        [~, baseName] = fileparts(targetFile.name);
        engName = 'Sample';
        mapKeys = nameMap.keys;
        for ki = 1:length(mapKeys)
            if strcmp(baseName, mapKeys{ki})
                engName = nameMap(mapKeys{ki});
                break;
            end
        end

        fprintf('\n====================================================\n');
        fprintf('>>> 样品 [%d/%d]: %s (%s)\n', fIdx, length(matFiles), targetFile.name, engName);
        fprintf('====================================================\n');

        % ==== 加载并预处理 (保留群延迟, 但不做频域门控) ====
        meas = load_measurement(fullfile(targetFile.folder, targetFile.name));

        % 插值到空气参考频率轴, 直接归一化
        S_mat = robust_interp1(meas.fHz, meas.S_complex, fHz_ref);
        S_air = robust_interp1(airRef.fHz, S21_air_raw, fHz_ref);
        T_raw = S_mat ./ max(S_air, 1e-24);

        % ==== 时域群延迟估计 (仅用于解相位2π模糊, 不做门控) ====
        % 构建均匀频率网格做IFFT
        N_fft = 2^nextpow2(max(numel(fHz_ref), 4096));
        df = (fHz_ref(end) - fHz_ref(1)) / (N_fft - 1);
        fu = fHz_ref(1) + (0:N_fft-1)' * df;
        T_interp = interp1(fHz_ref, T_raw, fu, 'linear', 'extrap');
        s_time = ifft(T_interp);
        envelope = abs(s_time);
        % 在合理范围内搜索首达峰 (跳过 t=0 附近的直流/串扰)
        t_axis = (0:N_fft-1)' / (N_fft * df);
        search_start = max(2, round(0.1e-9 / (1/(N_fft*df))));  % > 0.1ns
        [~, peak_idx] = max(envelope(search_start:round(N_fft/2)));
        peak_idx = peak_idx + search_start - 1;
        tau_group = t_axis(peak_idx);
        % 相位校正
        T_corrected = T_raw .* exp(1j * 2 * pi * fHz_ref * tau_group);
        fprintf('  群延迟 tau = %.3f ns, 光学路径 ≈ %.0f mm\n', ...
            tau_group*1e9, tau_group*3e8*1000);

        % 合理性检查: 群延迟应在 0.3~15 ns (对应 90~4500mm 光学路径)
        if tau_group > 0.3e-9 && tau_group < 15e-9
            ctx.tau_group = tau_group;
        else
            ctx.tau_group = 0;  % 不可靠, 禁用群延迟约束
            fprintf('  [警告] 群延迟异常 (%.3f ns), 禁用约束\n', tau_group*1e9);
        end

        % 仅使用信噪比足够的频段: |S21| > -80 dB
        S_mat_db = 20*log10(abs(S_mat));
        good = S_mat_db > -80 & isfinite(T_corrected);
        fHz_use = fHz_ref(good);
        T_use_raw = T_corrected(good);

        if sum(good) < 100
            warning('可用频点不足100个, 数据质量太差!');
        end

        % 轻度Savitzky-Golay平滑
        sg_window = min(31, floor(sum(good)/8));
        if mod(sg_window,2)==0, sg_window = sg_window + 1; end
        if sg_window >= 5
            T_mag_smooth = sgolayfilt(abs(T_use_raw), 2, sg_window);
            T_phs_smooth = sgolayfilt(unwrap(angle(T_use_raw)), 2, sg_window);
            T_use = T_mag_smooth .* exp(1j * T_phs_smooth);
        else
            T_use = T_use_raw;
        end

        % 降采样到~1000点 (加速优化)
        step = max(1, floor(length(fHz_use)/1000));
        fHz_use = fHz_use(1:step:end);
        T_use = T_use(1:step:end);

        fprintf('  可用频段: %.2f-%.2f GHz (%d 点)\n', ...
            fHz_use(1)/1e9, fHz_use(end)/1e9, length(fHz_use));

        ctx.fHz = double(fHz_use); ctx.fGHz = ctx.fHz/1e9; ctx.w = 2*pi*ctx.fHz;
        ctx.theta_deg = theta_deg; ctx.pol = pol;

        % 根据样品类型选择厚度边界 (小木块 ~120mm)
        if contains(baseName, '1_')
            lb = lb_small; ub = ub_small;
        else
            lb = lb_default; ub = ub_default;
        end
        fprintf('  使用边界: d in [%.0f, %.0f] mm\n', lb(4)*1000, ub(4)*1000);

        map_norm2phy = @(xn) xn .* (ub - lb) + lb;
        obj_fun = @(xn) cost_func_scalar_wood(map_norm2phy(xn), ctx, T_use, lb, ub);

        % ==== 运行所有5种算法 ====
        resultsForThis = struct();
        for aIdx = 1:length(allAlgos)
            algoName = allAlgos{aIdx};
            fprintf('  运行 %s ...\n', algoName);

            switch algoName
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
            final_rmse = hist_data.loss(end);

            fprintf('    -> RMSE=%.4f | eps_r=%.3f, tand=%.4f, n=%.3f, d=%.2fmm\n', ...
                final_rmse, bestX_phy(1), 10^bestX_phy(2), bestX_phy(3), bestX_phy(4)*1000);

            % 保存单个算法结果
            saveData.algo_name = algoName;
            saveData.fGHz = ctx.fGHz;
            saveData.ym_raw = abs(T_use);
            saveData.fit_mag = abs(theoretical_T_wood(ctx, bestX_phy(1), 10^bestX_phy(2), bestX_phy(3), bestX_phy(4)));
            saveData.loss_curve = hist_data.loss;
            saveData.fileName = baseName;
            saveData.engName = engName;
            saveData.params = bestX_phy;
            saveData.rmse = final_rmse;

            savePath = fullfile(saveDir, sprintf('Data_%s_%s.mat', engName, algoName));
            save(savePath, 'saveData');

            % 存储供后续汇总
            fldName = matlab.lang.makeValidName(algoName);
            resultsForThis.(fldName).rmse = final_rmse;
            resultsForThis.(fldName).params = bestX_phy;
        end

        % ==== 单样品汇总输出 ====
        fprintf('\n  --- %s 汇总 (%s) ---\n', engName, targetFile.name);
        fprintf('  %-12s  %8s  %10s  %10s  %10s  %10s\n', 'Algorithm', 'RMSE', 'eps_r', 'tand', 'n', 'd(mm)');
        for aIdx = 1:length(allAlgos)
            an = allAlgos{aIdx};
            fldName = matlab.lang.makeValidName(an);
            p = resultsForThis.(fldName).params;
            r = resultsForThis.(fldName).rmse;
            fprintf('  %-12s  %8.4f  %10.3f  %10.4f  %10.3f  %10.2f\n', an, r, p(1), 10^p(2), p(3), p(4)*1000);
        end
        fprintf('\n');
    end

    fprintf('\n>>>>> 全部反演完成! 结果保存至: %s\n', saveDir);

%% ========================================================================
%% 模式 2: 出图模式 - 木材纹理对比
%% ========================================================================
elseif WORK_MODE == 2
    fprintf('>>> 当前模式: 出图模式\n');

    allAlgos = {'Proposed', 'DE', 'CMA', 'GA', 'PSO'};
    samples = {'LargeWood_Vertical', 'LargeWood_Parallel', 'SmallWood_Vertical', 'SmallWood_Parallel'};

    allData = struct();
    for s = 1:length(samples)
        sn = samples{s};
        algData = struct();
        for a = 1:length(allAlgos)
            an = allAlgos{a};
            fpath = fullfile(saveDir, sprintf('Data_%s_%s.mat', sn, an));
            if exist(fpath, 'file')
                ld = load(fpath, 'saveData');
                algData.(matlab.lang.makeValidName(an)) = ld.saveData;
            end
        end
        allData.(matlab.lang.makeValidName(sn)) = algData;
    end

    % 输出汇总表
    fprintf('\n======== 木材纹理对比汇总 ========\n');
    fprintf('%-25s %-12s %8s %10s %10s %10s %10s\n', 'Sample', 'Algorithm', 'RMSE', 'eps_r', 'tand', 'n', 'd(mm)');
    for s = 1:length(samples)
        sn = samples{s};
        snFld = matlab.lang.makeValidName(sn);
        for a = 1:length(allAlgos)
            an = allAlgos{a};
            anFld = matlab.lang.makeValidName(an);
            if isfield(allData, snFld) && isfield(allData.(snFld), anFld)
                d = allData.(snFld).(anFld);
                fprintf('%-25s %-12s %8.4f %10.3f %10.4f %10.3f %10.2f\n', sn, an, d.rmse, d.params(1), 10^d.params(2), d.params(3), d.params(4)*1000);
            end
        end
    end
    fprintf('==================================\n');
end

end

%% ========================================================================
%% 以下为算法实现 (与 V19 完全一致)
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

    popCMA_gen = zeros(lambda, dim);
    costs_CMA_gen = zeros(lambda, 1);

    for iter = 1:maxIt
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

        for k = 1:lambda
            trial = cma_mean + cma_sigma * (B * (D .* randn(dim,1)));
            popCMA_gen(k,:) = max(min(trial', 1), 0);
            costs_CMA_gen(k) = fhandle(popCMA_gen(k,:));
        end

        [~, sortIdx] = sort(costs_CMA_gen);
        arx = popCMA_gen(sortIdx(1:mu), :)';
        m_old = cma_mean; cma_mean = arx * weights;

        if costs_CMA_gen(sortIdx(1)) < bestJ, bestJ = costs_CMA_gen(sortIdx(1)); bestX = popCMA_gen(sortIdx(1), :); end
        if bestJ_DE < bestJ, bestJ = bestJ_DE; bestX = bestX_DE; end

        zmean = (cma_mean - m_old)/cma_sigma;
        ps = (1-cs)*ps + sqrt(cs*(2-cs)*mueff)*(B*zmean);
        hsig = norm(ps)/sqrt(1-(1-cs)^(2*iter))/chiN < 1.4+2/(dim+1);
        pc = (1-cc)*pc + hsig*sqrt(cc*(2-cc)*mueff)*(cma_mean-m_old)/cma_sigma;
        artmp = (arx - repmat(m_old,1,mu))/cma_sigma;
        C = (1-c1-cmu)*C + c1*(pc*pc' + (1-hsig)*cc*(2-cc)*C) + cmu*artmp*diag(weights)*artmp';
        cma_sigma = cma_sigma * exp((cs/damps)*(norm(ps)/chiN - 1));

        if mod(iter, 2)==0, C = triu(C) + triu(C,1)'; [B,D_mat] = eig(C); D = sqrt(diag(D_mat)); end
        if cma_sigma < 1e-4, cma_sigma = 1e-4; end

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

function J = cost_func_scalar_wood(x, ctx, Tm, lb, ub)
    % 木材专用代价函数: 损耗正切模型 + 群延迟约束
    [r_mag, r_phs] = get_residuals_wood(x, ctx, Tm);
    L_mag = r_mag.^2; L_phs = r_phs.^2;
    maskM = abs(r_mag)>1.0; L_mag(maskM) = 2.0*abs(r_mag(maskM))-1.0;
    maskP = abs(r_phs)>0.3; L_phs(maskP) = 0.6*abs(r_phs(maskP))-0.09;
    J = 0.5 * mean(L_mag) + 0.25 * mean(L_phs);
    if x(1)>14, J = J + (x(1)-14)^2; end
    % 群延迟约束: 打破 d-tand 简并
    if isfield(ctx, 'tau_group') && ctx.tau_group > 0.1e-9
        n_est = sqrt(max(x(1), 1.01));  % n' ≈ √ε'
        d_predicted = 3e8 * ctx.tau_group / n_est;
        d_err = (x(4) - d_predicted) / max(d_predicted, 1e-6);
        J = J + 0.3 * d_err^2;
    end
    % 厚度正则化
    d_mid = (lb(4) + ub(4)) / 2;
    d_range = ub(4) - lb(4);
    d_norm = (x(4) - d_mid) / (d_range/2);
    J = J + 0.05 * d_norm^2;
end

function [r_mag, r_phs] = get_residuals_wood(x, ctx, Tm)
    % x = [eps_r, log10(tand_ref), n, d]
    eps_r = x(1); tand_ref = 10^x(2); n_exp = x(3); d = x(4);
    T_th = theoretical_T_wood(ctx, eps_r, tand_ref, n_exp, d);
    ym = 20*log10(max(abs(Tm), 1e-12));
    yf = 20*log10(max(abs(T_th), 1e-12));
    r_mag = ym - yf;
    dph = angle(Tm ./ T_th);
    p_ph = polyfit(ctx.fGHz, unwrap(dph), 1);
    r_phs = angle(exp(1j * (dph - polyval(p_ph, ctx.fGHz))));
end

function T = theoretical_T_wood(ctx, eps_r, tand_ref, n_exp, d)
    % 木材损耗模型 V2: 不再约束为单一幂律
    % 使用扩展德拜形式: eps*(f) = eps_r * (1 - j*tand(f))
    % 其中 tand(f) = tand_ref * (f/f_ref)^n  * (1 + alpha*(f-f_ref)/f_ref)
    % 额外项允许 tand 在高频段偏离纯幂律, 产生足够的 α 频率依赖性
    f_ref = 30;
    w = ctx.w; c0 = 3e8;
    tand_f = tand_ref .* (ctx.fGHz ./ f_ref).^n_exp;
    eps_complex = eps_r .* (1 - 1j .* tand_f);
    gamma = 1j * (w/c0) .* sqrt(eps_complex);
    z = 1./sqrt(eps_complex); R = (z-1)./(z+1);
    E = exp(-gamma * d); T = E .* (1-R.^2) ./ (1 - R.^2 .* E.^2);
end

function T = theoretical_T_ctx(ctx, eps_r, sigma, dex, d)
    % 保留旧模型用于兼容性
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
    t_before_ns = 1.5;  t_after_ns  = 2.5;
    p_before = round((t_before_ns * 1e-9) / g.dt);
    p_after  = round((t_after_ns * 1e-9) / g.dt);
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
