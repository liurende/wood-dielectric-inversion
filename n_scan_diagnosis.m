function n_scan_diagnosis()
%% 扫描 tand 频率指数 n, 寻找最优模型
c0 = 3e8;

ld = load('D:/woods/lzmwoods/results_Wood/Data_LargeWood_Vertical_Proposed.mat', 'saveData');
d = ld.saveData;
fGHz = d.fGHz(:);
w = 2*pi*fGHz*1e9;
T_meas = d.ym_raw(:);

lb3 = [1.2, -3, 5e-3];
ub3 = [8, -0.1, 500e-3];

n_vals = -1.5:0.1:0.5;
results = zeros(length(n_vals), 5);

for ni = 1:length(n_vals)
    n_fix = n_vals(ni);

    best_cost = inf;
    best_x = zeros(1,3);
    for restart = 1:20
        x0 = [lb3(1)+rand*(ub3(1)-lb3(1)), ...
              lb3(2)+rand*(ub3(2)-lb3(2)), ...
              lb3(3)+rand*(ub3(3)-lb3(3))];
        try
            opts = optimset('Display','off','MaxIter',500,'TolX',1e-8);
            [x, cost] = fmincon(@(p) cost_n(p, n_fix, fGHz, w, T_meas), ...
                x0, [],[],[],[], lb3, ub3, [], opts);
            if cost < best_cost
                best_cost = cost;
                best_x = x;
            end
        catch
        end
    end

    results(ni,:) = [n_fix, best_x, best_cost];
    fprintf('n=%.1f: eps_r=%.2f, tand=%.4f, d=%.0fmm, cost=%.4f\n', ...
        n_fix, best_x(1), 10^best_x(2), best_x(3)*1000, best_cost);
end

[~, idx] = min(results(:,5));
fprintf('\n最优: n=%.1f, eps_r=%.2f, tand=%.4f, d=%.0fmm, cost=%.4f\n', ...
    results(idx,1), results(idx,2), 10^results(idx,3), results(idx,4)*1000);

figure('Color','w');
subplot(1,2,1);
plot(results(:,1), results(:,5), 'bo-', 'LineWidth', 2);
xlabel('n (frequency exponent)'); ylabel('Cost'); grid on;
title('Cost vs tand exponent n');

subplot(1,2,2);
yyaxis left;
plot(results(:,1), results(:,2), 'b-o', 'LineWidth', 2);
ylabel('\epsilon_r');
yyaxis right;
plot(results(:,1), results(:,4)*1000, 'r-s', 'LineWidth', 2);
ylabel('d (mm)');
xlabel('n'); grid on;
title('Optimal \epsilon_r and d vs n');
legend('\epsilon_r', 'd', 'Location','best');

saveas(gcf, 'D:/woods/lzmwoods/results_Wood/n_scan.png');
fprintf('n扫描图已保存\n');
end

function J = cost_n(p, n_fix, fGHz, w, T_meas)
    eps_r = p(1); tand_ref = 10^p(2); d_th = p(3);
    f_ref = 30; c0 = 3e8;
    tand_f = tand_ref .* (fGHz ./ f_ref).^n_fix;
    eps_c = eps_r .* (1 - 1j .* tand_f);
    gamma = 1j * (w/c0) .* sqrt(eps_c);
    z = 1./sqrt(eps_c); R = (z-1)./(z+1);
    E = exp(-gamma * d_th); T_th = E .* (1-R.^2) ./ (1 - R.^2 .* E.^2);
    r_mag = 20*log10(max(abs(T_meas),1e-12)) - 20*log10(max(abs(T_th),1e-12));
    L_mag = r_mag.^2;
    maskM = abs(r_mag)>1.0; L_mag(maskM) = 2.0*abs(r_mag(maskM))-1.0;
    J = mean(L_mag);
end
