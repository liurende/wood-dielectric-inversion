function utils_plot(fHz, T_meas, T_fit, matName)
% 单图：纵轴=透射系数幅值(线性)|T|，横轴=频率(GHz)
    fHz    = fHz(:);
    T_meas = T_meas(:);
    T_fit  = T_fit(:);

    if numel(fHz)~=numel(T_meas) || numel(fHz)~=numel(T_fit)
        error('utils_plot: 输入长度不一致');
    end

    fGHz = fHz/1e9;

    figure('Color','w','Position',[100 100 820 480]);
    plot(fGHz, abs(T_meas), 'b-',  'LineWidth',1.6); hold on;
    plot(fGHz, abs(T_fit),  'r--', 'LineWidth',1.8);
    xlabel('频率 (GHz)');
    ylabel('透射系数 |T|（线性）');
    title(sprintf('拟合结果（%s）', strrep(matName,'_','\_')));
    legend('实测','拟合','Location','best');
    grid on; box on;
end
