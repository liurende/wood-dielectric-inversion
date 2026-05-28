function yq = robust_interp1(x, y, xq)
% 复数安全插值：清洗 NaN/Inf、排序去重，再对实部/虚部分别插值
    x  = double(x(:));
    y  = y(:);
    xq = double(xq(:));

    % 只保留有限值
    m = isfinite(x) & isfinite(real(y)) & isfinite(imag(y));
    x = x(m);  y = y(m);

    % 去重并按 x 排序
    if isempty(x)
        error('robust_interp1: 输入为空或全是非有限值');
    end
    [x, idx] = sort(x, 'ascend');
    y = y(idx);

    % 如果仍有重复 x，取第一次出现（或也可改成 mean 汇总）
    [xu, ia] = unique(x, 'stable'); 
    yu = y(ia);

    % 点数太少时兜底
    if numel(xu) < 2
        yq = repmat(yu(end), size(xq));
        return;
    end

    yr = interp1(xu, real(yu), xq, 'linear', 'extrap');
    yi = interp1(xu, imag(yu), xq, 'linear', 'extrap');
    yq = yr + 1j*yi;
end
