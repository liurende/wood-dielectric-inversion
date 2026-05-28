function Xout = opposition_init(X, lb, ub)
% 把对立点加入再筛选前 N
    [N, dim] = size(X);
    X_opp = bsxfun(@plus, lb, bsxfun(@minus, ub, X));
    Xcat  = [X; X_opp];
    % 打乱返回，外部会根据适应度取前 N
    idx = randperm(size(Xcat,1));
    Xout = Xcat(idx,:);
end
