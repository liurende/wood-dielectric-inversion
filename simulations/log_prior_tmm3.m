function lp = log_prior_tmm3(x, lb, ub)
    % Log-prior for TMM 3-parameter model.
    % Uniform on eps_inf_E, eps_inf_L within [lb(1:2), ub(1:2)].
    % Log-uniform on f_relax within [lb(3), ub(3)]: p(x) ~ 1/x.

    if any(x < lb) || any(x > ub)
        lp = -inf;
        return;
    end
    lp = -log(x(3));  % log-uniform prior on f_relax
end
