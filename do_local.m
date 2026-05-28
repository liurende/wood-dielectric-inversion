function [xBest, jBest] = do_local(x0, j0, fhandle, lb, ub)
% 独立文件，parfor 可安全调用
    if ~isscalar(j0) || ~isfinite(j0)
        j0 = fwrap(x0, fhandle);
    end
    xBest = x0; jBest = j0;
    try
        opts = optimoptions('fmincon','Display','none', ...
            'MaxFunctionEvaluations',2000,'StepTolerance',1e-10);
        [x1, j1] = fmincon(@(x) fwrap(x,fhandle), x0, [],[],[],[], lb, ub, [], opts);
        if isfinite(j1) && j1 < jBest, xBest=x1; jBest=j1; end
    catch
    end
    try
        opts2 = optimset('Display','off','MaxFunEvals',4000,'TolX',1e-10,'TolFun',1e-10);
        [x2, j2] = fminsearch(@(z) fwrap(project(z,lb,ub),fhandle), x0, opts2);
        x2 = project(x2,lb,ub);
        if isfinite(j2) && j2 < jBest, xBest=x2; jBest=j2; end
    catch
    end
end

function x = project(x, lb, ub), x = min(max(x, lb), ub); end
function J = fwrap(x, fhandle), J = fhandle(x); if ~isfinite(J), J = realmax; end, end
