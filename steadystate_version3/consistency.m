function [res, jac] = consistency(m, g, s)
    res = mean(exp(m * g'));
    if nargout > 1
        jac = mean(m .* repmat(exp(m * g'), 1, s), 1);
    end
end
