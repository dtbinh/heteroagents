function res = Chebyshev_sym(x, order, xmin, xmax)
    % Calculate from 0 to order - 1, counting #order values
    % Input for column vector only
    % Output is a matrix mat * order
    x = 2 .* (x - xmin) ./ (xmax - xmin) - 1;
    res = sym(zeros(numel(x), order));
    res(:, 1) = 1;
    res(:, 2) = x;
    for i = 3:order
        res(:, i) = 2 .* x .* res(:, i - 1) - res(:, i - 2);
    end
end
