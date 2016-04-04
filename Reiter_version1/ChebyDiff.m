function res = ChebyDiff(x, order, xmin, xmax)
    % Calculate analytic derivative for Chebyshev polynomial
    % from 0 to order - 1, counting #order values
    x = 2 .* (x - xmin) ./ (xmax - xmin) - 1;
    poly = zeros(order, order);
    poly(1, 1) = 1;
    poly(2, 2) = 1;
    for i = 3:order
        poly(i, :) = 2 * [0, poly(i - 1, 1:order - 1)] -  poly(i - 2, :);
    end
    for i = 1:order
        poly(i, :) = linspace(0, order - 1, order) .* poly(i, :);
        poly(i, :) = [poly(i, 2:end), 0];
    end
    x_poly = repmat(x, 1, order);
    x_poly(:, 1) = 1;
    x_poly = cumprod(x_poly, 2);
    res = x_poly * poly' .* 2 ./ (xmax - xmin);
end
