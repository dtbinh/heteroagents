function [w, x] = GaussLegendre(order)
    poly = zeros(order + 1, order + 1);
    poly(1, 1) = 1;
    poly(2, 2) = 1;
    for i = 3:(order + 1)
        poly(i, :) = (2 * i - 3) * [0, poly(i - 1, 1:order)] - (i - 2) * poly(i - 2, :);
        poly(i, :) = poly(i, :) / (i - 1);
    end
    x = roots(fliplr(poly(order + 1, :)));
    x = sort(x); % Put in ascending order
    poly_dif = linspace(0, order, order + 1) .* poly(order + 1, :);
    poly_dif = poly_dif(2:end)';
    x_mat = cumprod([ones(order, 1), repmat(x, 1, order - 1)], 2);
    w = 2.0 ./ ((1 - x .^ 2) .* (x_mat * poly_dif) .^ 2);
end
