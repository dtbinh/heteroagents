function [w, x] = GaussHermite(order)
    poly = zeros(order + 1, order + 1);
    poly(1, 1) = 1;
    for i = 2:(order + 1)
        tmp_dif = linspace(0, order, order + 1) .* poly(i - 1, :);
        tmp_dif = [tmp_dif(2:end), 0];
        poly(i, :) = 2 * [0, poly(i - 1, 1:order)] - tmp_dif;
    end
    x = roots(fliplr(poly(order + 1, :)));
    x = sort(x); % Put in ascending order
    x_mat = cumprod([ones(order, 1), repmat(x, 1, order)], 2);
    w = 2 ^ (order - 1) * sqrt(pi) * factorial(order) ./ (order ^ 2 .* (x_mat * poly(order, :)') .^ 2);
end
