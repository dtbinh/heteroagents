function [g0_inv, g_value] = consistency2(g, moments, idx, para, crit)
    g_value = g(1) * (para.ggrid(:, 1) - moments(1)) + g(2) * (para.ggrid(:, 2) - moments(2));
    for i = 3:crit.n_g * (crit.n_g + 3) / 2
        g_value = g_value + g(i) * ((para.ggrid(:, 1) - moments(1)) .^ idx(i, 1) .* ((para.ggrid(:, 2) - moments(2)) .^ idx(i, 2)) - moments(i));
    end
    g0_inv = mean(exp(g_value));
end
