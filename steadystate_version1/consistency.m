function [g0_inv, jac, moments_new, err, g0, g_value] = consistency(g, moments, idx, para, crit)
    g_value = g(1) * (para.ggrid(:, 1) - moments(1)) + g(2) * (para.ggrid(:, 2) - moments(2));
    for i = 3:crit.n_g * (crit.n_g + 3) / 2
        g_value = g_value + g(i) * ((para.ggrid(:, 1) - moments(1)) .^ idx(i, 1) .* ((para.ggrid(:, 2) - moments(2)) .^ idx(i, 2)) - moments(i));
    end
    g0_inv = exp(g_value)' * para.tau_g;
    g0 = 1 / g0_inv;

    moments_new = zeros(crit.n_g * (crit.n_g + 3) / 2, 1);
    moments_new(1) = (g0 .* exp(g_value) .* para.ggrid(:, 1))' * para.tau_g;
    moments_new(2) = (g0 .* exp(g_value) .* para.ggrid(:, 2))' * para.tau_g;
    for i = 3:crit.n_g * (crit.n_g + 3) / 2
        moments_new(i) = (g0 .* exp(g_value) .* (para.ggrid(:, 1) - moments_new(1)) .^ idx(i, 1) .* ((para.ggrid(:, 2) - moments_new(2)) .^ idx(i, 2)))' * para.tau_g;
    end
    jac = (moments_new - moments) .* g0_inv;
    err = sum((moments_new - moments) .^ 2);
end
