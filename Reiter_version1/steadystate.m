% Main file
clc;
clear;

parameters;
lambda = 0.05;
partial;

ind = 0; % ind = 0 for distribution using bins, otherwise using parametrized
if (ind == 0)
    distribution_bin;
    % Calculate aggregate consumption
    C = exp(para.ggrid(:, 1)) .* (para.ggrid(:, 2) .^ para.alpha) + (1 - para.delta) .* para.ggrid(:, 2);
    C = C - F_c0_hat_vec .* (ka_vec + para.c1 .* para.ggrid(:, 2) .* (ka_vec ./ para.ggrid(:, 2) - 1 + para.delta) .^ 2 + E_c0_hat_vec .* para.ggrid(:, 2));
    C = C - (1 - F_c0_hat_vec) .* (kn_vec + para.c1 .* para.ggrid(:, 2) .* (kn_vec ./ para.ggrid(:, 2) - 1 + para.delta) .^ 2);
    C = C' * Dist;
    
    lambda = C ^ -para.eta;
    ka_vec_dist = ka_vec;
    partial; % Generate steady state coefficients again
    save('SteadyStateResults.mat', 'lambda', 'theta', 'ka_vec', 'ka_vec_dist', 'Trans', 'Dist');
else
    distribution;
    % Calculate aggregate consumption
    C = exp(para.ggrid(:, 1)) .* (para.ggrid(:, 2) .^ para.alpha) + (1 - para.delta) .* para.ggrid(:, 2);
    C = C - F_c0_hat_vec .* (ka_vec + para.c1 .* para.ggrid(:, 2) .* (ka_vec ./ para.ggrid(:, 2) - 1 + para.delta) .^ 2 + E_c0_hat_vec .* para.ggrid(:, 2));
    C = C - (1 - F_c0_hat_vec) .* (kn_vec + para.c1 .* para.ggrid(:, 2) .* (kn_vec ./ para.ggrid(:, 2) - 1 + para.delta) .^ 2);
    C = (g0 .* exp(g_value) .* C)' * para.tau_g;
end
