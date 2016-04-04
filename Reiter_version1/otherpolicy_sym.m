function [ka, kn, c0_hat, F_c0_hat, E_c0_hat] = otherpolicy_sym(ka_vec_ss, theta_ss, lambda_ss, ka_vec, cur_sgrid, cur_kgrid, cur_Pi, lambda, theta, para, crit)
    % Calculate all the other policy variables from unconstrained capital
    [~, kn_ss, c0_hat_ss, ~, ~] = otherpolicy(ka_vec_ss, cur_sgrid, cur_kgrid, cur_Pi, lambda_ss, theta_ss, para, crit);
    k_dim = numel(cur_kgrid);
    s_dim = numel(cur_sgrid);
    % ka and kn
    tmp = reshape(ka_vec, k_dim, s_dim);
    ka = tmp.';
    k_mat = repmat(cur_kgrid', s_dim, 1);
    % Calculate SS choices
    ind_l = (kn_ss == (1 - para.delta - para.a) * k_mat);
    ind_r = (kn_ss == (1 - para.delta + para.a) * k_mat);
    kn = ind_l .* (1 - para.delta - para.a) .* k_mat + ind_r .* (1 - para.delta + para.a) .* k_mat ...
         + (1 - ind_l) .* (1 - ind_r) .* ka;
    % c0_hat
    c0_hat = sym(zeros(s_dim, k_dim));
    Cheby_s = Chebyshev(cur_sgrid, crit.n_s, crit.sbound(1), crit.sbound(2));
    for i = 1:s_dim
        for j = 1:k_dim
            % disp((i - 1) * k_dim + j);
            c0_hat(i, j) = -lambda * (ka(i, j) - kn(i, j) + para.c1 * cur_kgrid(j) * ((ka(i, j) / cur_kgrid(j) - 1 + para.delta) ^ 2 - (kn(i, j) / cur_kgrid(j) - 1 + para.delta) ^ 2));
            Cheby_k = Chebyshev_sym(ka(i, j), crit.n_k, crit.kbound(1), crit.kbound(2)) - ...
                      Chebyshev_sym(kn(i, j), crit.n_k, crit.kbound(1), crit.kbound(2));
            for k = 1:s_dim
                tmp = Cheby_s(k, :).' * Cheby_k;
                c0_hat(i, j) = c0_hat(i, j) + para.beta * cur_Pi(i, k) * vpa(theta * reshape(tmp.', s_dim * k_dim, 1));
            end
            c0_hat(i, j) = c0_hat(i, j) / (lambda * cur_kgrid(j));
        end
    end

    disp('Simplify');
    c0_hat = simplify(c0_hat, 'IgnoreAnalyticConstraints', true);
    disp('Simplification done!');

%{
% This part is for log-normal
    c0_hat = max(c0_hat, crit.eps);
    F_c0_hat = logncdf(c0_hat, para.mu_c, para.sigma_c);
    f = @(x)lognpdf(x, para.mu_c, para.sigma_c) .* x;
    E_c0_hat = zeros(s_dim, k_dim);
    for i = 1:s_dim
        for j = 1:k_dim
            E_c0_hat(i, j) = integral(f, 0, c0_hat(i, j));
        end
    end
%}
% Try a uniform one
    ind_l = (c0_hat_ss == 0);
    ind_r = (c0_hat_ss == para.mu_c);
    c0_hat = ind_r .* para.mu_c + (1 - ind_l) .* (1 - ind_r) .* c0_hat;
    F_c0_hat = vpa(c0_hat ./ para.mu_c);
    E_c0_hat = vpa(c0_hat ./ 2.0);
end
