function [ka, kn, c0_hat, F_c0_hat, E_c0_hat] = otherpolicy(ka_vec, cur_sgrid, cur_kgrid, cur_Pi, lambda, theta, para, crit)
    % Calculate all the other policy variables from unconstrained capital
    k_dim = numel(cur_kgrid);
    s_dim = numel(cur_sgrid);
    % ka and kn
    tmp = reshape(ka_vec, k_dim, s_dim);
    ka = tmp';
    k_mat = repmat(cur_kgrid', s_dim, 1);
    kn = min(ka, (1 - para.delta + para.a) * k_mat);
    kn = max(kn, (1 - para.delta - para.a) * k_mat);
    % c0_hat
    c0_hat = zeros(s_dim, k_dim);
    Cheby_s = Chebyshev(cur_sgrid, crit.n_s, crit.sbound(1), crit.sbound(2));
    for i = 1:s_dim
        for j = 1:k_dim
            c0_hat(i, j) = -lambda * (ka(i, j) - kn(i, j) + para.c1 * cur_kgrid(j) * ((ka(i, j) / cur_kgrid(j) - 1 + para.delta) ^ 2 - (kn(i, j) / cur_kgrid(j) - 1 + para.delta) ^ 2));
            Cheby_k = Chebyshev(ka(i, j), crit.n_k, crit.kbound(1), crit.kbound(2)) - ...
                      Chebyshev(kn(i, j), crit.n_k, crit.kbound(1), crit.kbound(2));
            for k = 1:s_dim
                c0_hat(i, j) = c0_hat(i, j) + para.beta * cur_Pi(i, k) * sum(sum(theta .* (Cheby_s(k, :)' * Cheby_k)));
            end
            c0_hat(i, j) = c0_hat(i, j) / (lambda * cur_kgrid(j));
        end
    end
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
    c0_hat = min(max(c0_hat, 0), para.mu_c);
    E_c0_hat = c0_hat ./ 2;
    F_c0_hat = c0_hat ./ para.mu_c;
end
