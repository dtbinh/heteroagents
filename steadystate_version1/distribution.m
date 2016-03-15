% TODO: Simplify the redundant codes here
f = @(ka_vec)policy2(ka_vec, lambda, theta, para, crit, 'Simul');
ka_vec = goldenx(f, repmat(crit.kbound(1), crit.m_g(3), 1), repmat(crit.kbound(2), crit.m_g(3), 1));
tmp = reshape(ka_vec, crit.m_g(2), crit.m_g(1));
ka = tmp';

k_mat = repmat(para.moment_kgrid', crit.m_g(1), 1);
kn = min(ka, (1 - para.delta) * k_mat + para.a);
kn = max(kn, max(crit.kbound(1), (1 - para.delta) * k_mat - para.a));
kn_vec = reshape(kn', crit.m_g(1) * crit.m_g(2), 1);

c0_hat = zeros(crit.m_g(1), crit.m_g(2));
for i = 1:crit.m_g(1)
    tmp_s = para.rho_s .* para.moment_sgrid(i) + para.sigma_s .* para.w_s;
    Cheby_s = Chebyshev(tmp_s, crit.n_s, crit.sbound(1), crit.sbound(2));
    for j = 1:crit.m_g(2)
        c0_hat(i, j) = -lambda * (ka(i, j) - kn(i, j) + para.c1 * para.moment_kgrid(j) * ((ka(i, j) / para.moment_kgrid(j) - 1 + para.delta) ^ 2 - (kn(i, j) / para.moment_kgrid(j) - 1 + para.delta) ^ 2));
        Cheby_k = Chebyshev(ka(i, j), crit.n_k, crit.kbound(1), crit.kbound(2)) - ...
                  Chebyshev(kn(i, j), crit.n_k, crit.kbound(1), crit.kbound(2));
        for k = 1:crit.m_s
            c0_hat(i, j) = c0_hat(i, j) + para.beta * para.tau_s(k) * sum(sum(theta .* (Cheby_s(k, :)' * Cheby_k)));
        end
        c0_hat(i, j) = c0_hat(i, j) / (lambda * para.moment_kgrid(j));
    end
end
c0_hat = max(c0_hat, crit.eps);
F_c0_hat = logncdf(c0_hat, para.mu_c, para.sigma_c);
f = @(x)lognpdf(x, para.mu_c, para.sigma_c) .* x;
E_c0_hat = zeros(crit.m_g(1), crit.m_g(2));
for i = 1:crit.m_g(1)
    for j = 1:crit.m_g(2)
        E_c0_hat(i, j) = integral(f, 0, c0_hat(i, j));
    end
end
F_c0_hat_vec = reshape(F_c0_hat', crit.m_g(1) * crit.m_g(2), 1);
E_c0_hat_vec = reshape(E_c0_hat', crit.m_g(1) * crit.m_g(2), 1);

sp_mat = para.rho_s * repmat(para.ggrid(:, 1), 1, crit.m_s) + para.sigma_s * repmat(para.w_s', crit.m_g(3), 1);

% Generate index
idx = zeros(crit.n_g * (crit.n_g + 3) / 2, 2);
tmp = 0;
for i = 1:crit.n_g
    for j = 0:i
        tmp = tmp + 1;
        idx(tmp, 1) = i - j;
        idx(tmp, 2) = j;
    end
end

options = optimoptions('fminunc', 'TolX', crit.eps, 'Display', 'off', 'GradObj', 'on');
moments = zeros(crit.n_g * (crit.n_g + 3) / 2, 1);
last_moments = zeros(crit.n_g * (crit.n_g + 3) / 2, 1);
total_err = 1e5;
iter = 0;
while (total_err > crit.eps && iter < 100)
    f = @(g) consistency(g, moments, idx, para, crit);
    g = fminunc(f, zeros(crit.n_g * (crit.n_g + 3) / 2, 1), options);
    err = 1e5;
    sub_iter = 0;
    while (err > crit.eps && sub_iter < 200)
        [g0_inv, ~, moments_new, err, g0, g_value] = consistency(g, moments, idx, para, crit);
        moments = moments_new * (1 - crit.dampen) + moments * crit.dampen;
        sub_iter = sub_iter + 1;
    end
    % Generate next period moments
    moments_new = zeros(crit.n_g * (crit.n_g + 3) / 2, 1);
    moments_new(1) = (g0 .* exp(g_value) .* (sp_mat * para.tau_s))' * para.tau_g;
    moments_new(2) = (g0 .* exp(g_value) .* (F_c0_hat_vec .* ka_vec + (1 - F_c0_hat_vec) .* kn_vec))' * para.tau_g;
    for i = 3:crit.n_g * (crit.n_g + 3) / 2
        moments_new(i) = (g0 .* exp(g_value) .* ((sp_mat - moments_new(1)) .^ idx(i, 1) * para.tau_s) .* (F_c0_hat_vec .* (ka_vec - moments_new(2)) .^ idx(i, 2) + (1 - F_c0_hat_vec) .* (kn_vec - moments_new(2)) .^ idx(i, 2)))' * para.tau_g;
    end

    total_err = sum((moments_new - last_moments) .^ 2);
    disp(total_err);
    moments = moments_new;
    last_moments = moments_new;
    iter = iter + 1;
end
