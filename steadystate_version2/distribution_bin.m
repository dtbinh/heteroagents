% TODO: simplify the redundant codes here
% Calculate policy and probability for simlulation
f = @(ka_vec)policy2(ka_vec, lambda, theta, para, crit, 'Simul');
ka_vec = goldenx(f, repmat(crit.kbound(1), crit.m_g(3), 1), repmat(crit.kbound(2), crit.m_g(3), 1));
tmp = reshape(ka_vec, crit.m_g(2), crit.m_g(1));
ka = tmp';

k_mat = repmat(para.moment_kgrid', crit.m_g(1), 1);
kn = min(ka, (1 - para.delta + para.a) * k_mat);
kn = max(kn, (1 - para.delta - para.a) * k_mat);
kn_vec = reshape(kn', crit.m_g(1) * crit.m_g(2), 1);

c0_hat = zeros(crit.m_g(1), crit.m_g(2));
Cheby_s = Chebyshev(para.moment_sgrid, crit.n_s, crit.sbound(1), crit.sbound(2));
for i = 1:crit.m_g(1)
    for j = 1:crit.m_g(2)
        c0_hat(i, j) = -lambda * (ka(i, j) - kn(i, j) + para.c1 * para.moment_kgrid(j) * ((ka(i, j) / para.moment_kgrid(j) - 1 + para.delta) ^ 2 - (kn(i, j) / para.moment_kgrid(j) - 1 + para.delta) ^ 2));
        Cheby_k = Chebyshev(ka(i, j), crit.n_k, crit.kbound(1), crit.kbound(2)) - ...
                  Chebyshev(kn(i, j), crit.n_k, crit.kbound(1), crit.kbound(2));
        for k = 1:crit.m_g(1)
            c0_hat(i, j) = c0_hat(i, j) + para.beta * para.moment_Pi_s(i, k) * sum(sum(theta .* (Cheby_s(k, :)' * Cheby_k)));
        end
        c0_hat(i, j) = c0_hat(i, j) / (lambda * para.moment_kgrid(j));
    end
end
c0_hat = max(c0_hat, crit.eps);
F_c0_hat = logncdf(c0_hat, para.mu_c, para.sigma_c);
F_c0_hat_vec = reshape(F_c0_hat', crit.m_g(1) * crit.m_g(2), 1);
f = @(x)lognpdf(x, para.mu_c, para.sigma_c) .* x;
E_c0_hat = zeros(crit.m_g(1), crit.m_g(2));
for i = 1:crit.m_g(1)
    for j = 1:crit.m_g(2)
        E_c0_hat(i, j) = integral(f, 0, c0_hat(i, j));
    end
end
E_c0_hat_vec = reshape(F_c0_hat', crit.m_g(1) * crit.m_g(2), 1);

% Build the transition matrix
Qk = sparse(crit.m_g(1) * crit.m_g(2), crit.m_g(2));
Qka = funbas(fundef({'spli', para.moment_kgrid, 0, 1}), ka_vec);
Qkn = funbas(fundef({'spli', para.moment_kgrid, 0, 1}), kn_vec);
for i = 1:crit.m_g(1) * crit.m_g(2)
    Qk(i, :) = F_c0_hat_vec(i) * Qka(i, :) + (1 - F_c0_hat_vec(i)) * Qkn(i, :);
end
Qs = kron(para.moment_Pi_s, ones(crit.m_g(2), 1));
Trans = dprod(Qs, Qk);

[tmp, ~] = eigs(Trans.');
Dist = tmp(:, 1) ./ sum(tmp(:, 1));
figure(1);
surf(reshape(Dist, crit.m_g(2), crit.m_g(1))');

% find a parametrized distribution
idx = zeros(crit.n_g * (crit.n_g + 3) / 2, 2);
tmp = 0;
for i = 1:crit.n_g
    for j = 0:i
        tmp = tmp + 1;
        idx(tmp, 1) = i - j;
        idx(tmp, 2) = j;
    end
end
moments = zeros(crit.n_g * (crit.n_g + 3) / 2, 1);
moments(1) = para.ggrid(:, 1)' * Dist;
moments(2) = para.ggrid(:, 2)' * Dist;
for i = 3:crit.n_g * (crit.n_g + 3) / 2
    moments(i) = ((para.ggrid(:, 1) - moments(1)) .^ idx(i, 1) .* ((para.ggrid(:, 2) - moments(2)) .^ idx(i, 2)))' * Dist;
end

options = optimoptions('fminunc', 'TolX', crit.eps, 'Display', 'off', 'MaxFunEvals', 10000);
f = @(g) consistency2(g, moments, idx, para, crit);
g = fminunc(f, zeros(crit.n_g * (crit.n_g + 3) / 2, 1), options);
[g0_inv, g_value] = consistency2(g, moments, idx, para, crit);
g_value = exp(g_value) ./ g0_inv ./ numel(g_value);
figure(2);
surf(reshape(g_value, crit.m_g(2), crit.m_g(1))');

% Compare simulation moments and parametrized moments
moments_para = zeros(crit.n_g * (crit.n_g + 3) / 2, 1);
moments_para(1) = para.ggrid(:, 1)' * g_value;
moments_para(2) = para.ggrid(:, 2)' * g_value;
for i = 3:crit.n_g * (crit.n_g + 3) / 2
    moments_para(i) = ((para.ggrid(:, 1) - moments_para(1)) .^ idx(i, 1) .* ((para.ggrid(:, 2) - moments_para(2)) .^ idx(i, 2)))' * g_value;
end
