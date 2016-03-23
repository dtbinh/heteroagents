% Calculate policy and probability for simlulation
f = @(ka_vec)policy2(ka_vec, para.moment_sgrid, para.moment_kgrid, para.moment_Pi_s, lambda, theta, para, crit);
ka_vec = goldenx(f, repmat(crit.kbound(1), crit.m_g(3), 1), repmat(crit.kbound(2), crit.m_g(3), 1));
[ka, kn, c0_hat, F_c0_hat, E_c0_hat] = otherpolicy(ka_vec, para.moment_sgrid, para.moment_kgrid, para.moment_Pi_s, lambda, theta, para, crit);

kn_vec = reshape(kn', crit.m_g(1) * crit.m_g(2), 1);
F_c0_hat_vec = reshape(F_c0_hat', crit.m_g(1) * crit.m_g(2), 1);
E_c0_hat_vec = reshape(E_c0_hat', crit.m_g(1) * crit.m_g(2), 1);

% Build the transition matrix
Qk = sparse(crit.m_g(3), crit.m_g(2));
Qka = funbas(fundef({'spli', para.moment_kgrid, 0, 1}), ka_vec);
Qkn = funbas(fundef({'spli', para.moment_kgrid, 0, 1}), kn_vec);
for i = 1:crit.m_g(3)
    Qk(i, :) = F_c0_hat_vec(i) * Qka(i, :) + (1 - F_c0_hat_vec(i)) * Qkn(i, :);
end
Qs = kron(para.moment_Pi_s, ones(crit.m_g(2), 1));
Trans = dprod(Qs, Qk);

[tmp, ~] = eigs(Trans.', 1);
tmp(abs(tmp) < 1e-10) = 0;
Dist = tmp ./ sum(tmp);
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

% Compare parametrized moments and moments after updating
tau_mat = reshape(repmat(para.moment_Pi_s', crit.m_g(2), 1), crit.m_g(1), crit.m_g(1) * crit.m_g(2));
tau_mat = tau_mat';

moments_new = zeros(crit.n_g * (crit.n_g + 3) / 2, 1);
moments_new(1) = (tau_mat * para.moment_sgrid)' * g_value;
moments_new(2) = (F_c0_hat_vec .* ka_vec + (1 - F_c0_hat_vec) .* kn_vec)' * g_value;
for i = 3:crit.n_g * (crit.n_g + 3) / 2
    moments_new(i) = ((tau_mat * ((para.moment_sgrid - moments_new(1))) .^ idx(i, 1)) .* (F_c0_hat_vec .* (ka_vec - moments_new(2)) .^ idx(i, 2) + (1 - F_c0_hat_vec) .* (kn_vec - moments_new(2)) .^ idx(i, 2)))' * g_value;
end
