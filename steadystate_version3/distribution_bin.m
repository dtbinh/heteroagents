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
w_s_vec = reshape(repmat(para.w_s', crit.m_g(2), 1), crit.m_g(3), 1);

moments = zeros(crit.m_g(1), crit.n_g);
for i = 1:crit.m_g(1)
    moments(i, 1) = para.moment_kgrid' * Dist((i - 1) * crit.m_g(2)+1:i*crit.m_g(2)) ./ para.w_s(i);
    for j = 2:crit.n_g
        moments(i, j) = ((para.moment_kgrid' - moments(i, 1)) .^ j) * Dist((i - 1) * crit.m_g(2)+1:i*crit.m_g(2)) ./ para.w_s(i);
    end
end

options = optimoptions('fminunc', 'TolX', crit.eps, 'Display', 'off', 'MaxFunEvals', 50000, 'GradObj', 'on');
for i = 1:crit.m_g(1)
    tmp_para = zeros(crit.m_g(2), crit.n_g);
    tmp_para(:, 1) = para.moment_kgrid - moments(i, 1);
    for j = 2:crit.n_g
        tmp_para(:, j) = (para.moment_kgrid - moments(i, 1)) .^ j - moments(i, j);
    end
    f = @(x) consistency(tmp_para, x, crit.n_g);
    [g(i, :), g0(i)] = fminunc(f, zeros(1, crit.n_g), options);
    g0(i) = 1 ./ g0(i);
    g_value((i - 1) * crit.m_g(2) + 1:i * crit.m_g(2)) = exp(tmp_para * g(i, :)') .* g0(i);
end
figure(2);
surf(reshape(g_value .* w_s_vec, crit.m_g(2), crit.m_g(1))');
