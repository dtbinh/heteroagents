% TODO: simplify the redundant codes here
% Calculate policy and probability for simlulation
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
F_c0_hat_vec = reshape(F_c0_hat', crit.m_g(1) * crit.m_g(2), 1);

% Build the transition matrix
Qk = sparse(crit.m_g(1) * crit.m_g(2), crit.m_g(2));
Qka = funbas(fundef({'spli', para.moment_kgrid, 0, 1}), ka_vec);
Qkn = funbas(fundef({'spli', para.moment_kgrid, 0, 1}), kn_vec);
for i = 1:crit.m_g(1) * crit.m_g(2)
    Qk(i, :) = F_c0_hat_vec(i) * Qka(i, :) + (1 - F_c0_hat_vec(i)) * Qkn(i, :);
end
tmp_s = para.rho_s .* repmat(para.moment_sgrid, 1, crit.m_s) + para.sigma_s .* repmat(para.w_s', crit.m_g(1), 1);
Qsp_vec_scat = funbas(fundef({'spli', para.moment_sgrid, 0, 1}), tmp_s(:));
Qsp_vec = sparse(crit.m_g(1), crit.m_g(1));
for i = 1:crit.m_s
    Qsp_vec = Qsp_vec + Qsp_vec_scat((i - 1) * crit.m_g(1) + 1:i * crit.m_g(1), :) .* para.tau_s(i);
end
clear Qka Qkn Qsp_vec_scat tmp_s;

Trans = repmat(Qk, 1, crit.m_g(1)) .* kron(Qsp_vec, ones(crit.m_g(2)));
[tmp, ~] = eigs(Trans.');
Dist = tmp(:, 1) ./ sum(tmp(:, 1));
surf(reshape(Dist, crit.m_g(2), crit.m_g(1)));
