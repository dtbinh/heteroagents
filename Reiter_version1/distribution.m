f = @(ka_vec)policy2(ka_vec, para.moment_sgrid, para.moment_kgrid, para.moment_Pi_s, lambda, theta, para, crit);
ka_vec = goldenx(f, repmat(crit.kbound(1), crit.m_g(3), 1), repmat(crit.kbound(2), crit.m_g(3), 1));
[ka, kn, c0_hat, F_c0_hat, E_c0_hat] = otherpolicy(ka_vec, para.moment_sgrid, para.moment_kgrid, para.moment_Pi_s, lambda, theta, para, crit);

kn_vec = reshape(kn', crit.m_g(1) * crit.m_g(2), 1);
F_c0_hat_vec = reshape(F_c0_hat', crit.m_g(1) * crit.m_g(2), 1);
E_c0_hat_vec = reshape(E_c0_hat', crit.m_g(1) * crit.m_g(2), 1);

% Moments convergence
tau_mat = reshape(repmat(para.moment_Pi_s', crit.m_g(2), 1), crit.m_g(1), crit.m_g(1) * crit.m_g(2));
tau_mat = tau_mat';
w_s_vec = reshape(repmat(para.w_s', crit.m_g(2), 1), crit.m_g(3), 1);

options = optimoptions('fminunc', 'TolX', crit.eps, 'Display', 'off', 'MaxFunEvals', 50000, 'GradObj', 'on');
% Initialize moments with a flat distribution
moments = zeros(crit.m_g(1), crit.n_g);
for i = 1:crit.m_g(1)
    moments(i, 1) = mean(para.moment_kgrid);
    for j = 2:crit.n_g
        moments(i, j) = mean((para.moment_kgrid - moments(i, 1)) .^ j);
    end
end
g = zeros(crit.m_g(1), crit.n_g);
g0 = zeros(crit.m_g(1), 1);
g_value = zeros(crit.m_g(3), 1);
total_err = 1e5;
iter = 0;
% TODO: I don't like my code. Ugly.
while (total_err > 2e-2)
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

    % Generate next period moments
    % TODO: ugly code
    moments_new = zeros(crit.m_g(1), crit.n_g);
    for i = 1:crit.m_g(1)
        moments_new(i, 1) = crit.m_g(1) * mean(tau_mat(:, i) .* (F_c0_hat_vec .* ka_vec + (1 - F_c0_hat_vec) .* kn_vec) .* w_s_vec .* g_value) / para.w_s(i);
        for j = 2:crit.n_g
            moments_new(i, j) = crit.m_g(1) * mean(tau_mat(:, i) .* (F_c0_hat_vec .* (ka_vec - moments_new(i, 1)) .^ j + (1 - F_c0_hat_vec) .* (kn_vec - moments_new(i, 1)) .^ j) .* w_s_vec .* g_value) / para.w_s(i);
        end
    end

    total_err = max(max(abs(moments_new - moments)));
    disp(total_err);
    moments = moments_new;
    iter = iter + 1;
end
figure(1);
surf(reshape(g_value .* w_s_vec, crit.m_g(2), crit.m_g(1))');
figure(2);
plot(para.moment_kgrid, reshape(g_value .* w_s_vec, crit.m_g(2), crit.m_g(1)));
legend('show');