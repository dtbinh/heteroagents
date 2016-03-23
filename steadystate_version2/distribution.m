f = @(ka_vec)policy2(ka_vec, para.moment_sgrid, para.moment_kgrid, para.moment_Pi_s, lambda, theta, para, crit);
ka_vec = goldenx(f, repmat(crit.kbound(1), crit.m_g(3), 1), repmat(crit.kbound(2), crit.m_g(3), 1));
[ka, kn, c0_hat, F_c0_hat, E_c0_hat] = otherpolicy(ka_vec, para.moment_sgrid, para.moment_kgrid, para.moment_Pi_s, lambda, theta, para, crit);

kn_vec = reshape(kn', crit.m_g(1) * crit.m_g(2), 1);
F_c0_hat_vec = reshape(F_c0_hat', crit.m_g(1) * crit.m_g(2), 1);
E_c0_hat_vec = reshape(E_c0_hat', crit.m_g(1) * crit.m_g(2), 1);

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

% TODO: I feel like it is not correct... to check
tau_mat = reshape(repmat(para.moment_Pi_s', crit.m_g(2), 1), crit.m_g(1), crit.m_g(1) * crit.m_g(2));
tau_mat = tau_mat';

options = optimoptions('fminunc', 'TolX', crit.eps, 'Display', 'off', 'GradObj', 'on');
moments = zeros(crit.n_g * (crit.n_g + 3) / 2, 1);
last_moments = zeros(crit.n_g * (crit.n_g + 3) / 2, 1);
total_err = 1e5;
iter = 0;
while (total_err > crit.eps && iter < 500)
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
    moments_new(1) = (g0 .* exp(g_value) .* (tau_mat * para.moment_sgrid))' * para.tau_g;
    moments_new(2) = (g0 .* exp(g_value) .* (F_c0_hat_vec .* ka_vec + (1 - F_c0_hat_vec) .* kn_vec))' * para.tau_g;
    for i = 3:crit.n_g * (crit.n_g + 3) / 2
        moments_new(i) = (g0 .* exp(g_value) .* (tau_mat * ((para.moment_sgrid - moments_new(1))) .^ idx(i, 1)) .* (F_c0_hat_vec .* (ka_vec - moments_new(2)) .^ idx(i, 2) + (1 - F_c0_hat_vec) .* (kn_vec - moments_new(2)) .^ idx(i, 2)))' * para.tau_g;
    end

    total_err = sum((moments_new - last_moments) .^ 2);
    disp(total_err);
    moments = moments_new;
    last_moments = moments_new;
    iter = iter + 1;
end

surf(reshape(exp(g_value) .* g0, crit.m_g(2), crit.m_g(1)));
