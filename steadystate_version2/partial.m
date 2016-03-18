lambda = 0.05;
theta = zeros(crit.n_s, crit.n_k);

err = 1e6;
iter = 0;
while (err > crit.eps)
%{
    ka = zeros(crit.n_s, crit.n_k);
    ka_vec = reshape(ka', numel(ka), 1);
    options = optimoptions('fsolve', 'TolX', crit.eps, 'Display', 'off');
    f = @(ka_vec)policy(ka_vec, lambda, theta, para, crit);
    ka_vec = fsolve(f, ka_vec, options);

    ka_vec = max(min(ka_vec, crit.kbound(2)), crit.kbound(1)); % Restrict the choice set
%}
    f = @(ka_vec)policy2(ka_vec, lambda, theta, para, crit, 'Calc');
    ka_vec = goldenx(f, repmat(crit.kbound(1), crit.n_s * crit.n_k, 1), repmat(crit.kbound(2), crit.n_s * crit.n_k, 1));

    tmp = reshape(ka_vec, crit.n_k, crit.n_s);
    ka = tmp';
    % disp(ka);

    k_mat = repmat(para.kgrid', crit.n_s, 1);
    kn = min(ka, (1 - para.delta + para.a) * k_mat);
    kn = max(kn, (1 - para.delta - para.a) * k_mat);

    c0_hat = zeros(crit.n_s, crit.n_k);
    Cheby_s = Chebyshev(para.sgrid, crit.n_s, crit.sbound(1), crit.sbound(2));
    for i = 1:crit.n_s
        for j = 1:crit.n_k
            c0_hat(i, j) = -lambda * (ka(i, j) - kn(i, j) + para.c1 * para.kgrid(j) * ((ka(i, j) / para.kgrid(j) - 1 + para.delta) ^ 2 - (kn(i, j) / para.kgrid(j) - 1 + para.delta) ^ 2));
            Cheby_k = Chebyshev(ka(i, j), crit.n_k, crit.kbound(1), crit.kbound(2)) - ...
                      Chebyshev(kn(i, j), crit.n_k, crit.kbound(1), crit.kbound(2));
            for k = 1:crit.n_s
                c0_hat(i, j) = c0_hat(i, j) + para.beta * para.Pi_s(i, k) * sum(sum(theta .* (Cheby_s(k, :)' * Cheby_k)));
            end
            c0_hat(i, j) = c0_hat(i, j) / (lambda * para.kgrid(j));
        end
    end
    c0_hat = max(c0_hat, crit.eps);

%{
    options = optimoptions('fsolve', 'TolX', crit.eps, 'Display', 'off');
    f = @(theta_new)value(theta_new, lambda, theta, ka, kn, c0_hat, para, crit);
    theta_new = fsolve(f, ones(crit.n_s, crit.n_k), options);
%}
    theta_new = value(lambda, theta, ka, kn, c0_hat, para, crit);

    err = sum(sum((theta_new - theta) .^ 2));
    disp(err);
    theta = theta_new; % Could add a dampen coefficient here
    iter = iter + 1;
end
