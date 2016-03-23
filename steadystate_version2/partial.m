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
    f = @(ka_vec)policy2(ka_vec, para.sgrid, para.kgrid, para.Pi_s, lambda, theta, para, crit);
    ka_vec = goldenx(f, repmat(crit.kbound(1), crit.n_s * crit.n_k, 1), repmat(crit.kbound(2), crit.n_s * crit.n_k, 1));
%{
    options = optimoptions('fsolve', 'TolX', crit.eps, 'Display', 'off');
    f = @(theta_new)value(theta_new, lambda, theta, ka, kn, c0_hat, para, crit);
    theta_new = fsolve(f, ones(crit.n_s, crit.n_k), options);
%}
    theta_new = value(lambda, theta, ka_vec, para, crit);

    err = sum(sum((theta_new - theta) .^ 2));
    disp(err);
    theta = theta_new; % Could add a dampen coefficient here
    iter = iter + 1;
end
