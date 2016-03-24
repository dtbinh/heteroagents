function res = policy2(ka_vec, cur_sgrid, cur_kgrid, cur_Pi, lambda, theta, para, crit)
    % Return the maximization function
    res = zeros(numel(ka_vec), 1);
    k_dim = numel(cur_kgrid);
    s_dim = numel(cur_sgrid);
    Cheby_ka = Chebyshev(ka_vec, crit.n_k, crit.kbound(1), crit.kbound(2));
    Cheby_s = Chebyshev(cur_sgrid, crit.n_s, crit.sbound(1), crit.sbound(2));
    for i = 1:s_dim
        for j = 1:k_dim
            idx = (i - 1) * k_dim + j;
            res(idx) = -lambda * (ka_vec(idx) + para.c1 * cur_kgrid(j) * (ka_vec(idx) / cur_kgrid(j) - (1 - para.delta)) ^ 2);
            for k = 1:s_dim
                res(idx) = res(idx) + para.beta * cur_Pi(i, k) * sum(sum(theta .* (Cheby_s(k, :)' * Cheby_ka(idx, :))));
            end
        end
    end
end
