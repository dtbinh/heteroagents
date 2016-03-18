function res = policy2(ka_vec, lambda, theta, para, crit, req)
    % Return the maximization function
    res = zeros(numel(ka_vec), 1);
    switch req
        case 'Calc'
            Cheby_ka = Chebyshev(ka_vec, crit.n_k, crit.kbound(1), crit.kbound(2));
            Cheby_s = Chebyshev(para.sgrid, crit.n_s, crit.sbound(1), crit.sbound(2));
            for i = 1:crit.n_s
                for j = 1:crit.n_k
                    idx = (i - 1) * crit.n_k + j;
                    res(idx) = -lambda * (ka_vec(idx) + para.c1 * para.kgrid(j) * (ka_vec(idx) / para.kgrid(j) - (1 - para.delta)) ^ 2);
                    for k = 1:crit.n_s
                        res(idx) = res(idx) + para.beta * para.Pi_s(i, k) * sum(sum(theta .* (Cheby_s(k, :)' * Cheby_ka(idx, :))));
                    end
                end
            end
        case 'Simul'
            Cheby_ka = Chebyshev(ka_vec, crit.n_k, crit.kbound(1), crit.kbound(2));
            Cheby_s = Chebyshev(para.moment_sgrid, crit.n_s, crit.sbound(1), crit.sbound(2));
            for i = 1:crit.m_g(1)
                for j = 1:crit.m_g(2)
                    idx = (i - 1) * crit.m_g(2) + j;
                    res(idx) = -lambda * (ka_vec(idx) + para.c1 * para.moment_kgrid(j) * (ka_vec(idx) / para.moment_kgrid(j) - (1 - para.delta)) ^ 2);
                    for k = 1:crit.m_g(1)
                        res(idx) = res(idx) + para.beta * para.moment_Pi_s(i, k) * sum(sum(theta .* (Cheby_s(k, :)' * Cheby_ka(idx, :))));
                    end
                end
            end
    end
end
