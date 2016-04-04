function res = policy(ka_vec, lambda, theta, para, crit)
    % Calculate ka
    LHS = lambda * (1 + 2 * para.c1 * (ka_vec ./ repmat(para.kgrid, crit.n_s, 1) - (1 - para.delta)));
    RHS = zeros(crit.n_s, crit.n_k);

    Cheby_s = Chebyshev(para.sgrid, crit.n_s, crit.sbound(1), crit.sbound(2));
    for i = 1:crit.n_s
        for j = 1:crit.n_k
            ChebyDiff_k = ChebyDiff(ka_vec((i - 1) * crit.n_k + j), crit.n_k, crit.kbound(1), crit.kbound(2));
            for k = 1:crit.n_s
                RHS(i, j) = RHS(i, j) + para.Pi_s(i, k) * sum(sum(theta .* (Cheby_s(k, :)' * ChebyDiff_k)));
            end
        end
    end
    RHS = reshape(RHS', numel(RHS), 1);
    res = LHS - para.beta * RHS;
end
