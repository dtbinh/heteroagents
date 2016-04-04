function res = value(lambda, theta, ka_vec, para, crit)
    % Solve out new value function
    % Construct intermediate variables
    [ka, kn, ~, F_c0_hat, E_c0_hat] = otherpolicy(ka_vec, para.sgrid, para.kgrid, para.Pi_s, lambda, theta, para, crit);

    coeff = zeros(crit.n_s * crit.n_k, crit.n_s * crit.n_k);
    const = zeros(crit.n_s * crit.n_k, 1);
    Cheby_s = Chebyshev(para.sgrid, crit.n_s, crit.sbound(1), crit.sbound(2));
    Cheby_k = Chebyshev(para.kgrid, crit.n_k, crit.kbound(1), crit.kbound(2));
    for i = 1:crit.n_s
        for j = 1:crit.n_k
            tmp_coeff = Cheby_s(i, :)' * Cheby_k(j, :);
            coeff((i - 1) * crit.n_k + j, :) = reshape(tmp_coeff', 1, crit.n_s * crit.n_k);
        end
    end
    for i = 1:crit.n_s
        for j = 1:crit.n_k
            const((i - 1) * crit.n_k + j) = lambda * (exp(para.sgrid(i)) * para.kgrid(j) ^ para.alpha - para.xi * para.kgrid(j) ^ para.nu + (1 - para.delta) * para.kgrid(j)) ...
            - lambda * F_c0_hat(i, j) * (ka(i, j) + para.c1 * para.kgrid(j) * (ka(i, j) / para.kgrid(j) - (1 - para.delta)) ^ 2 + E_c0_hat(i, j) * para.kgrid(j)) ...
            - lambda * (1 - F_c0_hat(i, j)) * (kn(i, j) + para.c1 * para.kgrid(j) * (kn(i, j) / para.kgrid(j) - (1 - para.delta)) ^ 2);
            Cheby_ka = Chebyshev(ka(i, j), crit.n_k, crit.kbound(1), crit.kbound(2));
            Cheby_kn = Chebyshev(kn(i, j), crit.n_k, crit.kbound(1), crit.kbound(2));

            % This one treats right hand side as unknown theta
            for k = 1:crit.n_s
                tmp_coeff = Cheby_s(k, :)' * Cheby_ka;
                coeff((i - 1) * crit.n_k + j, :) = coeff((i - 1) * crit.n_k + j, :) - para.beta * F_c0_hat(i, j) * para.Pi_s(i, k) * reshape(tmp_coeff', 1, crit.n_s * crit.n_k);
                tmp_coeff = Cheby_s(k, :)' * Cheby_kn;
                coeff((i - 1) * crit.n_k + j, :) = coeff((i - 1) * crit.n_k + j, :) - para.beta * (1 - F_c0_hat(i, j)) * para.Pi_s(i, k) * reshape(tmp_coeff', 1, crit.n_s * crit.n_k);
            end
%{
            % This one treats right hand side as given last time
            for k = 1:crit.n_s
                const((i - 1) * crit.n_k + j) = const((i - 1) * crit.n_k + j) + para.beta * F_c0_hat(i, j) * para.Pi_s(i, k) * sum(sum(Cheby_s(k, :)' * Cheby_ka .* theta));
                const((i - 1) * crit.n_k + j) = const((i - 1) * crit.n_k + j) + para.beta * (1 - F_c0_hat(i, j)) * para.Pi_s(i, k) * sum(sum(Cheby_s(k, :)' * Cheby_kn .* theta));
            end
%}
        end
    end
    % coeff * theta = const
    res = reshape(coeff \ const, crit.n_k, crit.n_s);
    res = res';

end
