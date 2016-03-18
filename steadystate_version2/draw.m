N = 100;

sgrid_dense = para.sgrid;
kgrid_dense = linspace(crit.kbound(1), crit.kbound(2), N);
V = zeros(crit.n_s, N);

for i = 1:crit.n_s
    tmp_cheby_s = Chebyshev(sgrid_dense(i), crit.n_s, crit.sbound(1), crit.sbound(2));
    for j = 1:N
        tmp_cheby_k = Chebyshev(kgrid_dense(j), crit.n_k, crit.kbound(1), crit.kbound(2));
        V(i, j) = sum(sum(tmp_cheby_s' * tmp_cheby_k .* theta));
    end
end

surf(V);
