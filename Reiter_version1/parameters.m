% Import all parameters

% variable para contains all model specific parameters
para.rho_z = 0.95;
para.sigma_z = 0.1;
para.rho_s = 0.9;
para.sigma_s = 0.1;

para.alpha = 0.5;
para.xi = 0.02;
para.nu = 0.7;
para.mu_c = 1.0;
para.sigma_c = 2.0;
para.c1 = 1.5;
para.a = 0.2;

para.beta = 0.95;
para.delta = 0.1;

para.rho_b = 0.9;
para.bbar = 1.0;
para.lambda_b = 1.0;
para.eta = 2.0;

crit.eps = 1e-8;
crit.dampen = 0.1;
crit.kbound = [1, 20.0]; % TODO: calculate the upper bound
crit.m_g = [7, 20]; % s, k
crit.n_s = 5;
crit.n_k = 5;
crit.n_g = 4;

% Compute Rouwenhorst approximations
[para.sgrid, para.Pi_s] = Rouwenhorst(para.rho_s, para.sigma_s, crit.n_s);
crit.sbound = [para.sgrid(1), para.sgrid(end)];

% Compute grids for Chebyshev collocation
tmp = linspace(crit.n_k, 1, crit.n_k)';
para.kgrid = sec(pi / (2 * crit.n_k)) .* cos(pi / (2 * crit.n_k) * (2 * tmp - 1));
para.kgrid = (para.kgrid + 1) / 2 * (crit.kbound(2) - crit.kbound(1)) + crit.kbound(1);
%{
% Alternative choice of capital grid
curv = 0.4;
para.kgrid = linspace(0, (crit.kbound(2) - crit.kbound(1)) ^ curv, crit.n_k) .^ (1 / curv) + crit.kbound(1);
para.kgrid = para.kgrid';
clear curv;
%}

[para.moment_sgrid, para.moment_Pi_s] = Rouwenhorst(para.rho_s, para.sigma_s, crit.m_g(1));
[para.w_s, ~] = eigs(para.moment_Pi_s.', 1);
para.w_s = para.w_s ./ sum(para.w_s);
%{
% Compute Gauss-Legendre weights and nodes
[w_k, para.moment_kgrid] = GaussLegendre(crit.m_g(2));
para.moment_kgrid = (para.moment_kgrid + 1) / 2 * (crit.kbound(2) - crit.kbound(1)) + crit.kbound(1);
para.tau_g = w_k * w_s'; % Notice that matlab collapses by column
para.tau_g = para.tau_g(:);
%}
para.moment_kgrid = linspace(crit.kbound(1), crit.kbound(2), crit.m_g(2))';
para.tau_g = ones(crit.m_g(1) * crit.m_g(2), 1) ./ crit.m_g(2);
crit.m_g(3) = crit.m_g(1) * crit.m_g(2);
para.ggrid = zeros(crit.m_g(3), 2);
tmp = repmat(para.moment_sgrid', crit.m_g(2), 1);
para.ggrid(:, 1) = tmp(:);
para.ggrid(:, 2) = repmat(para.moment_kgrid, crit.m_g(1), 1);
clear tmp w_s w_k;

save Parameters.mat;
