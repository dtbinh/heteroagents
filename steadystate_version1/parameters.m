% Import all parameters

% variable para contains all model specific parameters
para.rho_z = 0.95;
para.sigma_z = 0.1;
para.rho_s = 0.9;
para.sigma_s = 0.1;

para.alpha = 0.5;
para.xi = 0.02;
para.nu = 0.7;
para.mu_c = -5.0;
para.sigma_c = 0.01;
para.c1 = 0.1;
para.a = 25.0;

para.beta = 0.95;
para.delta = 0.1;

para.rho_b = 0.9;
para.bbar = 1.0;
para.lambda_b = 1.0;
para.eta = 2.0;

crit.eps = 1e-8;
crit.dampen = 0.1;
crit.kbound = [0.01, 25.0]; % TODO: calculate the upper bound
% crit.sbound = [-1.0, 1.0];
crit.m_g = [10, 20]; %s, k
crit.m_s = 3;
crit.n_s = 5;
crit.n_k = 15;
crit.n_g = 3;

% Compute Gauss-Hermite weights and nodes
[para.tau_s, para.w_s] = GaussHermite(crit.m_s);
% Adjust for Normal
para.w_s = sqrt(2) .* para.w_s;
para.tau_s = para.tau_s / sqrt(pi);

% Compute the upper and lower bound of idio. shock
tmp = para.sigma_s * para.w_s(end) / (1 - para.rho_s);
crit.sbound = [-tmp, tmp];

% Compute grids for Chebyshev collocation
tmp = linspace(crit.n_s, 1, crit.n_s)';
para.sgrid = sec(pi / (2 * crit.n_s)) .* cos(pi / (2 * crit.n_s) * (2 * tmp - 1));
para.sgrid = (para.sgrid + 1) / 2 * (crit.sbound(2) - crit.sbound(1)) + crit.sbound(1);
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

% Compute Gauss-Legendre weights and nodes
[w_s, para.moment_sgrid] = GaussLegendre(crit.m_g(1));
[w_k, para.moment_kgrid] = GaussLegendre(crit.m_g(2));
para.moment_sgrid = (para.moment_sgrid + 1) / 2 * (crit.sbound(2) - crit.sbound(1)) + crit.sbound(1);
para.moment_kgrid = (para.moment_kgrid + 1) / 2 * (crit.kbound(2) - crit.kbound(1)) + crit.kbound(1);
para.tau_g = w_k * w_s'; % Notice that matlab collapses by column
para.tau_g = para.tau_g(:);
crit.m_g(3) = crit.m_g(1) * crit.m_g(2);
para.ggrid = zeros(crit.m_g(3), 2);
tmp = repmat(para.moment_sgrid', crit.m_g(2), 1);
para.ggrid(:, 1) = tmp(:);
para.ggrid(:, 2) = repmat(para.moment_kgrid, crit.m_g(1), 1);
clear tmp w_s w_k;
