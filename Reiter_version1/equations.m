% General Equilibrium
% State variables (z, L), Controls (theta, ka, lambda)
% y = g(x): policy, x'=h(x): law of motion of states
n_exo = 1 + crit.m_g(1) * crit.m_g(2); % states
n_shocks = 1; % shocks (usually, the same number as exo states)
n_endo = 2 * crit.n_s * crit.n_k + 1; % controls
n_equ = n_exo + n_endo;

g1 = zeros(n_endo, n_exo); % 1st-order
h1 = zeros(n_exo, n_exo); % 1st-order
g2 = zeros(n_endo, n_exo, n_exo); % 2nd-order
h2 = zeros(n_exo, n_exo, n_exo); % 2nd-order

X = sym('x', [1 n_exo]);
XP = sym('xp', [1 n_exo]);
XSS = [0 Dist'];
Y = sym('y', [1 n_endo]);
YP = sym('yp', [1 n_endo]);
YSS = [reshape(theta.', 1, crit.n_s * crit.n_k) ka_vec.' lambda];
SHOCK = sym('eps', [1 n_shocks]);
SHOCKSS = [0];

% For notation simplicity
THETA = Y(1:crit.n_s * crit.n_k);
KA = Y(crit.n_s * crit.n_k + 1:end);
THETAP = YP(1:crit.n_s * crit.n_k);

% Manually put equations
EQU = sym(zeros(n_equ, 1));

% Euler equations
LHS = vpa(Y(end) * (1 + 2 * para.c1 * (KA.' ./ repmat(para.kgrid, crit.n_s, 1) - (1 - para.delta))));
RHS = sym(zeros(crit.n_s, crit.n_k));

Cheby_s = Chebyshev(para.sgrid, crit.n_s, crit.sbound(1), crit.sbound(2));
for i = 1:crit.n_s
    for j = 1:crit.n_k
        ChebyDiff_k = vpa(ChebyDiff(KA((i - 1) * crit.n_k + j), crit.n_k, crit.kbound(1), crit.kbound(2)));
        for k = 1:crit.n_s
            tmp = Cheby_s(k, :)' * ChebyDiff_k;
            RHS(i, j) = RHS(i, j) + para.Pi_s(i, k) * vpa(THETAP(1:crit.n_s * crit.n_k) * reshape(tmp.', crit.n_s * crit.n_k, 1));
        end
    end
end
RHS = reshape(RHS.', numel(RHS), 1);
EQU(1:crit.n_s * crit.n_k) = vpa(LHS - para.beta * RHS);
%{
% Check SS Euler
EQUSS = double(subs(EQU(1:crit.n_s * crit.n_k), [X, XP, Y, YP, SHOCK], [XSS, XSS, YSS, YSS, SHOCKSS]));
%}
disp('Euler done');

% Generate kn, c0_hat, F_c0_hat and E_c0_hat
[ka, kn, c0_hat, F_c0_hat, E_c0_hat] = otherpolicy_sym(ka_vec, theta, lambda, KA, para.sgrid, para.kgrid, para.Pi_s, Y(end), THETAP, para, crit);
%{
% Check SS policy
[ka_ss, kn_ss, c0_hat_ss, F_c0_hat_ss, E_c0_hat_ss] = otherpolicy(ka_vec, para.sgrid, para.kgrid, para.Pi_s, lambda, theta, para, crit);
tmp = E_c0_hat_ss - double(subs(E_c0_hat_ss, [X, XP, Y, YP, SHOCK], [XSS, XSS, YSS, YSS, SHOCKSS]));
%}
disp('Policy done');

% Value functions
Cheby_s = Chebyshev(para.sgrid, crit.n_s, crit.sbound(1), crit.sbound(2));
Cheby_k = Chebyshev(para.kgrid, crit.n_k, crit.kbound(1), crit.kbound(2));
for i = 1:crit.n_s
    for j = 1:crit.n_k
        tmp_coeff = Cheby_s(i, :)' * Cheby_k(j, :);
        EQU(crit.n_s * crit.n_k + (i - 1) * crit.n_k + j) = vpa(THETA * reshape(tmp_coeff.', crit.n_s * crit.n_k, 1));
    end
end
for i = 1:crit.n_s
    for j = 1:crit.n_k
        idx = crit.n_s * crit.n_k + (i - 1) * crit.n_k + j;
        % disp(idx);
        EQU(idx) = vpa(EQU(idx) - Y(end) * (exp(X(n_exo) + para.sgrid(i)) * para.kgrid(j) ^ para.alpha - para.xi * para.kgrid(j) ^ para.nu + (1 - para.delta) * para.kgrid(j)) ...
        + Y(end) * F_c0_hat(i, j) * (ka(i, j) + para.c1 * para.kgrid(j) * (ka(i, j) / para.kgrid(j) - (1 - para.delta)) ^ 2 + E_c0_hat(i, j) * para.kgrid(j)) ...
        + Y(end) * (1 - F_c0_hat(i, j)) * (kn(i, j) + para.c1 * para.kgrid(j) * (kn(i, j) / para.kgrid(j) - (1 - para.delta)) ^ 2));
        Cheby_ka = Chebyshev_sym(ka(i, j), crit.n_k, crit.kbound(1), crit.kbound(2));
        Cheby_kn = Chebyshev_sym(kn(i, j), crit.n_k, crit.kbound(1), crit.kbound(2));
        for k = 1:crit.n_s
            tmp = Cheby_s(k, :)' * Cheby_ka;
            EQU(idx) = vpa(EQU(idx) - para.beta * F_c0_hat(i, j) * para.Pi_s(i, k) * THETAP * reshape(tmp.', crit.n_s * crit.n_k, 1));
            tmp = Cheby_s(k, :)' * Cheby_kn;
            EQU(idx) = vpa(EQU(idx) - para.beta * (1 - F_c0_hat(i, j)) * para.Pi_s(i, k) * THETAP * reshape(tmp.', crit.n_s * crit.n_k, 1));
        end
    end
end
%{
% Check SS value
EQUSS = double(subs(EQU(crit.n_s * crit.n_k + 1:2 * crit.n_s * crit.n_k), [X, XP, Y, YP, SHOCK], [XSS, XSS, YSS, YSS, SHOCKSS]));
%}
disp('Value done');

% Distribution evolve

% Law of motion of aggregate shock
EQU(2 * crit.n_s * crit.n_k + 1) = XP(n_exo) - para.rho_z * X(n_exo) - para.sigma_z * SHOCK(1);
