% Solve the dynamic model
% Distribution using bins
clc;
clear;

% Load parameters
load('Parameters.mat');

% Load steady state results
load('SteadyStateResults.mat');

% Generate equations
equations;

EX = false(n_equ, n_exo);
EXP = false(n_equ, n_exo);
EY = false(n_equ, n_endo);
EYP = false(n_equ, n_endo);

for i = 1:n_equ
    tmp_vars = symvar(EQU(i));
    for j = 1:numel(tmp_vars)
        tmp = char(tmp_vars(j));
        if (tmp(1) == 'y')
            if (tmp(2) == 'p')
                EYP(i, eval(tmp(3:end))) = 1;
            else
                EY(i, eval(tmp(2:end))) = 1;
            end
        elseif (tmp(1) == 'x')
            if (tmp(2) == 'p')
                EXP(i, eval(tmp(3:end))) = 1;
            else
                EX(i, eval(tmp(2:end))) = 1;
            end
        end
    end
end

% Take first order derivatives
A = sym(zeros(n_equ, n_equ));
B = sym(zeros(n_equ, n_equ));
% It seems that C is not important...
% C = sym(zeros(n_equ, n_shocks));
idxx = linspace(1, n_exo, n_exo);
idxy = linspace(n_exo + 1, n_equ, n_endo);
for i = 1:n_equ
    A(i, idxx(EXP(i, :))) = jacobian(EQU(i), XP(EXP(i, :)));
    A(i, idxy(EYP(i, :))) = jacobian(EQU(i), YP(EYP(i, :)));
    B(i, idxx(EX(i, :))) = jacobian(EQU(i), X(EX(i, :)));
    B(i, idxy(EY(i, :))) = jacobian(EQU(i), Y(EY(i, :)));
end
% C = jacobian(EQU, SHOCK);

NA = zeros(n_equ, n_equ);
NB = zeros(n_equ, n_equ);
for i = 1:n_equ
    NA(i, [idxx(EXP(i, :)), idxy(EYP(i, :))]) = double(subs(A(i, [idxx(EXP(i, :)), idxy(EYP(i, :))]), [X(EX(i, :)), XP(EXP(i, :)), Y(EY(i, :)), YP(EYP(i, :)), SHOCK], [XSS(EX(i, :)), XSS(EXP(i, :)), YSS(EY(i, :)), YSS(EYP(i, :)), SHOCKSS]));
    NB(i, [idxx(EX(i, :)), idxy(EY(i, :))]) = double(subs(B(i, [idxx(EX(i, :)), idxy(EY(i, :))]), [X(EX(i, :)), XP(EXP(i, :)), Y(EY(i, :)), YP(EYP(i, :)), SHOCK], [XSS(EX(i, :)), XSS(EXP(i, :)), YSS(EY(i, :)), YSS(EYP(i, :)), SHOCKSS]));
end

% Schur decomposition
[AA, BB, Q, Z] = qz(NA, NB, 'real');
% TODO: check robustness here
[AA, BB, Q, Z] = ordqz(AA, BB, Q, Z, 'udo');
% Check B-K condition. Number of explosive roots should equal
% forward-looking varibles
% disp(sum(abs(ordeig(BB, AA)) > 1));

% Policies
Zp = Z';
g1 = -Zp(n_exo+1:end, n_exo+1:end) \ Zp(n_exo+1:end, 1:n_exo);
h1 = -(NA(1:n_exo, 1:n_exo) + NA(1:n_exo, n_exo+1:end) * g1) \ ...
     (NB(1:n_exo, 1:n_exo) + NB(1:n_exo, n_exo+1:end) * g1);

%{
% Construct second-order equations, linear system
% Put equation dimension in the end for the sake of MATLAB
Hyp = sym(zeros(n_equ, n_endo));
Hxp = sym(zeros(n_equ, n_exo));
Hy = sym(zeros(n_equ, n_endo));
Hx = sym(zeros(n_equ, n_exo));

Hypyp = sym(zeros(n_endo, n_endo, n_equ));
Hypy = sym(zeros(n_endo, n_endo, n_equ));
Hyy = sym(zeros(n_endo, n_endo, n_equ));

Hypxp = sym(zeros(n_endo, n_exo, n_equ));
Hypx = sym(zeros(n_endo, n_exo, n_equ));
Hyxp = sym(zeros(n_endo, n_exo, n_equ));
Hyx = sym(zeros(n_endo, n_exo, n_equ));

Hxpxp = sym(zeros(n_exo, n_exo, n_equ));
Hxpx = sym(zeros(n_exo, n_exo, n_equ));
Hxx = sym(zeros(n_exo, n_exo, n_equ));

for i = 1:n_equ
    disp(i);
    % X Jacobian
    for j = 1:n_exo
        Hxp(i, j) = A(i, j);
        Hx(i, j) = B(i, j);
    end
    % Y Jacobian
    for j = 1:n_endo
        Hyp(i, j) = A(i, j + n_exo);
        Hy(i, j) = B(i, j + n_exo);
    end
    % YYs
    for j = 1:n_endo
        for k = 1:n_endo
            Hypyp(j, k, i) = diff(Hyp(i, j), YP(k));
            Hypy(j, k, i) = diff(Hyp(i, j), Y(k));
            Hyy(j, k, i) = diff(Hy(i, j), Y(k));
        end
    end
    % YXs
    for j = 1:n_endo
        for k = 1:n_exo
            Hypxp(j, k, i) = diff(Hyp(i, j), XP(k));
            Hypx(j, k, i) = diff(Hyp(i, j), X(k));
            Hyxp(j, k, i) = diff(Hy(i, j), XP(k));
            Hyx(j, k, i) = diff(Hy(i, j), X(k));
        end
    end
    % XXs
    for j = 1:n_exo
        for k = 1:n_exo
            Hxpxp(j, k, i) = diff(Hxp(i, j), XP(k));
            Hxpx(j, k, i) = diff(Hxp(i, j), X(k));
            Hxx(j, k, i) = diff(Hx(i, j), X(k));
        end
    end
end

Hxp = double(subs(Hxp, [X, XP, Y, YP], [XSS, XSS, YSS, YSS]));
Hyp = double(subs(Hyp, [X, XP, Y, YP], [XSS, XSS, YSS, YSS]));
Hy = double(subs(Hy, [X, XP, Y, YP], [XSS, XSS, YSS, YSS]));

Hypyp = double(subs(Hypyp, [X, XP, Y, YP], [XSS, XSS, YSS, YSS]));
Hypy = double(subs(Hypy, [X, XP, Y, YP], [XSS, XSS, YSS, YSS]));
Hyy = double(subs(Hyy, [X, XP, Y, YP], [XSS, XSS, YSS, YSS]));

Hypxp = double(subs(Hypxp, [X, XP, Y, YP], [XSS, XSS, YSS, YSS]));
Hypx = double(subs(Hypx, [X, XP, Y, YP], [XSS, XSS, YSS, YSS]));
Hyxp = double(subs(Hyxp, [X, XP, Y, YP], [XSS, XSS, YSS, YSS]));
Hyx = double(subs(Hyx, [X, XP, Y, YP], [XSS, XSS, YSS, YSS]));

Hxpxp = double(subs(Hxpxp, [X, XP, Y, YP], [XSS, XSS, YSS, YSS]));
Hxpx = double(subs(Hxpx, [X, XP, Y, YP], [XSS, XSS, YSS, YSS]));
Hxx = double(subs(Hxx, [X, XP, Y, YP], [XSS, XSS, YSS, YSS]));

% Tranpose all the other matrices
for i = 1:n_equ
    Hxpyp = permute(Hypxp, [2 1 3]);
    Hxpy = permute(Hyxp, [2 1 3]);
    Hxyp = permute(Hypx, [2 1 3]);
    Hxy = permute(Hyx, [2 1 3]);
    Hyyp = permute(Hypy, [2 1 3]);
    Hxxp = permute(Hxpx, [2 1 3]);
end

% equations are (h,g) coeff * UNKNOWN + const = 0
const = zeros(n_equ * n_exo * n_exo, 1);
coeff = zeros(n_equ * n_exo * n_exo, n_equ * n_exo * n_exo);

for i = 1:n_equ
    for j = 1:n_exo
        for k = 1:n_exo
            idx = (i - 1) * n_exo * n_exo + (j - 1) * n_exo + k;
            t1 = (Hypyp(:, :, i) * g1 * h1(:, k) + Hypy(:, :, i) * g1(:, k) + Hypxp(:, :, i) * h1(:, k) + Hypx(:, k, i))' * (g1 * h1(:, j));
            t2 = (Hyyp(:, :, i) * g1 * h1(:, k) + Hyy(:, :, i) * g1(:, k) + Hyxp(:, :, i) * h1(:, k) + Hyx(:, k, i))' * g1(:, j);
            t3 = (Hxpyp(:, :, i) * g1 * h1(:, k) + Hxpy(:, :, i) * g1(:, k) + Hxpxp(:, :, i) * h1(:, k) + Hxpx(:, k, i))' * h1(:, j);
            t4 = Hxyp(j, :, i) * g1 * h1(:, k) + Hxy(j, :, i) * g1(:, k) + Hxxp(j, :, i) * h1(:, k) + Hxx(j, k, i);
            const(idx, 1) = t1 + t2 + t3 + t4;

            for tmp = 1:n_endo
                for t1 = 1:n_exo
                    for t2 = 1:n_exo
                        coeff(idx, (n_exo + tmp - 1) * n_exo * n_exo + (t1 - 1) * n_exo + t2) = ...
                        coeff(idx, (n_exo + tmp - 1) * n_exo * n_exo + (t1 - 1) * n_exo + t2) + Hyp(i, tmp) * h1(t2, k) * h1(t1, j);
                    end
                end
            end

            for tmp = 1:n_endo
                coeff(idx, (n_exo + tmp - 1) * n_exo * n_exo + (j - 1) * n_exo + k) = ...
                coeff(idx, (n_exo + tmp - 1) * n_exo * n_exo + (j - 1) * n_exo + k) + Hy(i, tmp);
            end

            for tmp = 1:n_exo
                tmpres = Hyp(i, :) * g1(:, tmp);
                coeff(idx, (tmp - 1) * n_exo * n_exo + (j - 1) * n_exo + k) = tmpres + Hxp(i, tmp);
            end
        end
    end
end

res = coeff \ (- const);
for i = 1:n_exo
    for j = 1:n_exo
        for k = 1:n_exo
            h2(i, j, k) = res((i - 1) * n_exo * n_exo + (j - 1) * n_exo + k);
        end
    end
end


for i = 1:n_endo
    for j = 1:n_exo
        for k = 1:n_exo
            g2(i, j, k) = res((n_exo + i - 1) * n_exo * n_exo + (j - 1) * n_exo + k);
        end
    end
end

% Construct second-order equations for sigma
eta = zeros(n_exo, n_shocks);
eta(n_exo - n_shocks + 1:end, :) = 1; % indicate (relative) variance
I = diag(n_shocks);

% Since cross terms are 0, we only need to solve gsigma and hsigma
const_sigma = zeros(n_equ, 1);
coeff_sigma = zeros(n_equ, n_equ);

% equations are (h,g) coeff * UNKNOWN + const = 0
for i1 = 1:n_equ
    for i2 = 1:n_endo
        for i3 = 1:n_exo
            coeff_sigma(i1, i3) = coeff_sigma(i1, i3) + Hyp(i1, i2) * g1(i2, i3);
        end
    end
    for i2 = 1:n_endo
        coeff_sigma(i1, n_exo + i2) = Hyp(i1, i2) + Hy(i1, i2);
    end
    for i2 = 1:n_exo
        coeff_sigma(i1, i2) = coeff_sigma(i1, i2) + Hxp(i1, i2);
    end

    for i2 = 1:n_endo
        for i3 = 1:n_endo
            for i4 = 1:n_exo
                for i5 = 1:n_shocks
                    for i6 = 1:n_exo
                        for i7 = 1:n_shocks
                            const_sigma(i1) = const_sigma(i1) + Hypyp(i2, i3, i1) * g1(i3, i4) * eta(i4, i5) * g1(i2, i6) * eta(i6, i7) * I(i7, i5);
                        end
                    end
                end
            end
        end
    end

    for i2 = 1:n_endo
        for i4 = 1:n_exo
            for i5 = 1:n_shocks
                for i6 = 1:n_exo
                    for i7 = 1:n_shocks
                        const_sigma(i1) = const_sigma(i1) + Hypxp(i2, i4, i1) * eta(i4, i5) * g1(i2, i6) * eta(i6, i7) * I(i7, i5);
                    end
                end
            end
        end
    end

    for i2 = 1:n_endo
        for i4 = 1:n_exo
            for i5 = 1:n_exo
                for i6 = 1:n_shocks
                    for i7 = 1:n_shocks
                        const_sigma(i1) = const_sigma(i1) + Hyp(i1, i2) * g2(i2, i4, i5) * eta(i5, i6) * eta(i4, i7) * I(i7, i6);
                    end
                end
            end
        end
    end

    for i2 = 1:n_exo
        for i3 = 1:n_endo
            for i4 = 1:n_exo
                for i5 = 1:n_shocks
                    for i6 = 1:n_shocks
                        const_sigma(i1) = const_sigma(i1) + Hxpyp(i2, i3, i1) * g1(i3, i4) * eta(i4, i5) * eta(i2, i6) * I(i6, i5);
                    end
                end
            end
        end
    end

    for i2 = 1:n_exo
        for i3 = 1:n_exo
            for i5 = 1:n_shocks
                for i6 = 1:n_shocks
                    const_sigma(i1) = const_sigma(i1) + Hxpxp(i2, i3, i1) * eta(i3, i5) * eta(i2, i6) * I(i6, i5);
                end
            end
        end
    end
end

res = coeff_sigma \ (- const_sigma);
h_sigma = res(1:n_exo, 1);
g_sigma = res(n_equ - n_endo + 1:end, 1);

%save RESULTS.mat g1 g2 h1 h2 g_sigma h_sigma;
%}
