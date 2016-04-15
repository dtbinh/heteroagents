% Solve the dynamic model
clc;
clear;

% Load parameters
load('Parameters.mat');

% Load steady state results
load('SteadyStateResults.mat');

% Generate equations
equations;

% Record all the variables used
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
    % disp(i);
    A(i, idxx(EXP(i, :))) = jacobian(EQU(i), XP(EXP(i, :)));
    A(i, idxy(EYP(i, :))) = jacobian(EQU(i), YP(EYP(i, :)));
    B(i, idxx(EX(i, :))) = jacobian(EQU(i), X(EX(i, :)));
    B(i, idxy(EY(i, :))) = jacobian(EQU(i), Y(EY(i, :)));
end
% C = jacobian(EQU, SHOCK);

NA = zeros(n_equ, n_equ);
NB = zeros(n_equ, n_equ);
for i = 1:n_equ
    % disp(i);
    NA(i, [idxx(EXP(i, :)), idxy(EYP(i, :))]) = double(subs(A(i, [idxx(EXP(i, :)), idxy(EYP(i, :))]), [X(EX(i, :)), XP(EXP(i, :)), Y(EY(i, :)), YP(EYP(i, :)), SHOCK], [XSS(EX(i, :)), XSS(EXP(i, :)), YSS(EY(i, :)), YSS(EYP(i, :)), SHOCKSS]));
    NB(i, [idxx(EX(i, :)), idxy(EY(i, :))]) = double(subs(B(i, [idxx(EX(i, :)), idxy(EY(i, :))]), [X(EX(i, :)), XP(EXP(i, :)), Y(EY(i, :)), YP(EYP(i, :)), SHOCK], [XSS(EX(i, :)), XSS(EXP(i, :)), YSS(EY(i, :)), YSS(EYP(i, :)), SHOCKSS]));
end

% Schur decomposition, Q'*AA*Z'=NA, Q'*BB*Z'=NB
[AA, BB, Q, Z] = qz(NA, NB);
% TODO: maybe sort here to satisfy B-K?
sel = (abs(diag(AA)) > (1 - crit.eps) .* abs(diag(BB)));
[AA, BB, Q, Z] = ordqz(AA, BB, Q, Z, sel);
% Check B-K condition. Number of explosive roots should equal
% forward-looking varibles
% disp(n_exo - sum(sel));

% Policies
Zp = Z';
g1 = -Zp(n_exo+1:end, n_exo+1:end) \ Zp(n_exo+1:end, 1:n_exo);
% disp(max(max(abs(imag(g1)))));
g1 = real(g1);
h1 = pinv(-(NA(1:n_exo, 1:n_exo) + NA(1:n_exo, n_exo+1:end) * g1)) * ...
     (NB(1:n_exo, 1:n_exo) + NB(1:n_exo, n_exo+1:end) * g1);
disp('First order done');
 
Hyp = sym(zeros(n_equ, n_endo));
Hxp = sym(zeros(n_equ, n_exo));
Hy = sym(zeros(n_equ, n_endo));
Hx = sym(zeros(n_equ, n_exo));

Hypyp = zeros(n_endo, n_endo, n_equ);
Hypy = zeros(n_endo, n_endo, n_equ);
Hyy = zeros(n_endo, n_endo, n_equ);

Hypxp = zeros(n_endo, n_exo, n_equ);
Hypx = zeros(n_endo, n_exo, n_equ);
Hyxp = zeros(n_endo, n_exo, n_equ);
Hyx = zeros(n_endo, n_exo, n_equ);

Hxpxp = zeros(n_exo, n_exo, n_equ);
Hxpx = zeros(n_exo, n_exo, n_equ);
Hxx = zeros(n_exo, n_exo, n_equ);

% Generate new index...
old_idxy = idxy;
idxy = linspace(1, n_endo, n_endo);
% For suck MATLAB high order matrix
anyYP = any(EYP, 2);
anyY = any(EY, 2);
anyXP = any(EXP, 2);
anyX = any(EX, 2);
for i = 1:n_equ
    disp(i);
    tmp_vars = [X(EX(i, :)), XP(EXP(i, :)), Y(EY(i, :)), YP(EYP(i, :)), SHOCK];
    tmp_ss = [XSS(EX(i, :)), XSS(EXP(i, :)), YSS(EY(i, :)), YSS(EYP(i, :)), SHOCKSS];
    % X Jacobian
    Hxp(i, idxx(EXP(i, :))) = A(i, idxx(EXP(i, :)));
    Hx(i, idxx(EX(i, :))) = B(i, idxx(EX(i, :)));
    % Y Jacobian
    Hyp(i, idxy(EYP(i, :))) = A(i, old_idxy(EYP(i, :)));
    Hy(i, idxy(EY(i, :))) = B(i, old_idxy(EY(i, :)));
    if anyYP(i)
        Hypyp(idxy(EYP(i, :)), idxy(EYP(i, :)), i) = double(subs(jacobian(Hyp(i, idxy(EYP(i, :))), YP(EYP(i, :))), tmp_vars, tmp_ss));
        if anyY(i) Hypy(idxy(EYP(i, :)), idxy(EY(i, :)), i) = double(subs(jacobian(Hyp(i, idxy(EYP(i, :))), Y(EY(i, :))), tmp_vars, tmp_ss)); end
        if anyXP(i) Hypxp(idxy(EYP(i, :)), idxx(EXP(i, :)), i) = double(subs(jacobian(Hyp(i, idxy(EYP(i, :))), XP(EXP(i, :))), tmp_vars, tmp_ss)); end
        if anyX(i) Hypx(idxy(EYP(i, :)), idxx(EX(i, :)), i) = double(subs(jacobian(Hyp(i, idxy(EYP(i, :))), X(EX(i, :))), tmp_vars, tmp_ss)); end
    end
    if anyY(i)
        Hyy(idxy(EY(i, :)), idxy(EY(i, :)), i) = double(subs(jacobian(Hy(i, idxy(EY(i, :))), Y(EY(i, :))), tmp_vars, tmp_ss));
        if anyXP(i) Hyxp(idxy(EY(i, :)), idxx(EXP(i, :)), i) = double(subs(jacobian(Hy(i, idxy(EY(i, :))), XP(EXP(i, :))), tmp_vars, tmp_ss)); end
        if anyX(i) Hyx(idxy(EY(i, :)), idxx(EX(i, :)), i) = double(subs(jacobian(Hy(i, idxy(EY(i, :))), X(EX(i, :))), tmp_vars, tmp_ss)); end
    end
    if anyXP(i)
        Hxpxp(idxx(EXP(i, :)), idxx(EXP(i, :)), i) = double(subs(jacobian(Hxp(i, idxx(EXP(i, :))), XP(EXP(i, :))), tmp_vars, tmp_ss));
        if anyX(i) Hxpx(idxx(EXP(i, :)), idxx(EX(i, :)), i) = double(subs(jacobian(Hxp(i, idxx(EXP(i, :))), X(EX(i, :))), tmp_vars, tmp_ss)); end
    end
    if anyX(i)
        Hxx(idxx(EX(i, :)), idxx(EX(i, :)), i) = double(subs(jacobian(Hx(i, idxx(EX(i, :))), X(EX(i, :))), tmp_vars, tmp_ss));
    end
    Hxp(i, idxx(EXP(i, :))) = subs(Hxp(i, idxx(EXP(i, :))), tmp_vars, tmp_ss);
    Hyp(i, idxy(EYP(i, :))) = subs(Hyp(i, idxy(EYP(i, :))), tmp_vars, tmp_ss);
    Hy(i, idxy(EY(i, :))) = subs(Hy(i, idxy(EY(i, :))), tmp_vars, tmp_ss);
end
clear old_idxy;

Hxp = double(Hxp); Hyp = double(Hyp); Hy = double(Hy);

% Tranpose all the other matrices
Hxpyp = permute(Hypxp, [2 1 3]);
Hxpy = permute(Hyxp, [2 1 3]);
Hxyp = permute(Hypx, [2 1 3]);
Hxy = permute(Hyx, [2 1 3]);
Hyyp = permute(Hypy, [2 1 3]);
Hxxp = permute(Hxpx, [2 1 3]);

% equations are (h,g) coeff * UNKNOWN + const = 0
const = zeros(n_equ * n_exo * n_exo, 1);
coeff = sparse(n_equ * n_exo * n_exo, n_equ * n_exo * n_exo);

for i = 1:n_equ
    disp(i);
    YPIND = idxy(EYP(i, :));
    YIND = idxy(EY(i, :));
    XPIND = idxx(EXP(i, :));
    XIND = idxx(EX(i, :));
    for j = 1:n_exo
        for k = 1:n_exo
            idx = (i - 1) * n_exo * n_exo + (j - 1) * n_exo + k;
            tmpgh = g1(YPIND, :) * h1(:, k);
            t1 = (Hypyp(YPIND, YPIND, i) * tmpgh + Hypy(YPIND, YIND, i) * g1(YIND, k) + Hypxp(YPIND, XPIND, i) * h1(XPIND, k) + Hypx(YPIND, k, i))' * (g1(YPIND, :) * h1(:, j));
            t2 = (Hyyp(YIND, YPIND, i) * tmpgh + Hyy(YIND, YIND, i) * g1(YIND, k) + Hyxp(YIND, XPIND, i) * h1(XPIND, k) + Hyx(YIND, k, i))' * g1(YIND, j);
            t3 = (Hxpyp(XPIND, YPIND, i) * tmpgh + Hxpy(XPIND, YIND, i) * g1(YIND, k) + Hxpxp(XPIND, XPIND, i) * h1(XPIND, k) + Hxpx(XPIND, k, i))' * h1(XPIND, j);
            t4 = Hxyp(j, YPIND, i) * tmpgh + Hxy(j, YIND, i) * g1(YIND, k) + Hxxp(j, XPIND, i) * h1(XPIND, k) + Hxx(j, k, i);
            const(idx) = t1 + t2 + t3 + t4;

            for tmp = idxy(EYP(i, :))
                coeff(idx, (n_exo + tmp - 1) * n_exo * n_exo + 1:(n_exo + tmp) * n_exo * n_exo) = ...
                coeff(idx, (n_exo + tmp - 1) * n_exo * n_exo + 1:(n_exo + tmp) * n_exo * n_exo) + Hyp(i, tmp) * reshape(h1(:, k) * h1(:, j)', 1, n_exo * n_exo);
            end

            for tmp = idxy(EY(i, :))
                coeff(idx, (n_exo + tmp - 1) * n_exo * n_exo + (j - 1) * n_exo + k) = ...
                coeff(idx, (n_exo + tmp - 1) * n_exo * n_exo + (j - 1) * n_exo + k) + Hy(i, tmp);
            end

            for tmp = 1:n_exo
                tmpres = Hyp(i, YPIND) * g1(YPIND, tmp);
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
