# Solve second-order equations given data
function read_and_define(filename::AbstractString)
    # Read data from MATLAB
    ml = matread(filename)
    for nm in keys(ml)
        snm = symbol(nm)
        eval(current_module(), :($snm = $(ml[nm])))
    end
end

function second_order(Hxp, Hyp, Hy, Hypyp, Hypy, Hypxp, Hypx, Hyy, Hyxp, Hyx, Hxpxp, Hxpx, Hxx, EYP, EY, EXP, EX, g1, h1, n_exo, n_endo, n_shocks, n_equ)
    n_exo = Int64(n_exo)
    n_endo = Int64(n_endo)
    n_shocks = Int64(n_shocks)
    n_equ = Int64(n_equ)
    idxx = collect(1:1:n_exo)
    idxy = collect(1:1:n_endo)
    # Deal with sigular...
    if (n_exo == 1)
        h1 = [h1]
        if (n_endo == 1) g1 = [g1] end
    end
    # Tranpose all the other matrices
    Hxpyp = permutedims(Hypxp, [2 1 3])
    Hxpy = permutedims(Hyxp, [2 1 3])
    Hxyp = permutedims(Hypx, [2 1 3])
    Hxy = permutedims(Hyx, [2 1 3])
    Hyyp = permutedims(Hypy, [2 1 3])
    Hxxp = permutedims(Hxpx, [2 1 3])

    # equations are (h,g) coeff * UNKNOWN + cons = 0
    cons = zeros(n_equ * n_exo * n_exo, 1);
    coef = spzeros(n_equ * n_exo * n_exo, n_equ * n_exo * n_exo);

    for i = 1:n_equ
        println(i)
        YPIND = idxy[EYP[i, :][:]]
        YIND = idxy[EY[i, :][:]]
        XPIND = idxx[EXP[i, :][:]]
        XIND = idxx[EX[i, :][:]]
        for j = 1:n_exo
            for k = 1:n_exo
                idx = (i - 1) * n_exo * n_exo + (j - 1) * n_exo + k
                tmpgh = g1[YPIND, :] * h1[:, k]
                t1 = (Hypyp[YPIND, YPIND, i] * tmpgh + Hypy[YPIND, YIND, i] * g1[YIND, k] + Hypxp[YPIND, XPIND, i] * h1[XPIND, k] + Hypx[YPIND, k, i])' * (g1[YPIND, :] * h1[:, j])
                t2 = (Hyyp[YIND, YPIND, i] * tmpgh + Hyy[YIND, YIND, i] * g1[YIND, k] + Hyxp[YIND, XPIND, i] * h1[XPIND, k] + Hyx[YIND, k, i])' * g1[YIND, j]
                t3 = (Hxpyp[XPIND, YPIND, i] * tmpgh + Hxpy[XPIND, YIND, i] * g1[YIND, k] + Hxpxp[XPIND, XPIND, i] * h1[XPIND, k] + Hxpx[XPIND, k, i])' * h1[XPIND, j]
                t4 = Hxyp[j, YPIND, i] * tmpgh + Hxy[j, YIND, i] * g1[YIND, k] + Hxxp[j, XPIND, i] * h1[XPIND, k] + Hxx[j, k, i]
                cons[idx] = -(t1 + t2 + t3 + t4)[1]
                for tmp = YPIND
                    coef[idx, (n_exo + tmp - 1) * n_exo * n_exo + 1:(n_exo + tmp) * n_exo * n_exo] = coef[idx, (n_exo + tmp - 1) * n_exo * n_exo + 1:(n_exo + tmp) * n_exo * n_exo] + Hyp[i, tmp] * reshape(h1[:, k] * h1[:, j]', 1, n_exo * n_exo)
                end
                for tmp = YIND
                    coef[idx, (n_exo + tmp - 1) * n_exo * n_exo + (j - 1) * n_exo + k] = coef[idx, (n_exo + tmp - 1) * n_exo * n_exo + (j - 1) * n_exo + k] + Hy[i, tmp];
                end
                for tmp = 1:n_exo
                    coef[idx, (tmp - 1) * n_exo * n_exo + (j - 1) * n_exo + k] = (Hyp[i, YPIND] * g1[YPIND, tmp] + Hxp[i, tmp])[1];
                end
            end
        end
    end

    res = coef \ cons
    h2 = reshape(res[1:n_exo * n_exo * n_exo, 1], n_exo, n_exo, n_exo)
    g2 = reshape(res[n_exo * n_exo * n_exo + 1:end, 1], n_exo, n_exo, n_endo)
    h2 = permutedims(h2, [3 2 1])
    g2 = permutedims(g2, [3 2 1])

    return h2, g2
end
