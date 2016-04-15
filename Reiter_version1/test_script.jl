using MAT
include("second_order.jl")
@time read_and_define("WARMUP.mat")
@time second_order(Hxp, Hyp, Hy, Hypyp, Hypy, Hypxp, Hypx, Hyy, Hyxp, Hyx, Hxpxp, Hxpx, Hxx, EYP, EY, EXP, EX, g1, h1, n_exo, n_endo, n_shocks, n_equ)
@time read_and_define("NECESSARY.mat");
@time h2, g2 = second_order(Hxp, Hyp, Hy, Hypyp, Hypy, Hypxp, Hypx, Hyy, Hyxp, Hyx, Hxpxp, Hxpx, Hxx, EYP, EY, EXP, EX, g1, h1, n_exo, n_endo, n_shocks, n_equ);
