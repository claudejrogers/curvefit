set terminal png
set output "examples/ic50_ex_plot.png"
f(x) = 1 - (0.817390/(1 + (0.185800/x)**1.682300))
set xrange [0.001563:102.400002]
set yrange [-0.017317:1.199684]
set log x
plot "examples/ic50_ex_input.dat" using 1:2 with points 4, f(x) with lines 22
