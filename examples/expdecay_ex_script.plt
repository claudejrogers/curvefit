set terminal png
set output "examples/expdecay_ex_plot.png"
f(x) = 1.011894 + (5.009268*exp(-0.100549*x))
set xrange [0.000000:39.000000]
set yrange [-0.293783:7.396519]
plot "examples/expdecay_ex_input.dat" using 1:2 with points 4, f(x) with lines 22
