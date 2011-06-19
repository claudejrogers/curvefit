set terminal png
set output "examples/mm_ex_plot.png"
f(x) = (10.265398 * x) / (0.336936 + x)
set xrange [0.025158:2.000000]
set yrange [-1.219506:10.774654]
plot "examples/mm_ex_input.dat" using 1:2 with points 4, f(x) with lines 22
