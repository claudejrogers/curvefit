set terminal png
set output "examples/gaussian_ex_plot.png"
f(x) = -0.000561 + (1.002816*exp(-((x - 50.004804)**2)/(10.014453**2)))
set xrange [20.000000:79.900002]
set yrange [-0.256838:1.247147]
plot "examples/gaussian_ex_input.dat" using 1:2 with points 4, f(x) with lines 22
