set terminal png
set output "examples/modsin_ex_plot.png"
f(x) = 2.345573*sin(pi*(x - 0.220712)/1.926802)
set xrange [0.000000:19.900000]
set yrange [-3.289913:2.998643]
plot "examples/modsin_ex_input.dat" using 1:2 with points 4, f(x) with lines 22
