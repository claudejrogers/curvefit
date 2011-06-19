set terminal png
set output "examples/hill_ex_plot.png"
f(x) = (0.826000/(1 + (0.123000/x)**1.677999))
set xrange [0.001563:102.400002]
set yrange [-0.164655:0.991188]
set log x
plot "examples/hill_ex_input.dat" using 1:2 with points 4, f(x) with lines 22
