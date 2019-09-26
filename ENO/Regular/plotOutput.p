set terminal png
set output 'enoReconstruction.png'

set key font ",10"
#np=1000

set title "Piecewise Polynomial Reconstruction"
set xlabel "x"
set ylabel "v(x)"
set yrange [-0.01:1.2]
plot "output.txt" using 1:2 with lines lc rgb "blue" title "reconstructed" ,\
"output.txt" using 1:3 with lines lc rgb "red" title "exact"



unset output
