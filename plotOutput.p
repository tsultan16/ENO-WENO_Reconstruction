set terminal png
set output 'enoReconstruction.png'

set key font ",10"
nx=50

set title "Piecewise Cubic Reconstruction"
set xlabel "x"
set ylabel "v(x)"
set yrange [-1.5:1.5]
plot "output.txt" using 1:2 with points pointtype 6 lc rgb "blue" title "reconstructed" ,\
"output.txt" using 1:3 with lines lc rgb "red" title "exact"



unset output
