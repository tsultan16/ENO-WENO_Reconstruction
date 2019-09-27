set terminal png
set output 'enoReconstruction.png'

set key font ",10"

set title "Piecewise Polynomial Reconstruction"
set xlabel "x"
set ylabel "v(x)"
set pointsize 0.2 
#set yrange [-0.01:1.2]
plot "output.txt" using 1:2 with linespoints pointtype 6 lc rgb "blue" title "reconstructed" ,\
"output.txt" using 1:3 with lines lc rgb "red" title "exact"

#plot "output.txt" using 1:2 with lines lc rgb "blue" title "reconstructed" ,\
#"output.txt" using 1:3 with lines lc rgb "red" title "exact"



unset output
