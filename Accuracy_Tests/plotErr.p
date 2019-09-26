set terminal png
set output 'L1.png'

set key font ",10"
#np=1000

set title "Piecewise Polynomial Reconstruction"
set xlabel "log(dx)"
set ylabel "log(L1 error)"
plot "output.txt" using 1:2 with lines lc rgb "blue" notitle ,\
"output.txt" using 1:3 with points pointtype 7 lc rgb "red" notitle

unset output
