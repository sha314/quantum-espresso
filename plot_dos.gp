set terminal pngcairo size 1000,700 enhanced font 'Arial,14'
set output 'dos.png'

set title "Density of States"
set xlabel "Energy (eV)"
set ylabel "DOS"
set grid
set key outside
set zeroaxis

set yrange [0:4];
set xrange [-5:5];

Ef = 17.7746   # adjust this to your Fermi energy
FILE_NAME="dos.dat"

plot FILE_NAME using ($1-Ef):2 with lines lw 2 lc rgb "red" title "DOS"

