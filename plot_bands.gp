# Set output
set terminal pngcairo enhanced font "Arial,12" size 1200,800
set output 'band_structure.png'

# Title and labels
set title "Silicon Band Structure"
set xlabel "k-path"
set ylabel "Energy (eV)"

# Fermi energy (replace 6.123 with your value)
Ef = 6.123
FILE_NAME="si_bands.dat.gnu"



# Plot bands (column 1: k-path, column 2: energy)
plot  FILE_NAME u 1:($2-Ef) w l lc rgb "blue" notitle

# Add Fermi energy line
set arrow from graph 0,0 to graph 1,0 nohead lc rgb "red" dt 2
set label "E_F" at graph 0.02,0.03 center





# High-symmetry points (positions from si_bands.dat.gnu)
set xtics ("Γ" 0.0, "X" 1.0, "K" 1.5, "Γ" 2.0)
# Vertical lines at high-symmetry points
set xtics add ("" 0.0, "" 1.0, "" 1.5, "" 2.0)
set grid xtics lt -1 lw 1.5 lc rgb "gray"

# Replot
plot FILE_NAME u 1:($2-Ef) w l lc rgb "blue" notitle




plot for [i=2:8] FILE_NAME u 1:(column(i)-Ef) w l lc i lw 1.5 notitle



set terminal postscript eps enhanced color
set output 'bands.pdf'
replot


