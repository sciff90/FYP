# Gnuplot script for plotting 2d random walk data
set title "2D Random Walk"
set xlabel "X"
set ylabel "Y"
set grid
set terminal postscript color
set output "plots/random_walk_2D.ps"
plot "data/random_walk.dat" with lines
