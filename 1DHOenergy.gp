#Gnuplot script file for plotting the excitation energy of a 1D HO as a function of time
set term pdfcairo
set output "1DHO_energy.pdf"
set size 1, 1
set autoscale
unset log
unset label
set xtic auto
set ytic auto
set title "Energy vs Time for x(t)x(t) sink-source"
set ylabel "Energy"
set xlabel "Time" 
set dummy t
E(t) = a
a = 1
plot "harmonic_oscillator_metropolis.txt" using 1:3:4 title "numerically retrieved" pt 7 ps 0.25 with yerrorbars, E(t) title "asymptotic result" with lines;
set output
