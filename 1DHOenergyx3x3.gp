#Gnuplot script file for plotting the excitation energy of a 1D HO as a function of time
set term pdfcairo
set output "1DHO_energy_x3x3.pdf"
set size 1, 1
set xrange [0:3.5]
set yrange [0:2.5]
unset log
unset label
set xtic auto
set ytic auto
set title "Energy vs Time for x^3(t)x^3(t) sink-source"
set ylabel "Energy"
set xlabel "Time" 
set dummy t
E(t) = a
a = 1
plot "harmonic_oscillator_metropolis_x3x3.txt" every ::0::5  using 1:3:4 title "numerically retrieved" pt 7 ps 0.25 with yerrorbars, E(t) title "asymptotic result" with lines;
set output
