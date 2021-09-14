set output "1.eps"
set terminal postscript eps enhanced solid color "Helvetica" 20
#set terminal postscript eps enhanced dl 2.0 lw 2.0 size 3.25,2.275 "Helvetica" 14

set multiplot

set size 0.55, 0.55
set origin 0.4, 0.33

set xlabel 'Shaping time [{/Symbol m}s]'
set xrange [0:500]
set ylabel 'ENC [e^-]'
set yrange [26:32]
set ytics 2

plot '600umCEE_1000e_1m_Noise_Shaping.dat' u ($1*2.0*0.1):($3*1000.0*1000.0/36.9):($3*1000.0*1000.0/36.9/sqrt(10000)) w yerrorline lc 0 t ''

set size 1,1
set origin 0,0

set key at first 4.2,1.049

set xlabel 't [ms]'
set xrange [0:6.5]
set ylabel '[V]'
set format y "%5.3f"
set yrange [0.993:1.050]
set ytics auto

plot '600umCEE_1000e_1m_wave.dat' every 10 u ($1*1000.0):2 w l lc 0 lw 0.1 t 'CSA output waveform, RC=1ms', \
     '' every 10 u ($1*1000.0):3 w step lt 1 lw 2.0 t 'After shaping (200{/Symbol m}s)'

unset multiplot

unset output
set term x11
