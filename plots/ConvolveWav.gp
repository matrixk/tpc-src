dfname="ConvolveWav.dat"

set output "1.eps"
set terminal postscript eps enhanced solid color "Helvetica" 20
#set terminal postscript eps enhanced dl 2.0 lw 2.0 size 3.25,2.275 "Helvetica" 14


set xlabel 't [ms]'
set xrange [0.63:0.88]
set ylabel '[V]'
set yrange [0.995:1.055]
set y2label 'Input number of charges per 1{/Symbol m}s slice' textcolor lt 1
set y2range [-5:55]
set y2tics

plot dfname u ($0*0.1*1e-3):($2+1.0) w l lc 0 lw 0.1 t 'CSA output waveform, RC=1ms', \
     '' u (($0-50)*0.1*1e-3):($3*1.0+1.0) w l lc 0 lw 4.0 t 'After shaping (10{/Symbol m}s)', \
     '' u (($0+605.0)*1e-3):1 axes x1y2 w step lt 1 t 'Charge input signal'

unset output
set term x11
