set terminal png
set output "RisultatiLTBb.png"
set autoscale
set title "Scale factors evolution"
unset key
set xlabel " t/t0"
set ylabel " b(t/t0)"
plot "RisultatiLTBMatterCostant.txt" using 1:3 with lines title "b(t)"
