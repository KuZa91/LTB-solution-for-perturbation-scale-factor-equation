set terminal png
set output "RisultatiLTBeta.png"
set autoscale
set title "Scale factors evolution"
unset key
set xlabel " t/t0"
set ylabel " eta(t/t0)"
plot "RisultatiLTBMatterCostant.txt" using 1:4 with lines title "eta(t)"
