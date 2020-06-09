set terminal png
set output "RisultatiLTBa.png"
set autoscale
set title "Scale factors evolution"
unset key
set xlabel " t/t0"
set ylabel " a(t/t0)"
plot "RisultatiLTBMatterCostant.txt" using 1:2 with lines title "a(t)"
