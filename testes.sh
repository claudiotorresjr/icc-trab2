#!/bin/bash

#executar como ./testes.sh <core> <grupo>

LIKWID_CMD="likwid-perfctr -C $1 -g $2 -m -O"

cache="L2 miss ratio"
cachebw="L3 bandwidth"
flopsDP="Scalar MUOPS/s"
flopsAVX="Packed DP MFLOP/s"

passo="$3"
inicio="$4"
fim="$5"

#shift 2

#${LIKWID_CMD} $*

#modprobe msr

#export PATH=/home/soft/likwid/bin:/home/soft/likwid/sbin:$PATH
#export LD_LIBRARY_PATH=/home/soft/likwid/lib:$LD_LIBRARY_PATH

# Para obter lista de grupos de indicadores de performance:
#      likwid-perfctr -a

# Para obter topologia dos processadores
#      likwid-topology -c -g

#---Confere se a opçao é FLOPS_AVX. Se for, complia com a opçao -o3---
#if [ "$2" == "FLOPS_AVX" ]
#then
#	make avx
#else
#	make
#fi

#---Confere os indicadores e salva o padrao que procuramos em cada um---
if [ "$2" == "L2CACHE" ]
then
	padrao="$cache"
elif [ "$2" == "L3" ]
then
	padrao="$cachebw"
elif [ "$2" == "FLOPS_DP" ]
then
	padrao="$flopsDP"
else
	padrao="$flopsAVX"
fi
#----------------------------V1-------------------------------------

for i in 32 64 128 256 512 1000 2000 4000 8000
do
	${LIKWID_CMD} ./trab2/cgSolver -n $i -k 7 -p 0.5 -i 10 -e -o teste$i.txt > temp.tmp

	printf "$(($i*8)) "
	printf "$(grep "$padrao" temp.tmp | awk -F"," '{print $2}' | tr "\n" " ")\n"

done > V1$2.tmp


#----------------------------V2-------------------------------------
for i in 32 64 128 256 512 1000 2000 4000 8000
do
	${LIKWID_CMD} ./gSolver -n $i -k 7 -p 0.5 -i 10 -e -o teste$i.txt > temp.tmp

	printf "$(($i*8)) "
	printf "$(grep "$padrao" temp.tmp | awk -F"," '{print $2}' | tr "\n" " ")\n"

done > V2$2.tmp

#-----------------------------------------------------------------------
#----------------gnuplot para fazer os graficos-------------------------
gnuplot <<- EOF
	set logscale x 2
	set xlabel "N (bytes)"
	set ylabel "$escolha"
	set title "Medição de performance para $2"   
	set style data point
	set style function line

	set style line 1 lc rgb "red" lw 2
	set style line 2 lc rgb "orange" lw 2
	set style line 3 lc rgb "green" lw 2
	set style line 4 lc rgb "blue" lw 2

	set term png
	set output "$2.png"
	plot "V1$2.tmp" using 1:2 ls 1 title 'OP1-V1' with lines, \
	"V1$2.tmp" using 1:3 ls 2 title 'OP2-V1' with lines, \
	"V2$2.tmp" using 1:2 ls 3 title 'OP1-V2' with lines, \
	"V2$2.tmp" using 1:3 ls 4 title 'OP2-V2' with lines
EOF

rm *.tmp