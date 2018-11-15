#!/bin/bash

#executar como ./testes.sh <core> <grupo>

cache="L2 miss ratio"
cachebw="L3 bandwidth"
flopsDP="DP MFLOP/s"
flopsAVX="Packed DP MFLOP/s"

#shift 2

#${LIKWID_CMD} $*

#modprobe msr

export PATH=/home/soft/likwid/bin:/home/soft/likwid/sbin:$PATH
export LD_LIBRARY_PATH=/home/soft/likwid/lib:$LD_LIBRARY_PATH

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
for grupo in L2CACHE L3 FLOPS_DP FLOPS_AVX
do
	LIKWID_CMD="likwid-perfctr -C $1 -g $grupo -m -O"
	echo "${LIKWID_CMD}"
	#---Confere os indicadores e salva o padrao que procuramos em cada um---
	if [ "$grupo" == "L2CACHE" ]
	then
		padrao="$cache"
	elif [ "$grupo" == "L3" ]
	then
		padrao="$cachebw"
	elif [ "$grupo" == "FLOPS_DP" ]
	then
		padrao="$flopsDP"
	else
		padrao="$flopsAVX"
	fi
	#----------------------------V1-------------------------------------

	for i in 32 64 128 256 512 1000 2000 4000 8000
	do
		${LIKWID_CMD} ./trab2/cgSolver -n $i -k 7 -p 0.5 -i 10 -e -o teste$i.txt > temp.tmp
		#salvando tempo
		if [ "$grupo" == "L2CACHE" ]
		then
			printf "$(($i*8)) "
			printf "$(grep "Tempo iter" teste$i.txt | awk -F"<" '{print $2}' | tr ">\n" " ")"
			printf "$(grep "Tempo residuo" teste$i.txt | awk -F"<" '{print $2}' | tr ">\n" " ")\n"  
		fi >> tempoV1.tmp
		#salvando likwid
		printf "$(($i*8)) "
		printf "$(grep "$padrao" temp.tmp | awk -F"," '{print $2}' | tr "\n" " ")\n"

	done > V1$grupo.tmp

	#----------------------------V2-------------------------------------
	for i in 32 64 128 256 512 1000 2000 4000 8000
	do
		${LIKWID_CMD} ./cgSolver -n $i -k 7 -p 0.5 -i 10 -e -o teste$i.txt > temp.tmp
		#salvando tempo
		if [ "$grupo" == "L2CACHE" ]
		then
			printf "$(($i*8)) "
			printf "$(grep "Tempo iter" teste$i.txt | awk -F"<" '{print $2}' | tr ">\n" " ")"
			printf "$(grep "Tempo residuo" teste$i.txt | awk -F"<" '{print $2}' | tr ">\n" " ")\n"  
		fi >> tempoV2.tmp
		#salvando likwid
		printf "$(($i*8)) "
		printf "$(grep "$padrao" temp.tmp | awk -F"," '{print $2}' | tr "\n" " ")\n"
	
	done > V2$grupo.tmp

	#-----------------------------------------------------------------------
	#----------------gnuplot para fazer os graficos do likwid---------------
	gnuplot <<- EOF
		set logscale x 2
		set xlabel "N (bytes)"
		set ylabel "$escolha"
		set title "Medição de performance para $grupo como: $padrao"   
		set style data point
		set style function line
	
		set style line 1 lc rgb "red" lw 2
		set style line 2 lc rgb "orange" lw 2
		set style line 3 lc rgb "green" lw 2
		set style line 4 lc rgb "blue" lw 2
	
		set term png
		set output "$grupo.png"
		plot "V1$grupo.tmp" using 1:2 ls 1 title 'OP1-V1' with lines, \
		"V1$grupo.tmp" using 1:3 ls 2 title 'OP2-V1' with lines, \
		"V2$grupo.tmp" using 1:2 ls 3 title 'OP1-V2' with lines, \
		"V2$grupo.tmp" using 1:3 ls 4 title 'OP2-V2' with lines
	EOF
done 
#----------------gnuplot para fazer os graficos do tempo medio-----------
gnuplot <<- EOF
	set logscale x 2
	set xlabel "N (bytes)"
	set ylabel "tempo (segundos)"
	set title "Medição de tempo gasto em segundos"   
	set style data point
	set style function line

	set style line 1 lc rgb "red" lw 2
	set style line 2 lc rgb "orange" lw 2
	set style line 3 lc rgb "green" lw 2
	set style line 4 lc rgb "blue" lw 2

	set term png
	set output "tempo.png"
	plot "tempoV1.tmp" using 1:2 ls 1 title 'OP1-V1' with lines, \
	"tempoV1.tmp" using 1:3 ls 2 title 'OP2-V1' with lines, \
	"tempoV2.tmp" using 1:2 ls 3 title 'OP1-V2' with lines, \
	"tempoV2.tmp" using 1:3 ls 4 title 'OP2-V2' with lines
EOF

mkdir graficos
mv *.png graficos

rm *.tmp
