#!/bin/bash
num_iter=10

########Generate and evaluate simulated data for various Repeat Length values########
L_values=(150 200 300 400 500 600 700 800 900 1000)
for L in "${L_values[@]}"; do
    echo "Executing Sim_L_$L..."
    for ((i=1; i<=num_iter; i++)); do
        echo "iteration $i..."

        # generate reference genome:
        echo "Generating reference genome ..."
        python insertRepeat.py --name "Sim_L_$L-$i" --L "$L"
        mkdir -p "Results/"Sim_L_$L-$i""

        # generate the reads:
        echo "Generating reads ..."
        wgsim -N 2000000 -1 101 -2 101 "Data/"Sim_L_$L-$i"/ref_genome.fasta" "Data/"Sim_L_$L-$i"/outRead1.fq" "Data/"Sim_L_$L-$i"/outRead2.fq" > "Results/"Sim_L_$L-$i"/wgsim.log"

        # generate the assembly graph:
        echo "Generating assembly graph ...."
        python mainGraph.py --name ""Sim_L_$L-$i"" > "Results/"Sim_L_$L-$i"/graph.log"

        # detect the repeats:
        echo "Detecting repeats ..."
        python mainRepeatDetection.py --p 50 --name ""Sim_L_$L-$i"" > "Results/"Sim_L_$L-$i"/repeatDetection.log"

    done
    echo ""
done

########Generate and evaluate simulated data for various Copy number values########
C_values=(10 20 30 40 50 60 70 80 90 100 125 150)
for C in "${C_values[@]}"; do
    echo "Executing Sim_C_$C..."
    for ((i=1; i<=num_iter; i++)); do
        echo "iteration $i..."

        # generate reference genome:
        echo "Generating reference genome ..."
        python insertRepeat.py --name "Sim_C_$C-$i" --C "$C"
        mkdir -p "Results/Sim_C_$C-$i"

        # generate the reads:
        echo "Generating reads ..."
        wgsim -N 2000000 -1 101 -2 101 "Data/Sim_C_$C-$i/ref_genome.fasta" "Data/Sim_C_$C-$i/outRead1.fq" "Data/Sim_C_$C-$i/outRead2.fq" > "Results/Sim_C_$C-$i/wgsim.log"

        # generate the assembly graph:
        echo "Generating assembly graph ...."
        python mainGraph.py --name "Sim_C_$C-$i" > "Results/Sim_C_$C-$i/graph.log"

        # detect the repeats:
        echo "Detecting repeats ..."
        python mainRepeatDetection.py --p 50 --name "Sim_C_$C-$i" > "Results/Sim_C_$C-$i/repeatDetection.log"

    done
    echo ""
done

########Generate and evaluate simulated data for various coverage values########
L_value=400
C_value=20
coverage_values=(5 10 15 20 25 30 35 40 45 50)
for cov in "${coverage_values[@]}"; do
    echo "Executing Sim_Cov_$cov..."
    N=$((cov * 5 * 10000))
    echo "N=$N"
    for ((i=1; i<=num_iter; i++)); do
        echo "iteration $i..."

        # generate reference genome:
        echo "Generating reference genome ..."
        python insertRepeat.py --name "Sim_Cov_$cov-$i" --C "$C_value" --L "$L_value"
        mkdir -p "Results/Sim_Cov_$cov-$i"

        # generate the reads:
        echo "Generating reads ..."
        wgsim -N "$N" -1 101 -2 101 "Data/Sim_Cov_$cov-$i/ref_genome.fasta" "Data/Sim_Cov_$cov-$i/outRead1.fq" "Data/Sim_Cov_$cov-$i/outRead2.fq" > "Results/Sim_Cov_$cov-$i/wgsim.log"

        # generate the assembly graph:
        echo "Generating assembly graph ...."
        python mainGraph.py --name "Sim_Cov_$cov-$i" > "Results/Sim_Cov_$cov-$i/graph.log"

        # detect the repeats:
        echo "Detecting repeats ..."
        python mainRepeatDetection.py --p 50 --name "Sim_Cov_$cov-$i" > "Results/Sim_Cov_$cov-$i/repeatDetection.log"

    done
    echo ""
done

