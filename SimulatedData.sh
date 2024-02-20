#!/bin/bash

########Generate and evaluate simulated data for various L values########
# Define the values of L
L_values=(100 200 300 400 500 600 700 800 900 1000)
for L in "${L_values[@]}"; do
    echo "Executing Sim_L_$L ..."

    # generate reference genome:
    echo "Generating reference genome for Sim_L_$L ..."
    python insertRepeat.py --name "Sim_L_$L" --L "$L"

    # generate the reads:
    echo "Generating reads for Sim_L_$L ..."
    wgsim -N 2000000 -1 101 -2 101 "Data/Sim_L_$L/ref_genome.fasta" "Data/Sim_L_$L/outRead1.fq" "Data/Sim_L_$L/outRead2.fq" > "Results/Sim_L_$L/wgsim.log"

    # generate the assembly graph:
    echo "Generating assembly graph for Sim_L_$L ..."
    python mainGraph.py --name "Sim_L_$L"

    # detect the repeats:
    echo "Detecting repeats for Sim_L_$L ..."
    python mainRepeatDetection.py --name "Sim_L_$L" > "Results/Sim_L_$L/repeatDetection.log"

    echo "Sim_L_$L Finished"
done

########Generate and evaluate simulated data for various C values########
C_values=(10 20 30 40 50 60 70 80 90 100 125 150)
for C in "${C_values[@]}"; do
    echo "Executing Sim_C_$C ..."

    # generate reference genome:
    echo "Generating reference genome for Sim_C_$C ..."
    python insertRepeat.py --name "Sim_C_$C" --C "$C"

    # generate the reads:
    echo "Generating reads for Sim_C_$C ..."
    wgsim -N 2000000 -1 101 -2 101 "Data/Sim_C_$C/ref_genome.fasta" "Data/Sim_C_$C/outRead1.fq" "Data/Sim_C_$C/outRead2.fq" > "Results/Sim_C_$C/wgsim.log"

    # generate the assembly graph:
    echo "Generating assembly graph for Sim_C_$C ..."
    python mainGraph.py --name "Sim_C_$C"

    # detect the repeats:
    echo "Detecting repeats for Sim_C_$C ..."
    python mainRepeatDetection.py --name "Sim_C_$C" > "Results/Sim_C_$C/repeatDetection.log"

    echo "Sim_C_$C Finished"
done

########Generate and evaluate simulated data for various coverage values########

L_value= 400
C_value= 25
coverage_values=(500000 750000 1000000 1250000 1500000 1750000 2000000 2250000 2500000)
for cov in "${coverage_values[@]}"; do
    echo "Executing Sim_Cov_$cov ..."

    # generate reference genome:
    echo "Generating reference genome for Sim_Cov_$cov ..."
    python insertRepeat.py --name "Sim_Cov_$cov" --C "$C_value" --L "$L_value"

    # generate the reads:
    echo "Generating reads for Sim_Cov_$cov ..."
    wgsim -N "$cov" -1 101 -2 101 "Data/Sim_Cov_$cov/ref_genome.fasta" "Data/Sim_Cov_$cov/outRead1.fq" "Data/Sim_Cov_$cov/outRead2.fq" > "Results/Sim_Cov_$cov/wgsim.log"

    # generate the assembly graph:
    echo "Generating assembly graph for Sim_Cov_$cov ..."
    python mainGraph.py --name "Sim_Cov_$cov"

    # detect the repeats:
    echo "Detecting repeats for Sim_Cov_$cov ..."
    python mainRepeatDetection.py --name "Sim_Cov_$cov" > "Results/Sim_Cov_$cov/repeatDetection.log"

    echo "Sim_Cov_$cov Finished"
done
