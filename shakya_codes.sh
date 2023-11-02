#!/bin/bash

'''shakya_1'''
# Generateing the error-free reads, assembling them, mapping them back to the unitigs, and finding the ground truth repeats: (Step 1)
echo "Running sequencing.py..."
python sequencing.py --name shakya_1 --N 40 --readError 0 --np 8 --level 1
# Constructin unitig graph (Step 2)
echo "Running unitigGraphFeatures.py..."
python unitigGraphFeatures.py --name shakya_1 --idy 100 --nL 1 --cN 2
# Detecting the repeats (Steps 3, 4, and 5)
echo "Running repeatDetection.py..."
python repeatDetection.py --name shakya_1 --p 35 --N_iter 10

'''shakya_2'''
# Generateing the error-free reads, assembling them, mapping them back to the unitigs, and finding the ground truth repeats: (Step 1)
echo "Running sequencing.py..."
python sequencing.py --name shakya_2 --np 8 --level 2
# Constructin unitig graph (Step 2)
echo "Running unitigGraphFeatures.py..."
python unitigGraphFeatures.py --name shakya_2 --idy 95 --nL 0.95 --cN 2
# Detecting the repeats (Steps 3, 4, and 5)
echo "Running repeatDetection.py..."
python repeatDetection.py --name shakya_2 --p 35 --N_iter 10
