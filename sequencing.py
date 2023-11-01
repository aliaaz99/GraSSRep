import argparse as ap
import subprocess
import os
import json

# parameters
parser = ap.ArgumentParser()
parser.add_argument('--name', help='Folder containing reads and reference genomes (in Data folder)', default="shakya_1")
parser.add_argument('--N', type=float, help='Number of read pairs in million', default=2)
parser.add_argument('--readError', type=float, help='Error in read generating', default=1.0)
parser.add_argument('--np', type=int, help='Number of cores inassembly and read mapping', default=8)
parser.add_argument('--level', type=int, help='Start evel of the assembly: 1 for wgsim(read generating) 2 for abyss (assembly) 3 for bowtie2(readmapping) 4 for nucmer(groundtruth repeats)', default=1)
args = parser.parse_args()

os.chdir('Data/' + args.name)
# Constant values:
readLength = 101
read_outerDist = 500
read_outerDist_std = 50
kmer = 51
max_map = 100

# name of input files:
genome_name = 'ref_genome.fasta'
readFile_name_1 = 'outRead1.fq'
readFile_name_2 = 'outRead2.fq'

# name of output files:
assemblies_name = 'assemblies'
unitig_name = assemblies_name + '-unitigs.fa'
index_name = 'unitigs_index'
myUnitig = 'unitigSeq.fa'
nucmer_name = 'unitig_mapping'
nucmer_repeat = nucmer_name + '.coords'
alignment1_name = 'alignment1.sam'
alignment2_name = 'alignment2.sam'

# name of the terminal files:
args_file_name_assembly = 'assembly_parameters.txt'
wgsim_terminal_file = "terminal_wgsim.txt"
abyss_terminal_file = "terminal_abyss.txt"
bowtie_terminal_index = "terminal_bowtie2_index.txt"
bowtie_terminal_1 = "terminal_bowtie2_read1.txt"
bowtie_terminal_2 = "terminal_bowtie2_read2.txt"
bowtie_terminal_Seq = "terminal_bowtie2_Seq.txt"
nucmer_terminal_file1= "terminal_nucmer1.txt"
nucmer_terminal_file2= "terminal_nucmer2.txt"


# check if genome_name is in the folder
if not os.path.isfile(genome_name):
    print("Genome file not found")
    exit()

# save assembly parameters:
with open(args_file_name_assembly, 'w+') as f:
    json.dump(args.__dict__, f, indent=2)
    json.dump({'readLength': readLength, 'read_outerDist': read_outerDist, 'read_outerDist_std': read_outerDist_std, 'kmer': kmer, 'max_map': max_map}, f, indent=2)


if args.level <= 1: # Wgsim (Read simulation): 
    if args.readError == 1.0: # generate reads with default values for error:
        command_wgsim = ['wgsim', '-N', str(1e6 * args.N), '-1', str(readLength), '-2', str(readLength), '-d',
                    str(read_outerDist), '-s', str(read_outerDist_std), '-S', str(7), genome_name, readFile_name_1, readFile_name_2]
    elif args.readError == 0.0: # generate reads with no error:
        command_wgsim = ['wgsim', '-N', str(1e6 * args.N), '-1', str(readLength), '-2', str(readLength), '-d',
                    str(read_outerDist), '-s', str(read_outerDist_std), '-e', str(0.0), '-r', str(0.0), '-R', str(0.0), '-X', str(0.0),
                    '-S', str(7), genome_name, readFile_name_1, readFile_name_2]
    else:
        print("readError value not correct")
    with open(wgsim_terminal_file, "w") as outfile:
        p_wgsim = subprocess.run(command_wgsim, stdout=outfile, stderr=outfile)

if args.level <= 2: # Abyss (Assembly):
    command_abyss = ['abyss-pe' + ' j='+ str(args.np) + ' np=' + str(args.np) + ' k=' + str(
        kmer) + ' name=' + assemblies_name + ' in=\'' + readFile_name_1 + ' ' + readFile_name_2 + '\'' + ' unitigs']
    with open(abyss_terminal_file, "w") as outfile:
        p_abyss = subprocess.run(command_abyss, shell=True, stdout=outfile, stderr=outfile)



if args.level <= 3: # Bowtie2 (Read mapping):
    trimValue = readLength - kmer

    command_bowtie_0 = 'bowtie2-build ' + unitig_name + ' ' + index_name

    command_bowtie_1 = 'bowtie2 ' + '-q ' + '-x ' + index_name + ' --seed ' + str(7) + ' --sensitive ' + '-k ' + str(max_map) + ' -3 ' + str(
        trimValue) + ' -U ' + readFile_name_1 + ' -p ' + str(args.np) + ' --no-unal ' + '> ' + alignment1_name

    command_bowtie_2 = 'bowtie2 ' + '-q ' + '-x ' + index_name + ' --seed ' + str(7) + ' --sensitive ' + '-k ' + str(max_map) + ' -5 ' + str(
        trimValue) + ' -U ' + readFile_name_2 + ' -p ' + str(args.np) + ' --no-unal ' + '> ' + alignment2_name

    command_bowtie_Seq = 'bowtie2-inspect -a 100 ' +  index_name +  ' > ' + myUnitig

    with open(bowtie_terminal_index, "w") as outfile:
        p_bowtie_0 = subprocess.run(command_bowtie_0, stdout=outfile, stderr=outfile, shell=True)

    with open(bowtie_terminal_1, "w") as outfile:
        p_bowtie_1 = subprocess.run(command_bowtie_1, stdout=outfile, stderr=outfile, shell=True)

    with open(bowtie_terminal_2, "w") as outfile:
        p_bowtie_2 = subprocess.run(command_bowtie_2, stdout=outfile, stderr=outfile, shell=True)
        
    with open(bowtie_terminal_Seq, "w") as outfile:
        p_bowtie_seq = subprocess.run(command_bowtie_Seq, stdout=outfile, stderr=outfile, shell=True)


if args.level <= 4: # NUCMER (Repeat detection):
    command_nucmer_1 = ['nucmer' + ' --maxmatch' + ' --nosimplify' + ' -l ' + str(kmer-1) + ' -c ' + str(kmer-1) + ' --prefix=' + nucmer_name + ' ' + genome_name + ' ' + myUnitig]
    command_nucmer_2 = ['show-coords' + ' -r ' + nucmer_name + '.delta' + ' > ' + nucmer_repeat]

    with open(nucmer_terminal_file1, "w") as outfile:
        p_nucmer1 = subprocess.run(command_nucmer_1, shell=True, stdout=outfile, stderr=outfile)

    with open(nucmer_terminal_file2, "w") as outfile:
        p_nucmer2 = subprocess.run(command_nucmer_2, shell=True, stdout=outfile, stderr=outfile)


# Remove the files that are not needed:

files_to_keep = [genome_name, wgsim_terminal_file, readFile_name_1, readFile_name_2, abyss_terminal_file, unitig_name, myUnitig, 
              bowtie_terminal_1, bowtie_terminal_2, alignment1_name, alignment2_name, nucmer_repeat, args_file_name_assembly]

current_dir = os.getcwd()
all_files = os.listdir(current_dir)
for file in all_files:
    if file not in files_to_keep:
        if not os.path.isdir(file):
            os.remove(os.path.join(current_dir, file))

