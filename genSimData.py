import os
import random
import numpy as np
import argparse as ap
import json

parser = ap.ArgumentParser()
parser.add_argument('--name', help='Dir', default="simulated")
parser.add_argument('--L', type=int, default=400,
                    help="Length of repeats")
parser.add_argument('--C', type=int, default=25,
                    help="Copy number of repeats ")
args = parser.parse_args()

args.name = args.name + '_L' + str(args.L) + '_C' + str(args.C)


def make_genome(fname, org):
    genome = "".join(random.choices(["A", "T", "G", "C"], [0.25, 0.25, 0.25, 0.25], k=(5000000-(args.L*args.C*2))))

    with open(fname.split(".")[0] + '_' + str(org) + ".fasta", "w+") as wf:
        wf.write(fname.split(".")[0] + "\n")
        wf.write(genome + "\n")
    return genome


def make_repeats(fnameR, N, probs, lens, org):
    Repeats = []
    with open(fnameR.split(".")[0] + '_' + str(org) + ".fasta", "w+") as wf:
        for i in range(N):
            prob_i = probs[i]
            len_i = lens[i]
            r_i = "".join(random.choices(["A", "T", "G", "C"], prob_i, k=len_i))
            wf.write(">" + str(i) + " " + fnameR.split(".")[0] + "_" + str(org) + "_" + str(i) + "\n")
            wf.write(r_i + "\n")
            Repeats.append(r_i)

    return Repeats


def generate(len_, count, mingap):
    n = count
    limit = len_
    mingap = mingap
    slack = limit - mingap * (n - 1)
    steps = random.randint(0, slack)
    increments = np.hstack([np.ones((steps,)), np.zeros((n,))])
    np.random.shuffle(increments)
    locs = np.argwhere(increments == 0).flatten()

    return np.cumsum(increments)[locs] + mingap * np.arange(0, n)


def insert_repeats(dir_, fnameR, genome, N_intra, Repeats_intra, copys_intra, N_inter, Repeats_inter, copy_inter, org):
    L = len(genome)
    N_copy = sum(copys_intra) + sum(copy_inter)
    mingap = int(0.5 * L / (N_copy - 1))
    pos_all = generate(len(genome), N_copy, mingap)

    with open(os.path.join(dir_, dir_ + '_inserted_offset_Repeats' + "_" + str(org) + '.fasta'), "w+") as wf:
        for i in range(N_intra):
            c_i = copys_intra[i]
            r_i = Repeats_intra[i]
            random_indices_i = random.sample(range(len(pos_all)), c_i)
            pos_i = [pos_all[i] for i in random_indices_i]
            pos_all = [x for x in pos_all if x not in pos_i]
            wf.write(">" + str(i) + " " + fnameR.split(".")[0] + "_intra_" + str(i) + "\n")
            for k in pos_i:
                k = int(k)
                genome = genome[:k] + r_i + genome[k:]
                wf.write(str(k) + "\n")

        for j in range(N_inter):
            c_j = copy_inter[j]
            r_j = Repeats_inter[j]
            random_indices_j = random.sample(range(len(pos_all)), c_j)
            pos_j = [pos_all[i] for i in random_indices_j]
            pos_all = [x for x in pos_all if x not in pos_j]
            wf.write(">" + str(j + N_intra) + " " + fnameR.split(".")[0] + "_inter_" + str(j) + "\n")
            for k in pos_j:
                k = int(k)
                genome = genome[:k] + r_j + genome[k:]
                wf.write(str(k) + "\n")

    return genome


os.chdir('Data')

try:
    os.makedirs(args.name, exist_ok=True)
    print("Directory '%s' created successfully" % args.name)
except OSError as error:
    print("Directory '%s' can not be created" % args.name)

cl_args = os.path.join(args.name, "commandline_args.txt")

with open(cl_args, 'w+') as f:
    json.dump(args.__dict__, f, indent=2)

fname = os.path.join(args.name, args.name + ".fasta")
fnameR_intra = os.path.join(args.name, args.name + "_intraRepeats.fasta")
fnameR_inter = os.path.join(args.name, args.name + "_interRepeats.fasta")


copys_intra = [[args.C],[args.C]]
copys_inter = [[args.C],[args.C]]

length_intra = [[args.L],[args.L]]
length_inter = [args.L]

probs_intra = [[0.2, 0.3, 0.3, 0.2], [0.4, 0.2, 0.3, 0.1]]
probs_inter = [[0.1, 0.3, 0.2, 0.4]]

Repeats_inter = make_repeats(fnameR_inter, 1, probs_inter, length_inter, 99)

with open(os.path.join(args.name, "ref_genome.fasta"), "w+") as wf:
    for i in range(2):
        genome_i = make_genome(fname, i)
        Repeats_intra_i = make_repeats(fnameR_intra, 1, probs_intra, length_intra[i], i)
        inserted_genome_i = insert_repeats(args.name, fnameR_intra, genome_i, 1, Repeats_intra_i,
                                           copys_intra[i], 1, Repeats_inter, copys_inter[i], i)
        
        with open(os.path.join(args.name, "ref_genome_" + str(i) + ".fasta"), "w+") as wf_2:
            wf_2.write(">ref_genome_" + str(i) + "\n")
            wf_2.write(inserted_genome_i + "\n")

        wf.write(">ref_genome_" + str(i) + "\n")
        wf.write(inserted_genome_i + "\n")


