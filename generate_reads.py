import random
import numpy as np
import matplotlib as mpl
import os
import re

from matplotlib import pyplot as plt
from matplotlib.colors import from_levels_and_colors

def rev_comp(genome_part):
    comp = {"T":"A", "G":"C", "C":"G", "A":"T"}
    return ''.join(list(map(lambda x: comp[x], genome_part)))

def plot_unitigs(genome, reads, read_len=30, filename="unitigs.txt"):
    mmers = {}
    mmer = ""
    line_no = 0
    char_no = 0
    unitigs = []

    with open(filename, "r") as fin:
        lines = fin.read().splitlines()

        while (line_no < len(lines)):
            # if line is clear mmer and skip line
            if not lines[line_no]:
                mmer = ""
                line_no += 1
                continue
            # if mmer is blank read new mmer
            if not mmer:
                mmer = lines[line_no]
                if mmer not in mmers:
                    mmers[mmer] = 0

                line_no += 1

            # if read mmer read kmer and read ids
            else:
                kmer = lines[line_no]
                mmers[mmer] += 1
                line_no += 1
                char_no = 0

                unitig = {}
                while (char_no < len(kmer)):
                    read_ids = list(map(int, lines[line_no].rstrip().split(" ")))
                    for read_id in read_ids:
                        if read_id not in unitig:
                            unitig[read_id] = ""
                    
                    for key in unitig.keys():
                        if key in read_ids:
                            unitig[key] = unitig[key] + kmer[char_no]
                        else:
                            unitig[key] = unitig[key] + " "

                    line_no += 1
                    char_no += 1

                # split across reads ids to find portions of genome that contribute to 
                unitigs.append(unitig)
            
    # generate graph for mmers
    mpl.rcParams['xtick.labelsize'] = 8
    plt.bar(range(len(mmers)), list(mmers.values()), align='center')
    plt.xticks(range(len(mmers)), list(mmers.keys()), rotation='vertical', ) # reduce font size
    plt.savefig("mmers.png")

    # generate graphs for unitigs
    matrix = np.zeros((len(unitigs), len(genome)), dtype=int)
    for i, unitig in enumerate(unitigs):
        for key, val in unitig.items():
            for part in val.split(" "):
                start = reads[key]
                index = start + genome[start : start + read_len].find(part)
                if not index:
                    index = start + rev_comp(genome[start:start+read_len]).find(part)
                for j in range(index, index + len(part)):
                    matrix[i][j] = 1
    # define the colors
    cmap, norm = from_levels_and_colors([0, 0.5, 1], ['r', 'k'])

    # create a normalize object the describes the limits of
    # each color
    bounds = [0., 0.5, 1.]
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    plt.figure(figsize=(20, 10))
    plt.imshow(matrix, interpolation='nearest', cmap=cmap, norm=norm)
    plt.savefig("kmers.png")
    plt.close()

def generate_reads(genome_len = 500, read_len = 30, read_num = 20, seed = 20):

    letters = ['A', 'C', 'G', 'T']
    genome = ''.join([random.choice(letters) for _ in range(genome_len)])
    random.seed(seed)
    mode = random.randint(0, genome_len-1-read_len)
    reads = []

    # use triangular probability to generate reads
    # this ensures high probability of intersection between reads
    for _ in range(read_num):
        mode = int(random.triangular(0, genome_len-1-read_len, mode))
        reads.append(mode)

    return genome, reads

def write_reads(genome, reads, read_len=30, filename="reads.txt"):
    with open(filename, "w") as f:
        for read in reads:
            f.write(genome[read:read + read_len] + "\n")


def plot_reads(reads, genome_len = 500, read_len = 30):

    # create binary matrix to represent read
    matrix = np.zeros((len(reads), genome_len), dtype=int)
    for i, read in enumerate(reads):
        for j in range(read, read + read_len):
            matrix[i][j] = 1

    # define the colors
    cmap, norm = from_levels_and_colors([0, 0.5, 1], ['r', 'k'])

    # create a normalize object the describes the limits of
    # each color
    bounds = [0., 0.5, 1.]
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

    # plot it
    plt.figure(figsize=(20, 10))
    plt.imshow(matrix, interpolation='nearest', cmap=cmap, norm=norm)
    plt.savefig("reads.png")
    plt.close()

if __name__ == "__main__":
    genome, reads = generate_reads()
    write_reads(genome, reads)
    plot_reads(reads)
    os.system("./a.out > unitigs.txt")
    plot_unitigs(genome, reads)
