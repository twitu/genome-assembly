import random
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.colors import from_levels_and_colors

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
    plt.show()

if __name__ == "__main__":
    genome, reads = generate_reads()
    plot_reads(reads)

