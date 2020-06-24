# This script takes a folder/single-file of CCMPred produced Coevolution Matrices, and a number N,
# to create N sized coevolution networks for the given proteins


import argparse

parser = argparse.ArgumentParser(description="Coevolution Network Creator Script")
# parser.add_argument("-p", "--path", help="Path to file/folder of the Coevolution Matrices", type=str, required=True)
parser.add_argument("-n", "--num", help="Number of nodes for the Coevolution Networks", type=int, required=True)
# parser.add_argument("-wd", "--workDir", help="Working directory", type=str, required=True)
parser.add_argument("-a", "--aln", help="Alignment file", type=str, default=False)
parser.add_argument("-cpu", '--threads', help="Number of threads to use. Default is 1", type=int, default=1)
parser.add_argument("-f", "--file", help="Path to file/folder of the Coevolution Matrices", type=str, required=True)


args = parser.parse_args()

# _PATH = args.path
_FILE = args.file
_NUM = args.num
# _WORKDIR = args.workDir
_ALN = args.aln
_THREADS = args.threads
_FILELIST = False
_ALNSEQS = False
_GRAPHPATH = False

########################################

########################################

import os
import numpy as np
import networkx as nx
from Bio import AlignIO
from multiprocessing import Queue, Process
from tqdm import tqdm
import time
from copy import deepcopy


def createNetwork(file):

    print(file)
    fileroot = file.split('.cv')[0]
    temp_mat = np.loadtxt(file)

    # print "Picking {} coevolving residue pairs from {}".format(_NUM, fileroot)

    indices = [[], []]

    for i in range(0, _NUM): # Picking the top N pairs

            index = np.where(temp_mat==np.max(temp_mat))
            indices[0] += [index[0][0]]
            indices[1] += [index[1][0]]

            temp_mat[index[0][0], index[1][0]] = 0
            temp_mat[index[1][0], index[0][0]] = 0

    zip_obj = zip(indices[0],indices[1]) # Reordering them
    indices[0] = [sorted(point)[0] for point in deepcopy(zip_obj)]
    indices[1] = [sorted(point)[1] for point in zip_obj]

    G = nx.Graph()

    for i in range(len(indices[0])):
        G.add_edge(fileroot + '-' + str(indices[0][i] + 1), fileroot + '-' + str(indices[1][i] + 1))

    nx.write_graphml(G, f'/Work/RRCoevNets/{fileroot}.graphml')

########################################

########################################

def main():

    createNetwork(_FILE)


if __name__ == '__main__':
    main()
