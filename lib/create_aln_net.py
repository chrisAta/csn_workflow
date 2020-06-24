# This script produces an Alignment network based on the Coevolution Networks produced by coev_net_creator.py

import argparse
from Bio import AlignIO
import os

parser = argparse.ArgumentParser(description="Alignment Network Creator Script")
parser.add_argument("-wd", "--workDir", help="RRCoevNets folder", type=str, required=True)
parser.add_argument("-f", "--filter", help="Filter number to ignore coevolving columns that occur too frequently", type=int, required=True)
parser.add_argument("-a", "--aln", help="Alignment file", type=str, default=False)
args = parser.parse_args()

_WORKDIR = args.workDir
_FILTER = args.filter
_ALNFILE = AlignIO.read(args.aln, "fasta")
_ALNFILE.sort()
_FILELIST = sorted(os.listdir(_WORKDIR))
_FILELIST = [os.path.join(_WORKDIR, i) for i in _FILELIST]

pair_dict = {}

########################################

########################################

import networkx as nx
from itertools import combinations
from tqdm import tqdm
import time
import os.path

def createAlnMap(index):

    aln_map = {}
    gapped_seq = str(_ALNFILE[index].seq)
    ungapped_seq = str(_ALNFILE[index].seq.ungap('-'))

    count = 1

    for i in range(len(gapped_seq)):

        if gapped_seq[i] == '-':
            continue

        aln_map[str(count)] = str(i + 1)
        count += 1

    return aln_map


def readCoevNetwork(file, aln_map):

    G = nx.read_graphml(file)

    for edge in G.edges():

        source = int(edge[0].split('-')[-1])
        sink = int(edge[1].split('-')[-1])

        if source > sink:
            temp = sink
            sink = source
            source = temp

        pair = aln_map[str(source)] + '-' + aln_map[str(sink)]

        if pair not in pair_dict.keys():

            if int(edge[0].split('-')[1]) < int(edge[1].split('-')[1]):
                pair_dict[pair] = [(edge[0], edge[1])]
            else:
                pair_dict[pair] = [(edge[1], edge[0])]
        else:

            if int(edge[0].split('-')[1]) < int(edge[1].split('-')[1]):
                pair_dict[pair] += [(edge[0], edge[1])]
            else:
                pair_dict[pair] += [(edge[1], edge[0])]


def read_partition_file():

    fr = open("./partitions.txt", "r")

    part_count_dict = {}
    prot_part_dict = {}

    for line in fr:

        line = line.strip()
        temp_arr = line.split(":")

        if temp_arr[0] == "PART":
            part_count_dict[temp_arr[1]] = temp_arr[2]
        elif temp_arr[0] == "PROT":
            prot_part_dict[temp_arr[1]] = temp_arr[2]

    return part_count_dict, prot_part_dict


def createALNGraph():

    global _FILTER

    G = nx.Graph()
    fw = open("./coev_freqs.txt", "w")

    prot_part_dict = {}
    part_count_dict = {}

    if os.path.isfile("./partitions.txt"):
        part_count_dict, prot_part_dict = read_partition_file()

    full_dict = {}

    counter = 0

    for key in tqdm(pair_dict.keys()):

        value = pair_dict[key]

        fw.write(key + ' ; ' + str(len(value)) + ' ; ' + str(value) + '\n')

        if prot_part_dict != {}:
            thr = 0.6
            temp_arr = []
            max_part = -1
            max_val = -1

            for prot in value:
                temp_arr += [prot_part_dict[prot[0].split('-')[0]]]

            for part in part_count_dict.keys():
                if temp_arr.count(part) > max_val:
                    max_part = part
                    max_val = temp_arr.count(part)

            _FILTER = thr * int(part_count_dict[max_part])

        if len(value) <= _FILTER:

            indicesA = [temp[0] for temp in value]
            indicesB = [temp[1] for temp in value]

            edgesA = list(combinations(indicesA, 2))
            edgesB = list(combinations(indicesB, 2))

            if not edgesA == [] and not edgesB == []:

                for edge in edgesA:
                    G.add_edge(edge[0], edge[1])

                for edge in edgesB:
                    G.add_edge(edge[0], edge[1])

        counter += 1

    fw.close()
    return G


def workflow():

    print("Creating Pair Dict")

    tqdm.write('')

    for i in tqdm(range(len(_FILELIST))):

        aln_map = createAlnMap(i)
        readCoevNetwork(_FILELIST[i], aln_map)

    tqdm.write('')
    time.sleep(1)

    print("Computing ALN Graph")
    G = createALNGraph()

    out_path = './aln_net.graphml'

    print(f"Writing ALN Graph to {out_path}")
    nx.write_graphml(G, out_path)


########################################

########################################


def main():
    workflow()

if __name__ == '__main__':
    main()
