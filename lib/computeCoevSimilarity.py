

import argparse
import os
from multiprocessing import Queue, Process, Manager


parser = argparse.ArgumentParser(description="Coev Similarity Compute Script")
parser.add_argument("-wd", "--workDir", help="Working directory", type=str, required=True)
parser.add_argument("-a", "--aln", help="Alignment Graph", type=str, required=True)
parser.add_argument("-cpu", '--threads', help="Number of threads to use. Default is 1", type=int, default=1)
args = parser.parse_args()

_THREADS = args.threads
_NETDIR = args.workDir
_ALNPATH = args.aln
_FILELIST = sorted(os.listdir(_NETDIR))
_FILELIST = [os.path.join(_NETDIR, i) for i in _FILELIST]

manager = Manager()

max_clique_dict = manager.dict()
full_score_dict = manager.dict()
_POSITION = 0

########################################

########################################

import networkx as nx
from tqdm import tqdm
from copy import deepcopy
import numpy as np
import pandas as pd
import time

def checkAnalysisDir():

    	if not os.path.exists(_ANALYSISDIR):
    		os.makedirs(_ANALYSISDIR)


def prepareScoreDict(curr_file):

    score_dict = {}

    for file in _FILELIST:

        file_root = file.split('/')[-1].split('.')[0]

        if file_root == curr_file:
            score_dict[file_root] = 0
            continue

        score_dict[file_root] = np.nan

    return score_dict



def computeOneFile(file, aln_net):

    global max_clique_dict
    global full_score_dict

    file_root = file.split('/')[-1].split('.')[0]

    score_dict = prepareScoreDict(file_root)

    coev_net = nx.read_graphml(file)
    cliques = list(nx.find_cliques(coev_net))

    max_clique_dict[file_root] = nx.graph_number_of_cliques(coev_net)
    # print file_root
    # print max_clique_dict[file_root]


    for clique in cliques:

        flag = False

        res_count_dict = {}

        for res in clique:

            if res in list(nx.nodes(aln_net)):

                for neighbor in list(nx.all_neighbors(aln_net, res)):

                    prot = neighbor.split('-')[0]

                    if prot not in res_count_dict.keys():
                        res_count_dict[prot] = 1
                    else:
                        res_count_dict[prot] += 1
            else:
                flag = True
                break

        if flag:
            max_clique_dict[file_root] -= 1
            continue

        for key, value in res_count_dict.items():

            if value == len(clique):

                    if np.isnan(score_dict[key]):
                        score_dict[key] = 1
                    else:
                        score_dict[key] += 1

        full_score_dict[file_root] = score_dict


def workflow(files, aln_net, position):

    for file in files:
        file_root = file.split('/')[-1].split('.')[0]
        computeOneFile(file, aln_net)

    time.sleep(1)


def parallel():

    global max_clique_dict
    global full_score_dict

    print("Reading ALN Net")
    aln_net = nx.read_graphml(_ALNPATH)

    file_chunks = np.array_split(np.array(_FILELIST), _THREADS)
    positions = [1,2,3,4]
    processes = [Process(target=workflow, args=(list(chunk), aln_net, positions)) for chunk in file_chunks]

    print("Comparing Coevolution Networks")
    for p in processes:
        p.start()

    for p in processes:
        p.join()


    time.sleep(5)


def createMatrix():


    global max_clique_dict
    global full_score_dict

    dist_mat = np.zeros([len(full_score_dict.keys()), len(full_score_dict.keys())])
    jac_mat = np.zeros([len(full_score_dict.keys()), len(full_score_dict.keys())])

    counter = 0

    for key in sorted(full_score_dict.keys()):

        dist_mat[counter,] = [full_score_dict[key][inner_key] for inner_key in sorted(full_score_dict.keys())]
        counter += 1


    for i in range(len(dist_mat)):
        for j in range(i+1, len(dist_mat[0])):

            if dist_mat[i,j] == np.nan or dist_mat[j,i] == np.nan:
                dist_mat[j,i] = 0
                dist_mat[i,j] = 0
                jac_mat[i,j] = 0
                jac_mat[j,i] = 0

            if dist_mat[j,i] < dist_mat[i,j]:
                dist_mat[i,j] = dist_mat[j,i]

            else:
                dist_mat[j,i] = dist_mat[i,j]

            prot1 = sorted(full_score_dict.keys())[i]
            prot2 = sorted(full_score_dict.keys())[j]

            jac_mat[i,j] = float(dist_mat[i,j]) / (max_clique_dict[prot1] + max_clique_dict[prot2] - float(dist_mat[i,j]))
            jac_mat[j,i] = float(dist_mat[i,j]) / (max_clique_dict[prot1] + max_clique_dict[prot2] - float(dist_mat[i,j]))


    dist_df = pd.DataFrame(dist_mat, columns=sorted(full_score_dict.keys()), index = sorted(full_score_dict.keys()))
    jac_df = pd.DataFrame(jac_mat, columns=sorted(full_score_dict.keys()), index = sorted(full_score_dict.keys()))

    dist_df.to_csv('./mat.csv')
    jac_df.to_csv('./jac.csv')



########################################

########################################


def main():

    # checkAnalysisDir()
    parallel()
    createMatrix()


if __name__ == '__main__':
    main()
