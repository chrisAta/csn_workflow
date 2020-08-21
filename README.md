# csn_workflow

## Requirements

This workflow uses a self-contained Nextflow pipeline that works in a Docker container. Nextflow and Docker are the only requirements for using this workflow.

## Usage

To create Coevolution Similarity Networks, you first need to produce coevolution matrices by CCMPred for a set of protein sequences ([instructions here](https://github.com/soedinglab/CCMpred/wiki/FAQ)). You can download an already produced zipped folder of such matrices [here](https://data.ncl.ac.uk/articles/dataset/Trans241CoevMatrices_tar_gz/12555620). Decompress it.

You also need a Multiple Sequence Alignment (MSA) for your sequences in FASTA format. I use [ClustalOmega](https://www.ebi.ac.uk/Tools/msa/clustalo/) for this. The MSA for the test dataset mentioned above is in this repo in 'example/trans241_aln.fasta'.

## Arguments
```
--path: Path to the directory containing all the Coevolution Matrices produced by CCMPred
--aln: Path to the MSA of the sequences
--n: Top N hits to take from each Coevolution Matrix to form a Residue-Residue Coevolution Network
--f: Filter out Coevolving Pairs that appear in at least F of the Proteins in the folder
--t: Coevolution similarity threshold to use to produce the final Coevolution Similarity Networks. Use 'all' to produce CSNs for all thresholds from 0.05 to 0.09 with increments of 0.01
--outdir: Path to output directory
--cpu: Number of threads to use
--stage: Number letting you stop the workflow at early stages. Stage 1 stops at just the Residue-Residue Coevolution Networks. Stage 2 or 3 go through the whole workflow. Stage 4 only does the last step of the workflow, which produces CSNs from a Coevolution Similarity Matrix. Can only be used with the --csm parameter
--csm: Path to Coevolution Similarity Matrix to be used for --stage=4
```

## Example command

`
nextflow run chrisAta/csn_workflow --path ./Trans241CoevMatrices --n 700 --cpu 4 --aln ./example/trans241_aln.fasta --f 150 --outdir ../test --stage 3 --t all
`

## Description of results

In the directory provided to the --outdir parameter, there are four resulting folders:

* aln_dir: This folder contains two files. 'aln_net.graphml' is the alignment metanetwork the workflow produces. 'coev_freqs.txt' contains the coevolving pairs  based on the columns of the MSA given, how many sequences contain each coevolving pair, and then a list of the sequences that contain each coevolving pair with the local positions of the residues for each sequence.
* coev_sim: This folder contains one file. 'jac.csv' is the Coevolution Similarity Matrix produced by the workflow. This file can be reused using the --csm and --stage=4 parameters.
* csn: This folder contains all the CSNs for the requested thresholds in Graphml format.
* graphml: This folder contains the Residue-Residue Coevolution Network for every sequence. The workflow can be stopped early to only produce these by using --stage=1
