# csn_workflow

## Requirements

This workflow uses a self-contained Nextflow pipeline that works in a Docker container. Nextflow and Docker are the only requirements for using this workflow.

## Usage

To create Coevolution Similarity Networks, you first need to produce coevolution matrices by CCMPred for a set of protein sequences ([instructions here](https://github.com/soedinglab/CCMpred/wiki/FAQ)). You can download an already produced zipped folder of such matrices [here](http://nextcloud.bioswarm.net/s/2nskyeNXCNiH6J5). Decompress it.

You also need a Multiple Sequence Alignment for your sequences in FASTA format. I use [ClustalOmega](https://www.ebi.ac.uk/Tools/msa/clustalo/) for this. The MSA for the test dataset mentioned above is in this repo in 'example/trans241_aln.fasta'.

Here is an example command:

`
nextflow run chrisAta/csn_workflow --path ./Trans241CoevMatrices --n 700 --cpu 4 --aln ./example/trans241_aln.fasta --f 150 --outdir ../test --stage 3 --t all
`
