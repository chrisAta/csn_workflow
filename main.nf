#!/usr/bin/env nextflow

params.path = null
params.aln = null
params.n = null
params.f = null
params.cpu = 1
params.outdir = null
params.t = null
params.csm = false
params.stage = null

enzymes = Channel.fromPath( params.path + "*")
aln = file( params.aln )
n = params.n
cpu = params.cpu
f = params.f
t = params.t
stage = params.stage

if (params.stage == 4){
  csm2 = Channel.fromPath ( params.csm )
} else {
  csm2 = Channel.empty()
}

outdir = file( params.outdir )
outdir.mkdirs()
graphml = file(outdir.name + '/graphml')
graphml.mkdirs()
aln_dir = file(outdir.name + '/aln_dir')
aln_dir.mkdirs()
coev_sim = file(outdir.name + '/coev_sim')
coev_sim.mkdirs()
csn = file(outdir.name + '/csn')
csn.mkdirs()

process createCoevNets{
    // echo true
    publishDir graphml, mode : "copy"

    input:
    file matrix from enzymes

    output:
    file "*.graphml" into coev_net

    when:
    stage != 4

    """
    python3 /usr/bin/coev_net_creator.py -f ${matrix} -n ${n} -a ${aln} -cpu ${cpu}
    cp /Work/RRCoevNets/* .
    """
}

coev_net
  .collect()
  .set { coev_net_set }

process createAlignmentNet{
  echo true
  publishDir aln_dir, mode : "copy"


  input:
  file coev from coev_net_set
  file aln
  output:
  file("aln_net.graphml") into aln_net
  file("coev_freqs.txt")
  when:
  stage !=4 && stage!=1

  """
  mkdir coev_nets
  mv *.graphml coev_nets
  python3 /usr/bin/create_aln_net.py -wd coev_nets -f ${f} -a ${aln}
  """

}

process createCoevSimMatrix {
  echo true
  publishDir coev_sim, mode : "copy"

  input:
  file aln_net
  file coev from coev_net_set
  output:
  file "jac.csv" into csm1
  when:
  stage !=4 && stage!=1

  """
  mkdir aln_dir
  mkdir coev_nets
  mv aln_net.graphml aln_dir
  mv *.graphml coev_nets
  python3 /usr/bin/computeCoevSimilarity.py -wd coev_nets -a ./aln_dir/aln_net.graphml -cpu ${cpu}
  """

}

if (t == "all"){
  t = Channel.from(5..90)
      .map{ it / 100.0 }
}

csm = csm2.mix(csm1).first()

process createCSN{
  publishDir csn, mode: "copy"

  input:
  file csm
  val threshold from t
  when:
  stage==2 || stage==3 || stage==4

  output:
  file "test_${threshold}.graphml" into csn_out

  """
  python3 /usr/bin/csnCreator.py -m ${csm} -t ${threshold} -o test_${threshold}.graphml
  """

}
