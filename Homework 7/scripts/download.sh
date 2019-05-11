#! /bin/bash
set -eu

sigfile="../data/sig_genes.txt"
allfile="../data/all_genes.txt"

rm -f "${sigfile}"
rm -f "${allfile}"

urlprefix="http://gdac.broadinstitute.org/runs/analyses__2016_01_28/reports/cancer/SKCM-TM/Pathway_GSEA_mRNAseq/"
prefix="SKCM-TM-mRNAseq_cNMF-clus"
suffix="___Class2-CanonicalPathway"

mkdir -p ../raw_data/
cd ../raw_data/
for i in $(seq 1 4); do
  wget -nc "${urlprefix}${prefix}${i}${suffix}.significant.genes.txt"
  wget -nc "${urlprefix}${prefix}${i}${suffix}.up.regulated.genes.txt"
  wget -nc "${urlprefix}${prefix}${i}${suffix}.down.regulated.genes.txt"

  if [[ ${i} = 1 ]]; then
    cat "${prefix}${i}${suffix}.significant.genes.txt" >> ${sigfile}
    cat "${prefix}${i}${suffix}.up.regulated.genes.txt" >> ${allfile}
    cat "${prefix}${i}${suffix}.down.regulated.genes.txt" >> ${allfile}
  else
    cat "${prefix}${i}${suffix}.significant.genes.txt" | \
      tail -n +2 >> ${sigfile}
    cat "${prefix}${i}${suffix}.up.regulated.genes.txt" | \
      tail -n +2 >> ${allfile}
    cat "${prefix}${i}${suffix}.down.regulated.genes.txt" | \
      tail -n +2 >> ${allfile}
  fi
done
