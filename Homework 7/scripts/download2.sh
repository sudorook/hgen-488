#! /bin/bash
set -eu

urlprefix="http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/HNSC/20160128/"
archive="gdac.broadinstitute.org_HNSC.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2016012800.0.0.tar.gz"
dirname="gdac.broadinstitute.org_HNSC.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2016012800.0.0"
filename="HNSC.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt"

mkdir -p raw_data
cd raw_data

wget -nc "${urlprefix}${archive}"
tar xf "${archive}"

mkdir -p data
mv "${dirname}/${filename}" data/data.txt
