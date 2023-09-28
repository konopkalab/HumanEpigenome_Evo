#!/bin/bash

source ~/load_modules.sh

cd MULTIOME_FETAL

for i in {1..178}
do
Rscript MULTIOME_FETAL/LinkPeaksToGenes.R 'MULTIOME_FETAL/seurObj.RDS' \
	'/home2/s422159/workdir/pr5/MULTIOME_FETAL/pr_gns.RDS' \
	$i
done

echo 'DONE'

