#!/bin/bash

source ~/load_modules.sh

cd MULTIOME_FETAL

for i in {1..178}
do
Rscript MULTIOME_FETAL/LinkPeaksToGenes.R 'MULTIOME_FETAL/seurObj.RDS' \
	'pr_gns.RDS' \
	$i
done

echo 'DONE'

