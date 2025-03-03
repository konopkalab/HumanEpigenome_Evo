#!/bin/bash

module load python/3.6.4-anaconda
source activate ldsc

cd ~/workdir/pr5/05_LDSC/TOP20K_EXPAND_50kb

# Loop through the directories per peak set
# Each peak set belongs to a node
for j in $(ls | grep 'ADULT')
do

name=$j
cd ~/workdir/pr5/05_LDSC/TOP20K_EXPAND_50kb/$name

# Liftover from hg38 to hg19
~/workdir/programs/liftOver peakset.bed ~/workdir/reference_genomes/chain_files/hg38ToHg19.over.chain.gz $name"_hg19.bed" "unlifted_"$name

name=$name"_hg19"

	# Create annot file
	for i in {1..22}
	do
	python ~/workdir/programs/ldsc/./make_annot.py --bed-file $name.bed --bimfile ../../1000G_EUR_Phase3_plink/1000G.EUR.QC.$i.bim --annot-file $name.$i.annot.gz &
	done

	wait
	echo 'DONE_1'

	# Perform LDSC
	for i in {1..22}
	do
	python ~/workdir/programs/ldsc/./ldsc.py \
	--l2 \
	--bfile ../../1000G_EUR_Phase3_plink/1000G.EUR.QC.$i \
	--ld-wind-cm 1 \
	--annot $name.$i.annot.gz \
	--thin-annot \
	--out $name.$i \
	--print-snps ../../hapmap3_snps/hm.$i.snp
	done

done

wait

echo "DONE_ALL"



# END OF SCRIPT
