#!/bin/bash

module load python/3.6.4-anaconda
source activate ldsc


type=TOP20K_EXPAND_50kb
cd $type

# Copy files required for ldsc
cp ../ADULT.ldcts ADULT.ldcts
mkdir RESULTS
cp -r ../BCG/ADULT_BCG .

# Run the regressions

for i in $(ls -d ADULT_* | grep 'ADULT' | grep -v -P 'BCG|.*ldcts.*')
do
	name=$i
	mkdir RESULTS/$name
	sed "s/CHANGE/$name/g" ADULT.ldcts > $name".ldcts"

	# ADHD (hg19)
	python ldsc/./ldsc.py \
	--h2-cts GWAS/DATASETS/ADHD/ADHD_2017_Modified.sumstats \
	--ref-ld-chr 1000G_EUR_Phase3_baseline/baseline. \
	--out $type/RESULTS/$name/$name"_ADHD" \
	--ref-ld-chr-cts $type/$name".ldcts" \
	--w-ld-chr weights_hm3_no_hla/weights. &

	# AUTISM (hg19)
	python ldsc/./ldsc.py \
	--h2-cts GWAS/DATASETS/AUTISM/iPSYCH-PGC_ASD_Nov2017_Modified.sumstats \
	--ref-ld-chr 1000G_EUR_Phase3_baseline/baseline. \
	--out $type/RESULTS/$name/$name"_AUT" \
	--ref-ld-chr-cts $type/$name".ldcts" \
	--w-ld-chr weights_hm3_no_hla/weights. &

	# BP 2021 (hg19)
	python ldsc/./ldsc.py \
	--h2-cts GWAS/DATASETS/BP/PGC_2021/pgc-bip2021-all_Modified.sumstats.gz \
	--ref-ld-chr 1000G_EUR_Phase3_baseline/baseline. \
	--out $type/RESULTS/$name/$name"_BP" \
	--ref-ld-chr-cts $type/$name".ldcts" \
	--w-ld-chr weights_hm3_no_hla/weights. &

	# SCZ_2022 (hg19)
	python ldsc/./ldsc.py \
	--h2-cts GWAS/DATASETS/SCZ/SCZ_TRUBETSKOY_2022/PGC3_SCZ_wave3_Modified.sumstats.gz \
	--ref-ld-chr 1000G_EUR_Phase3_baseline/baseline. \
	--out $type/RESULTS/$name/$name"_SCZ" \
	--ref-ld-chr-cts $type/$name".ldcts" \
	--w-ld-chr weights_hm3_no_hla/weights. &

	# MDD
	python ldsc/./ldsc.py \
	--h2-cts GWAS/DATASETS/MDD/HOWARD_2019/PGC_UKB_depression_genome-wide_Modified.sumstats.gz \
	--ref-ld-chr 1000G_EUR_Phase3_baseline/baseline. \
	--out $type/RESULTS/$name/$name"_MDD" \
	--ref-ld-chr-cts $type/$name".ldcts" \
	--w-ld-chr weights_hm3_no_hla/weights. &

	# AD
	python ldsc/./ldsc.py \
	--h2-cts GWAS/DATASETS/AD/JANSEN_2019/AD_Jansen2019_Modified.sumstats.gz \
	--ref-ld-chr 1000G_EUR_Phase3_baseline/baseline. \
	--out $type/RESULTS/$name/$name"_AD" \
	--ref-ld-chr-cts $type/$name".ldcts" \
	--w-ld-chr weights_hm3_no_hla/weights. &

	# OCD
	python ldsc/./ldsc.py \
	--h2-cts GWAS/DATASETS/PGC_OCD_Aug2017/ocd_aug2017_Modified.sumstats.gz \
	--ref-ld-chr 1000G_EUR_Phase3_baseline/baseline. \
	--out $type/RESULTS/$name/$name"_OCD" \
	--ref-ld-chr-cts $type/$name".ldcts" \
	--w-ld-chr weights_hm3_no_hla/weights. &

	# Anxiety
	python ldsc/./ldsc.py \
	--h2-cts GWAS/DATASETS/PANIC_DISORDER/PGC_2019/pgc-panic2019_Modified.sumstats.gz \
	--ref-ld-chr 1000G_EUR_Phase3_baseline/baseline. \
	--out $type/RESULTS/$name/$name"_ANXIETY" \
	--ref-ld-chr-cts $type/$name".ldcts" \
	--w-ld-chr weights_hm3_no_hla/weights. &

	# COGNITIVE_FUNCTION
	python ldsc/./ldsc.py \
	--h2-cts GWAS/DATASETS/COGNITIVE_FUNCTION/Cognitive_Modified.sumstats \
	--ref-ld-chr 1000G_EUR_Phase3_baseline/baseline. \
	--out $type/RESULTS/$name/$name"_COG" \
	--ref-ld-chr-cts $type/$name".ldcts" \
	--w-ld-chr weights_hm3_no_hla/weights. &

	# Intelligence
	python ldsc/./ldsc.py \
	--h2-cts GWAS/DATASETS/Intelligence/Intelligence_Modified.sumstats \
	--ref-ld-chr 1000G_EUR_Phase3_baseline/baseline. \
	--out $type/RESULTS/$name/$name"_INT" \
	--ref-ld-chr-cts $type/$name".ldcts" \
	--w-ld-chr weights_hm3_no_hla/weights. &

	# CAD
	python ldsc/./ldsc.py \
	--h2-cts GWAS/DATASETS/NON_BRAIN/CAD/CARDIoGRAM_GWAS_Modified.sumstats \
	--ref-ld-chr 1000G_EUR_Phase3_baseline/baseline. \
	--out $type/RESULTS/$name/$name"_CAD" \
	--ref-ld-chr-cts $type/$name".ldcts" \
	--w-ld-chr weights_hm3_no_hla/weights. &

	# OST
	python ldsc/./ldsc.py \
	--h2-cts GWAS/DATASETS/NON_BRAIN/OSTEOPOROSIS/GEFOS2_FNBMD_POOLED_Modified.sumstats \
	--ref-ld-chr 1000G_EUR_Phase3_baseline/baseline. \
	--out $type/RESULTS/$name/$name"_OST" \
	--ref-ld-chr-cts $type/$name".ldcts" \
	--w-ld-chr weights_hm3_no_hla/weights. &

	# HEIGHT
	python ldsc/./ldsc.py \
	--h2-cts GWAS/DATASETS/NON_BRAIN/HEIGHT/Height_UKBiobank_2018_Modified.sumstats \
	--ref-ld-chr 1000G_EUR_Phase3_baseline/baseline. \
	--out $type/RESULTS/$name/$name"_HEIGHT" \
	--ref-ld-chr-cts $type/$name".ldcts" \
	--w-ld-chr weights_hm3_no_hla/weights. &

	wait
	printf "\n\n\nLOOPING\n\n\n"

done




