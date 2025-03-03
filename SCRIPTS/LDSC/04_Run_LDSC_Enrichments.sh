#!/bin/bash

module load python/3.6.4-anaconda
source activate ldsc


type=TOP20K_EXPAND_50kb
cd /home2/s422159/workdir/pr5/05_LDSC/$type

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
	python ~/workdir/programs/ldsc/./ldsc.py \
	--h2-cts ~/workdir/GWAS/DATASETS/ADHD/ADHD_2017_Modified.sumstats \
	--ref-ld-chr ~/workdir/pr5/05_LDSC/1000G_EUR_Phase3_baseline/baseline. \
	--out ~/workdir/pr5/05_LDSC/$type/RESULTS/$name/$name"_ADHD" \
	--ref-ld-chr-cts ~/workdir/pr5/05_LDSC/$type/$name".ldcts" \
	--w-ld-chr ~/workdir/pr5/05_LDSC/weights_hm3_no_hla/weights. &

	# AUTISM (hg19)
	python ~/workdir/programs/ldsc/./ldsc.py \
	--h2-cts ~/workdir/GWAS/DATASETS/AUTISM/iPSYCH-PGC_ASD_Nov2017_Modified.sumstats \
	--ref-ld-chr ~/workdir/pr5/05_LDSC/1000G_EUR_Phase3_baseline/baseline. \
	--out ~/workdir/pr5/05_LDSC/$type/RESULTS/$name/$name"_AUT" \
	--ref-ld-chr-cts ~/workdir/pr5/05_LDSC/$type/$name".ldcts" \
	--w-ld-chr ~/workdir/pr5/05_LDSC/weights_hm3_no_hla/weights. &

	# BP 2021 (hg19)
	python ~/workdir/programs/ldsc/./ldsc.py \
	--h2-cts ~/workdir/GWAS/DATASETS/BP/PGC_2021/pgc-bip2021-all_Modified.sumstats.gz \
	--ref-ld-chr ~/workdir/pr5/05_LDSC/1000G_EUR_Phase3_baseline/baseline. \
	--out ~/workdir/pr5/05_LDSC/$type/RESULTS/$name/$name"_BP" \
	--ref-ld-chr-cts ~/workdir/pr5/05_LDSC/$type/$name".ldcts" \
	--w-ld-chr ~/workdir/pr5/05_LDSC/weights_hm3_no_hla/weights. &

	# SCZ_2022 (hg19)
	python ~/workdir/programs/ldsc/./ldsc.py \
	--h2-cts ~/workdir/GWAS/DATASETS/SCZ/SCZ_TRUBETSKOY_2022/PGC3_SCZ_wave3_Modified.sumstats.gz \
	--ref-ld-chr ~/workdir/pr5/05_LDSC/1000G_EUR_Phase3_baseline/baseline. \
	--out ~/workdir/pr5/05_LDSC/$type/RESULTS/$name/$name"_SCZ" \
	--ref-ld-chr-cts ~/workdir/pr5/05_LDSC/$type/$name".ldcts" \
	--w-ld-chr ~/workdir/pr5/05_LDSC/weights_hm3_no_hla/weights. &

	# MDD
	python ~/workdir/programs/ldsc/./ldsc.py \
	--h2-cts ~/workdir/GWAS/DATASETS/MDD/HOWARD_2019/PGC_UKB_depression_genome-wide_Modified.sumstats.gz \
	--ref-ld-chr ~/workdir/pr5/05_LDSC/1000G_EUR_Phase3_baseline/baseline. \
	--out ~/workdir/pr5/05_LDSC/$type/RESULTS/$name/$name"_MDD" \
	--ref-ld-chr-cts ~/workdir/pr5/05_LDSC/$type/$name".ldcts" \
	--w-ld-chr ~/workdir/pr5/05_LDSC/weights_hm3_no_hla/weights. &

	# AD
	python ~/workdir/programs/ldsc/./ldsc.py \
	--h2-cts ~/workdir/GWAS/DATASETS/AD/JANSEN_2019/AD_Jansen2019_Modified.sumstats.gz \
	--ref-ld-chr ~/workdir/pr5/05_LDSC/1000G_EUR_Phase3_baseline/baseline. \
	--out ~/workdir/pr5/05_LDSC/$type/RESULTS/$name/$name"_AD" \
	--ref-ld-chr-cts ~/workdir/pr5/05_LDSC/$type/$name".ldcts" \
	--w-ld-chr ~/workdir/pr5/05_LDSC/weights_hm3_no_hla/weights. &

	# OCD
	python ~/workdir/programs/ldsc/./ldsc.py \
	--h2-cts ~/workdir/GWAS/DATASETS/PGC_OCD_Aug2017/ocd_aug2017_Modified.sumstats.gz \
	--ref-ld-chr ~/workdir/pr5/05_LDSC/1000G_EUR_Phase3_baseline/baseline. \
	--out ~/workdir/pr5/05_LDSC/$type/RESULTS/$name/$name"_OCD" \
	--ref-ld-chr-cts ~/workdir/pr5/05_LDSC/$type/$name".ldcts" \
	--w-ld-chr ~/workdir/pr5/05_LDSC/weights_hm3_no_hla/weights. &

	# Anxiety
	python ~/workdir/programs/ldsc/./ldsc.py \
	--h2-cts ~/workdir/GWAS/DATASETS/PANIC_DISORDER/PGC_2019/pgc-panic2019_Modified.sumstats.gz \
	--ref-ld-chr ~/workdir/pr5/05_LDSC/1000G_EUR_Phase3_baseline/baseline. \
	--out ~/workdir/pr5/05_LDSC/$type/RESULTS/$name/$name"_ANXIETY" \
	--ref-ld-chr-cts ~/workdir/pr5/05_LDSC/$type/$name".ldcts" \
	--w-ld-chr ~/workdir/pr5/05_LDSC/weights_hm3_no_hla/weights. &

	# COGNITIVE_FUNCTION
	python ~/workdir/programs/ldsc/./ldsc.py \
	--h2-cts ~/workdir/GWAS/DATASETS/COGNITIVE_FUNCTION/Cognitive_Modified.sumstats \
	--ref-ld-chr ~/workdir/pr5/05_LDSC/1000G_EUR_Phase3_baseline/baseline. \
	--out ~/workdir/pr5/05_LDSC/$type/RESULTS/$name/$name"_COG" \
	--ref-ld-chr-cts ~/workdir/pr5/05_LDSC/$type/$name".ldcts" \
	--w-ld-chr ~/workdir/pr5/05_LDSC/weights_hm3_no_hla/weights. &

	# Intelligence
	python ~/workdir/programs/ldsc/./ldsc.py \
	--h2-cts ~/workdir/GWAS/DATASETS/Intelligence/Intelligence_Modified.sumstats \
	--ref-ld-chr ~/workdir/pr5/05_LDSC/1000G_EUR_Phase3_baseline/baseline. \
	--out ~/workdir/pr5/05_LDSC/$type/RESULTS/$name/$name"_INT" \
	--ref-ld-chr-cts ~/workdir/pr5/05_LDSC/$type/$name".ldcts" \
	--w-ld-chr ~/workdir/pr5/05_LDSC/weights_hm3_no_hla/weights. &

	# CAD
	python ~/workdir/programs/ldsc/./ldsc.py \
	--h2-cts ~/workdir/GWAS/DATASETS/NON_BRAIN/CAD/CARDIoGRAM_GWAS_Modified.sumstats \
	--ref-ld-chr ~/workdir/pr5/05_LDSC/1000G_EUR_Phase3_baseline/baseline. \
	--out ~/workdir/pr5/05_LDSC/$type/RESULTS/$name/$name"_CAD" \
	--ref-ld-chr-cts ~/workdir/pr5/05_LDSC/$type/$name".ldcts" \
	--w-ld-chr ~/workdir/pr5/05_LDSC/weights_hm3_no_hla/weights. &

	# OST
	python ~/workdir/programs/ldsc/./ldsc.py \
	--h2-cts ~/workdir/GWAS/DATASETS/NON_BRAIN/OSTEOPOROSIS/GEFOS2_FNBMD_POOLED_Modified.sumstats \
	--ref-ld-chr ~/workdir/pr5/05_LDSC/1000G_EUR_Phase3_baseline/baseline. \
	--out ~/workdir/pr5/05_LDSC/$type/RESULTS/$name/$name"_OST" \
	--ref-ld-chr-cts ~/workdir/pr5/05_LDSC/$type/$name".ldcts" \
	--w-ld-chr ~/workdir/pr5/05_LDSC/weights_hm3_no_hla/weights. &

	wait
	printf "\n\n\nLOOPING\n\n\n"

done



