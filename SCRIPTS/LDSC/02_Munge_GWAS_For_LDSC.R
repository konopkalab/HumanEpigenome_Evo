library(plyr)
library(dplyr)
library(tidyverse)
library(tidyr)
library(reshape2)
library(data.table)

# LDSC regression requires GWAS summary statistics to include certain columns.
# Below code modifies the summary stats for "munge_sumstats.py" script

# Before running this script, please do the following in bash.
# Download hapmap3 snps
#$wget https://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist.bz2
#$bunzip2 w_hm3.snplist.bz2

#$module load python/3.6.4-anaconda
#$source activate ldsc

####
## IN R
####

# ADHD
# PMID: 29325848. N=55374

adhd = fread('~/workdir/GWAS/DATASETS/ADHD/adhd_eur_jun2017', header = T)
adhd$N = 55374
fwrite(adhd, '~/workdir/GWAS/DATASETS/ADHD/adhd_eur_jun2017_Modified', col.names = T, row.names = F, sep = '\t')

python ~/workdir/programs/ldsc/./munge_sumstats.py \
--sumstats adhd_eur_jun2017_Modified \
--merge-alleles ../w_hm3.snplist \
--out ADHD_2017_Modified


# AUTISM
# PMID: 30804558. N=46350

aut = fread('~/workdir/GWAS/DATASETS/AUTISM/iPSYCH-PGC_ASD_Nov2017.bed', header = T)
aut$N = 46350
fwrite(aut, 'workdir/GWAS/DATASETS/AUTISM/iPSYCH-PGC_ASD_Nov2017_Modified.bed', col.names = T, row.names = F, sep = '\t')

python ~/workdir/programs/ldsc/./munge_sumstats.py \
--sumstats iPSYCH-PGC_ASD_Nov2017_Modified.bed \
--merge-alleles ../w_hm3.snplist \
--out iPSYCH-PGC_ASD_Nov2017_Modified


# BIPOLAR DISORDER 2021
# PMID:  34002096 (2021), N=41917 + 371549 = 413466

bd = fread('~/workdir/GWAS/DATASETS/BP/PGC_2021/pgc-bip2021-all.vcf.tsv', header = T, nThread = 23)
bd$N = bd$NCAS + bd$NCON
bd$DIRE = NULL
colnames(bd)[3] = 'SNP'
fwrite(bd, '~/workdir/GWAS/DATASETS/BP/PGC_2021/pgc-bip2021-all_Modified', col.names = T, row.names = F, sep = '\t')

python ~/workdir/programs/ldsc/./munge_sumstats.py \
--sumstats pgc-bip2021-all_Modified \
--merge-alleles ../../w_hm3.snplist \
--chunksize 500000 \
--out pgc-bip2021-all_Modified

# SCZ_2022
# PMID:  35396580, N=76755 + 243649 = 320404

scz = fread('~/workdir/GWAS/DATASETS/SCZ/SCZ_TRUBETSKOY_2022/PGC3_SCZ_wave3.core.autosome.public.v3.vcf.tsv', header = T)
colnames(scz)[2] = 'SNP'
scz$N = scz$NCAS + scz$NCON
scz$DIRE = NULL
fwrite(scz, '~/workdir/GWAS/DATASETS/SCZ/SCZ_TRUBETSKOY_2022/PGC3_SCZ_wave3.core.autosome.public.v3.vcf.tsv_Modified', col.names = T, row.names = F, sep = '\t')

python ~/workdir/programs/ldsc/./munge_sumstats.py \
--sumstats PGC3_SCZ_wave3.core.autosome.public.v3.vcf.tsv_Modified \
--merge-alleles ../../w_hm3.snplist \
--chunksize 500000 \
--out PGC3_SCZ_wave3_Modified

# MDD (2019)
# PMID: 30718901, N=1,306,354 (414,055 + 892,299)

mdd = fread('~/workdir/GWAS/DATASETS/MDD/HOWARD_2019/PGC_UKB_depression_genome-wide.txt', header = T, fill=TRUE)
mdd$Nca = 414055
mdd$Nco = 892299
mdd$N = mdd$Nca + mdd$Nco
colnames(mdd)[1] = 'SNP'
fwrite(mdd, '~/workdir/GWAS/DATASETS/MDD/HOWARD_2019/PGC_UKB_depression_genome-wide_Modified.txt', col.names = T, row.names = F, sep = '\t')

python ~/workdir/programs/ldsc/./munge_sumstats.py \
--sumstats PGC_UKB_depression_genome-wide_Modified.txt \
--merge-alleles ../../w_hm3.snplist \
--signed-sumstats LogOR,0 \
--frq Freq \
--chunksize 500000 \
--out PGC_UKB_depression_genome-wide_Modified


# Anxiety-panic disorder (data:2019, publication:2021)
# PMID:  31712720, N=10240 (2248 + 7992)

anx = fread('~/workdir/GWAS/DATASETS/PANIC_DISORDER/PGC_2019/pgc-panic2019.vcf.tsv', header = T, sep = '\t')
anx$N = anx$NCAS + anx$NCON
colnames(anx)[3] = 'SNP'
fwrite(anx, '~/workdir/GWAS/DATASETS/PANIC_DISORDER/PGC_2019/pgc-panic2019_Modified.tsv', col.names = T, row.names = F, sep = '\t')

python ~/workdir/programs/ldsc/./munge_sumstats.py \
--sumstats pgc-panic2019_Modified.tsv \
--merge-alleles ../../w_hm3.snplist \
--chunksize 500000 \
--out pgc-panic2019_Modified

# Obsessive Compulsive Disorder (data:2017, publication:2018)
# PMID: 28761083, N=10240 (2688 + 7037)

ocd = fread('~/workdir/GWAS/DATASETS/PGC_OCD_Aug2017/ocd_aug2017', header = T, sep = '\t')
ocd$NCAS = 2688
ocd$NCON = 7037
ocd$N = ocd$NCAS + ocd$NCON
fwrite(ocd, '~/workdir/GWAS/DATASETS/PGC_OCD_Aug2017/ocd_aug2017_Modified', col.names = T, row.names = F, sep = '\t')

python ~/workdir/programs/ldsc/./munge_sumstats.py \
--sumstats ocd_aug2017_Modified \
--merge-alleles ../w_hm3.snplist \
--chunksize 500000 \
--out ocd_aug2017_Modified


# AD_2019_Jansen (71,880 cases, 383,378 controls)
# PMID: 30617256, N=455258

ad = fread('~/workdir/GWAS/DATASETS/AD/JANSEN_2019/AD_sumstats_Jansenetal_2019sept.txt', header = T)
ad$N = 455258
fwrite(ad, '~/workdir/GWAS/DATASETS/AD/JANSEN_2019/AD_Jansen2019_Modified', col.names = T, row.names = F, sep = '\t')

python ~/workdir/programs/ldsc/./munge_sumstats.py \
--sumstats AD_Jansen2019_Modified \
--merge-alleles ../../w_hm3.snplist \
--chunksize 500000 \
--signed-sumstats Z,0 \
--out AD_Jansen2019_Modified


# COGNITIVE_FUNCTION
# PMID: 29844566, N=282014

cog = fread('~/workdir/GWAS/DATASETS/COGNITIVE_FUNCTION/Davies2018_OPEN_DATASET_summary_results.txt', header = T, fill=TRUE)
cog$N = 282014
fwrite(cog, '~/workdir/GWAS/DATASETS/COGNITIVE_FUNCTION/Cognitive_Modified', col.names = T, row.names = F, sep = '\t')

python ~/programs/ldsc/./munge_sumstats.py \
--sumstats Cognitive_Modified \
--merge-alleles ../w_hm3.snplist \
--out Cognitive_Modified

mv Cognitive_Modified.sumstats.gz Cognitive_Modified.sumstats


# Intelligence
# PMID: 30038396, N=269867

int = fread('~/workdir/GWAS/DATASETS/Intelligence/SavageJansen_2018_intelligence_metaanalysis.txt', header = T, fill=TRUE)
int$N = 269867
fwrite(int, '~/workdir/GWAS/DATASETS/Intelligence/Intelligence_Modified', col.names = T, row.names = F, sep = '\t')

python ~/programs/ldsc/./munge_sumstats.py \
--sumstats Intelligence_Modified \
--merge-alleles ../w_hm3.snplist \
--out Intelligence_Modified

mv Intelligence_Modified.sumstats.gz Intelligence_Modified.sumstats


# HEIGHT
# https://www.biorxiv.org/content/early/2018/07/09/355057, N=589762

height = fread('~/workdir/GWAS/DATASETS/NON_BRAIN/Meta-analysis_Wood_et_al+UKBiobank_2018.txt', header = T, fill=TRUE)

# N is already added, change column names.
colnames(height)[4:5] = c('A1', 'A2')
fwrite(height, '~/workdir/GWAS/DATASETS/NON_BRAIN/Height_UKBB_Modified', col.names = T, row.names = F, sep = '\t')

python ~/programs/ldsc/./munge_sumstats.py \
--sumstats Height_UKBB_Modified \
--merge-alleles ../w_hm3.snplist \
--out Height_UKBiobank_2018_Modified

mv Height_UKBiobank_2018_Modified.sumstats.gz Height_UKBiobank_2018_Modified.sumstats

# Coronary Artery Disease
# 21378990, N=253288

cad = fread('~/workdir/GWAS/DATASETS/NON_BRAIN/CAD/CARDIoGRAM_GWAS_RESULTS.txt', header = T, fill=TRUE)
cad$N = 84509
fwrite(cad, '~/workdir/GWAS/DATASETS/NON_BRAIN/CAD/CARDIoGRAM_GWAS_Modified', col.names = T, row.names = F, sep = '\t')

python ~/programs/ldsc/./munge_sumstats.py \
--sumstats CARDIoGRAM_GWAS_Modified \
--merge-alleles ../../w_hm3.snplist \
--out CARDIoGRAM_GWAS_Modified

mv CARDIoGRAM_GWAS_Modified.sumstats.gz CARDIoGRAM_GWAS_Modified.sumstats


# Osteoporosis
# 22504420, N=153377

ost = fread('~/workdir/GWAS/DATASETS/NON_BRAIN/OSTEOPOROSIS/GEFOS2_FNBMD_POOLED_GC.txt', header = T, fill=TRUE)
ost$N = 153377
fwrite(ost, '~/workdir/GWAS/DATASETS/NON_BRAIN/OSTEOPOROSIS/GEFOS2_FNBMD_POOLED_Modified', col.names = T, row.names = F, sep = '\t')

python ~/programs/ldsc/./munge_sumstats.py \
--sumstats GEFOS2_FNBMD_POOLED_Modified \
--merge-alleles ../../w_hm3.snplist \
--out GEFOS2_FNBMD_POOLED_Modified

mv GEFOS2_FNBMD_POOLED_Modified.sumstats.gz GEFOS2_FNBMD_POOLED_Modified.sumstats



