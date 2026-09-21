#!/bin/bash
# Copied and modified from /mnt/lvm_vol_1/txu/projects/mwas/chr22.region/process.with.pca.sh
# Modifications:
# M1: CpG sites from selected 30 to all 850k
# M2: SNPs change for each CpG, use SNPs 1M regions flanking the CpG
# Same:
# S1: covariate file, covFile='/mnt/lvm_vol_1/txu/projects/mwas/Methyl_covar_for_Hawaii_ECID.covariates.from.PC'$startPC'.tab';

# if [ "$#" -ne 2 ]; then
#     echo "Usage: $0 raceNumber.like.pAll5.or.from.10.to.14 CpG.list";
if [ "$#" -ne 5 ]; then
    echo "Usage: $0 CpG.list start.PC race.number output.directory covariates.file";
    exit;
fi

# parameters which may change
cpgList=$1;
startPC=$2; # 1, 2, 3, 4, or 5
raceNum=$3; # p10 p11 p12 p13 p14 or pAll5;
opd=$4; # $cwd/$raceNum.fromPC1;
covFile=$5; # '/mnt/lvm_vol_1/txu/projects/mwas/Methyl_covar_for_Hawaii_ECID.'$raceNum'.covariates.from.PC'$startPC'.tab';
cwd=/mnt/lvm_vol_1/txu/projects/mwas/all.850k.CpG.sites;
vcfDataDir=/mnt/lvm_vol_1/txu/projects/mwas/vcfData;
flank=1000000; # 1M up/downstream of each CpG

# Need run plink association analysis for each of 850716 CpG sites against SNPs in +/-1M flanking regions of the given CpG

# for each CpG site,
# A) ########METHYLATION######## need to generate a pheno file which has the methylation level and looks like
# FID IID pheno
# 0 EC002449_EC002449 0.0257587045710443
# 0 EC008620_EC008620 0.0536625027377206
# where column pheno is the methylation level, from files named as All_trans_100K_probes?_header.txt and All_trans_100K_probes?.txt; and IID needs convert between The patient ID in each of the 2110 lines is in the format of 202253620089_R02C01, and it can be converted to EC001952 from file
# /mnt/lvm_vol_1/langwu/MEC_methylation_genetic_data/Covariates/Methyl_covar_for_Hawaii_ECID.csv

# B) ########SNPs######## get SNPs in +/-1M flanking regions of the given CpG

# C) ########COVARIATES######## use covariate file
# /mnt/lvm_vol_1/txu/projects/mwas/vcfData/pAll5.chrA.PCA.analysis/Methyl_covar_for_Hawaii_ECID.covariates.tab has Sex/Age etc
# pcaFile='/mnt/lvm_vol_1/txu/projects/mwas/vcfData/pAll5.chrA.PCA.analysis/pAll5.chrA.ippca1.pca.evec' has PC1~10 from PCA evec file
# mnt/lvm_vol_1/txu/projects/mwas/idolcell.sas7bdat.tab has H_CD8T	H_CD4T	H_NK	H_BCell	H_Mono

# D) ########RACE########
# D1) for each race
# D2) for all race combined, Either SNPs in any race Or SNPs common in all races --- designated as 'any.' and 'com.' in file names


: <<EOF
########METHYLATION########
Methylation levels were parsed from
/mnt/lvm_vol_1/langwu/MEC_methylation_genetic_data/MethylationData/All_trans_100K_probes?.txt
/mnt/lvm_vol_1/langwu/MEC_methylation_genetic_data/MethylationData/All_trans_100K_probes?_header.txt

for i in `seq 1 9`; do ln -s /mnt/lvm_vol_1/langwu/MEC_methylation_genetic_data/MethylationData/All_trans_100K_probes$i.txt; ln -s /mnt/lvm_vol_1/langwu/MEC_methylation_genetic_data/MethylationData/All_trans_100K_probes${i}_header.txt; done
wc -l *header*
 100000 All_trans_100K_probes1_header.txt
 100000 All_trans_100K_probes2_header.txt
 100001 All_trans_100K_probes3_header.txt
 100001 All_trans_100K_probes4_header.txt
 100001 All_trans_100K_probes5_header.txt
 100001 All_trans_100K_probes6_header.txt
 100001 All_trans_100K_probes7_header.txt
 100001 All_trans_100K_probes8_header.txt
  53156 All_trans_100K_probes9_header.txt
 853162 total
Therefore total 853162 CpG sites (~850k CpG sites)

Each header file All_trans_100K_probes?_header.txt has 100000 lines, all have lines like
Methyl_SampleID
cg00000029
cgxxxxxxxx
cgyyyyyyyy
...continue...

while cg00000029 accidentally occurred in all header files.

Each probe file All_trans_100K_probes?.txt have 2110 lines, and each line is a sample/patient
Each line has 100000 columns (comma-delimited), and and each column represents the methylation level of the cgxxxxxxxx ID in the order of all header 
files, with cg00000029 duplicated in all headers (therefore all probes?.txt has the same values in column 2 which is for the first cgxxxxxxxx ID, and
 the first cgxxxxxxxx ID accidentally happen to be cg00000029

The patient ID in each of the 2110 lines is in the format of 202253620089_R02C01, and it can be converted to EC001952 from file
/mnt/lvm_vol_1/langwu/MEC_methylation_genetic_data/Covariates/Methyl_covar_for_Hawaii_ECID.csv

/mnt/lvm_vol_1/txu/projects/mwas/chr22.region/PheT.EWAS.meQTL.tab
Name	CHR	MAPINFO
cg06134331	22	24232799
cg04642813	22	24235569
cg01585852	22	24235823

EOF

# convert beta value to M, using convert.methylation.pheno.beta.to.M.pl
# still using path all.850k.CpG.sites/pheno, and move previous all.850k.CpG.sites to b4.Mvalue.Beta/
: <<EOF
Thanks.

We had a discussion with external colleagues. It seems that we need conduct some updated analyses as:

convert M value from beta value using this formula:

M=log2(beta/(1-beta))

Then please re-conduct the meQTL analyses we did.


>>> Lang Wu 11/08/21 5:56 PM >>>
I mean the beta values of DNA methylation levels not the beta from plink associations. The beta values are from files under /mnt/lvm_vol_1/langwu/MEC_methylation_genetic_data/MethylationDa
ta/

Basically, we need convert beta values to M values for the methylation levels, then rerun plink associations using the M values.
EOF

# cwd=/mnt/lvm_vol_1/txu/projects/mwas/all.850k.CpG.sites;
# opd=$cwd/$raceNum.fromPC1;
# vcfDataDir=/mnt/lvm_vol_1/txu/projects/mwas/vcfData;

# genPheno runBCF mergeRace anyGLM getCommonID comGLM etc are set to 1/0 to turn on/off the corresponding analysis

: <<EOF
########COVARIATES######## 
Covariates file: 3 parts, see get.covariates.pl
Part 1:
53	age at blood draw = AGE_COL
62	sex = SEX
60	Body Mass Index = BMI
61	smoking status = SMKSTAT
75	smoking packyears = PACKYRS
58	total nicotine equivalents = TNE_C

Part 2:
32-41	genetic principal components = PC1 - PC10 (note that there are such PCA infor in the covariate file Methyl_covar_for_Hawaii_ECID.csv generated by MEC programmer, however, in our association analysis we need run PCAs ourselves based on each ethnic group or combined for use. The existing PCAs are expected to generated from analysis of all ethnicity subjects)

Part 3:
126	estimated cell type composition = H_BCELL,
124	H_CD4T,
123	H_CD8T,
127	H_MONO,
125	H_NK

EOF


########COVARIATES########
echo 'beginning to get PC1-10';
# for i in 1 2 3 4 5; do perl get.covariates.pl $i Methyl_covar_for_Hawaii_ECID.covariates.from.PC$i.tab 1>get.covariates.from.PC$i.log 2>get.covariates.from.PC$i.err; echo $i; done
# covFile='/mnt/lvm_vol_1/txu/projects/mwas/Methyl_covar_for_Hawaii_ECID.'$raceNum'.covariates.from.PC'$startPC'.tab';
echo 'end of PC1-10';

########METHYLATION########
# STEP 2, generate pheno file for each CpG
# CpG methylation level from /mnt/lvm_vol_1/txu/projects/mwas/all.850k.CpG.sites/pheno/, which was generated by /mnt/lvm_vol_1/txu/projects/mwas/all.850k.CpG.sites/get.all.CpG.methylation.level.pl
phenoDir=$cwd/pheno;
# The corresponding pheno file is named like pheno/cgxxxxxxxx.methyLevel such as pheno/cg00000776.methyLevel

########SNPs########

echo 'beginning glm';
cpgListFN=`basename $cpgList`;
bimCountFile=$opd/$cpgListFN.bim.count.tab;
# plink2 --pfile hhh.dosage --covar  Methyl_covar_for_Hawaii_ECID.covariates.tab --pheno pheno/cg00024416.methyLevel --glm --out hhh.assoc --covar-variance-standardize
anyGLM=1;
if (( $anyGLM == 1 )); then
#for i in `seq $raceNum $raceNum`; do
    while read line; do
    cpgName=`echo "$line" | cut -f1`;
    cpgChr=`echo "$line" | cut -f2`;
    cpgPos=`echo "$line" | cut -f3`;
    pos1=$cpgPos;
    if (( $cpgPos > $flank )); then
	pos1=`echo $cpgPos - $flank | bc`;
    fi
    pos2=`echo $cpgPos + $flank | bc`;

    pfn=$cwd/pheno/$cpgName.methyLevel;
    if [ -f $pfn ]; then
	rtr=$opd/chr$cpgChr.$cpgName.UD1M;

	# extract SNPs for the CpG site
	if [[ $raceNum == "pAll5" ]]; then 
	    plink --bfile $vcfDataDir/pAll5.chr$cpgChr.pml --chr $cpgChr --from-bp $pos1 --to-bp $pos2 --make-bed --out $rtr # this is for pAll5
	else
	    plink --bfile $vcfDataDir/$raceNum.chr$cpgChr.dose --chr $cpgChr --from-bp $pos1 --to-bp $pos2 --make-bed --out $rtr # this is for ind.race
	fi
	# Need to catch errors like
	# Error: All variants excluded
	# which later generate .assoc. file without ADD line
	# Manually mv these files such as chr9.cg13407658.assoc.ADD.txt to chr9.cg13407658.no.variants file
	if (( $? != 0 )); then # plink make-bed fails
	    rm $rtr.log;
	    echo -n '' > $opd/chr$cpgChr.$cpgName.no.variants;
	else # plink make-bed succeeds
	    lcn=`wc -l $rtr.bim | awk '{print $1}'`;
	    echo -e "$lcn\t$line" >> $bimCountFile;
	# make dosage file
	# plink2 --bfile $rtr --make-pgen --out $rtr.dosage --const-fid && rm $rtr.nosex $rtr.fam $rtr.bed $rtr.bim && gzip $rtr.log;
	plink2 --bfile $rtr --make-pgen --out $rtr.dosage --const-fid && rm $rtr.nosex $rtr.fam $rtr.bed $rtr.bim $rtr.log;

	# plink2 --pfile $rtr.dosage --covar $covFile --pheno $pfn --glm --out $rtr.fromPC$startPC.assoc --covar-variance-standardize && gzip $rtr.dosage.psam $rtr.dosage.pvar $rtr.dosage.pgen $rtr.dosage.log $rtr.fromPC$startPC.assoc.pheno.glm.linear $rtr.fromPC$startPC.assoc.log;
	plink2 --pfile $rtr.dosage --covar $covFile --pheno $pfn --glm --out $rtr.fromPC$startPC.assoc --covar-variance-standardize && rm $rtr.dosage.log && rm $rtr.dosage.psam $rtr.dosage.pvar $rtr.dosage.pgen;

	# clean all other files if association has 'inflation factor for covariate 'PC1' is too high (VIF_TOO_HIGH)' message
	tooHigh=`grep -c 'VIF_TOO_HIGH' $rtr.fromPC$startPC.assoc.log`;
	if (( $tooHigh > 0 )); then
	    rm $rtr.* && echo -n '' > $opd/chr$cpgChr.$cpgName.collinearity.issue;
	else
	    awk '{if($7=="ADD") {print}}' $rtr.fromPC$startPC.assoc.pheno.glm.linear > $opd/chr$cpgChr.$cpgName.assoc.ADD.txt && gzip $opd/chr$cpgChr.$cpgName.assoc.ADD.txt && rm $rtr.fromPC$startPC.assoc.pheno.glm.linear $rtr.fromPC$startPC.assoc.log;
	fi

	fi # if (( $? == 0 ))
    else
	echo -n '' > $opd/chr$cpgChr.$cpgName.no.pheno.file;
    fi
    
    done < $cpgList;
#done
fi

exit
