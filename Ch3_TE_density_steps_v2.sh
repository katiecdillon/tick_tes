#!/bin/bash
#SBATCH --job-name=repeatmasker
#SBATCH --partition=iob_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=500gb
#SBATCH --export=NONE
#SBATCH --time=3-00:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-user=kcd88651@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL,ARRAY_TASKS

##########
### STEP 1: Download assemblies
##########
cd /scratch/kcd88651/ticks/ch3/karyoploter/iscap
curl -s https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/920/785/GCF_016920785.2_ASM1692078v2/GCF_016920785.2_ASM1692078v2_genomic.fna.gz | gunzip -c > iscap.fna
cd /scratch/kcd88651/ticks/ch3/karyoploter/iric
curl -s https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/964/199/275/GCA_964199275.3_IXRI_v3/GCA_964199275.3_IXRI_v3_genomic.fna.gz | gunzip -c > iric.fna
cd /scratch/kcd88651/ticks/ch3/karyoploter/dret
curl -s https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/057/332/255/GCA_057332255.1_DR_Devon_UK_G_1.2/GCA_057332255.1_DR_Devon_UK_G_1.2_genomic.fna.gz | gunzip -c > dret.fna

##########
## STEP 2: Unmask assemblies
##########
ml SeqKit/2.9.0

for i in iscap iric dret
do
     cd /scratch/kcd88651/ticks/ch3/karyoploter/${i}
     seqkit seq -u ${i}.fna > ${i}_unmasked.fna
done

module unload SeqKit/2.9.0

##########
## STEP 3: Mask assemblies with manual tick TE library
##########
ml RepeatMasker/4.2.3-foss-2023a

for i in iscap iric dret
do
     cd /scratch/kcd88651/ticks/ch3/karyoploter/${i}
     RepeatMasker -pa 32 \
     -e ncbi \
     -lib /scratch/kcd88651/ticks/NCBI/Dermacentor_reticulatus/UK/repeatmasker_20260806/Elska1004_09252025.fa \
     -xsmall \
     -gff \
     -dir ./ \
     /scratch/kcd88651/ticks/ch3/karyoploter/${i}/${i}_unmasked.fna

module unload RepeatMasker/4.2.3-foss-2023a

##########
## STEP 4 - Index assemblies
##########
ml SAMtools/1.23.1-GCC-13.3.0

for i in iscap iric dret
do
     cd /scratch/kcd88651/ticks/ch3/karyoploter/${i}
     samtools faidx ${i}_unmasked.fna
done

module unload SAMtools/1.23.1-GCC-13.3.0

##########
## STEP 5 - Make chromosome BED files
##########
# D. reticulatus - 11 chromosomes
cd /scratch/kcd88651/ticks/ch3/karyoploter/dret
head -n 11 dret.fna.fai | awk 'BEGIN{OFS="\t"} {print $1, 0, $2}' > dret.genome.bed
awk 'BEGIN{OFS="\t"} {print "chr"NR, $2, $3}' devon.genome.bed > devon.genome.renamed.bed
mv devon.genome.renamed.bed devon.genome.bed ## download to R script directory

# I. ricinus - 14 chromosomes
cd /scratch/kcd88651/ticks/ch3/karyoploter/iric
head -n 14 ${i}.fna.fai | awk 'BEGIN{OFS="\t"} {print $1, 0, $2}' > ${i}.genome.bed
awk 'BEGIN{OFS="\t"} {print "chr"NR, $2, $3}' ${i}.genome.bed > ${i}.genome.renamed.bed
mv ${i}.genome.renamed.bed ${i}.genome.bed ## download to R script directory

# Order I. scapularis scaffolds by length (they don't have chromosome numbers like the other assemblies)
cd /scratch/kcd88651/ticks/ch3/karyoploter/iscap
sort -k2,2 -rn iscap.fna.fai > scaffolds_sorted_by_length.txt
head -n 14 scaffolds_sorted_by_length.txt | awk 'BEGIN{OFS="\t"} {print $1, 0, $2}' > iscap.genome.bed
awk 'BEGIN{OFS="\t"} {print "chr"NR, $2, $3}' iscap.genome.bed > iscap.genome.renamed.bed
mv iscap.genome.renamed.bed iscap.genome.bed ## download to R script directory

##########
## STEP 6 - Filter RepeatMasker.out files
##########
# D. reticulatus - rename chromosomes & filter main TE classes
cd /scratch/kcd88651/ticks/ch3/karyoploter/dret
awk '$5=="CM170536.1"{$5="chr1"} 
     $5=="CM170537.1"{$5="chr2"} 
     $5=="CM170538.1"{$5="chr3"} 
     $5=="CM170539.1"{$5="chr4"} 
     $5=="CM170540.1"{$5="chr5"} 
     $5=="CM170541.1"{$5="chr6"} 
     $5=="CM170542.1"{$5="chr7"} 
     $5=="CM170543.1"{$5="chr8"} 
     $5=="CM170544.1"{$5="chr9"} 
     $5=="CM170545.1"{$5="chr10"} 
     $5=="CM170546.1"{$5="chr11"} 
     $5 ~ /^chr[0-9]+$/ {
        split($11, cls, "/")
        if (cls[1]=="Unknown" || cls[1]=="RC" || cls[1]=="DNA" || 
            cls[1]=="LTR" || cls[1]=="SINE" || cls[1]=="LINE") print
     }' dret_unmasked.fna.out > dret_TEs.out

# D. reticulatus - exclude all "iRic" and "iScap" unknowns
awk '{
    split($11, cls, "/")
    if (cls[1]=="Unknown" && ($10 ~ /^iRic/ || $10 ~ /^iSca/)) {
        next   # skip this line
    }
    print
}' dret_TEs.out > dret_TEs_filtered.out ##download to R folder

# I. ricinus - rename chromosomes & filter main TE classes
cd /scratch/kcd88651/ticks/ch3/karyoploter/iric
awk '$5=="OZ285876.1"{$5="chr1"} 
     $5=="OZ285877.1"{$5="chr2"} 
     $5=="OZ285878.1"{$5="chr3"} 
     $5=="OZ285879.1"{$5="chr4"} 
     $5=="OZ285880.1"{$5="chr5"} 
     $5=="OZ285881.1"{$5="chr6"} 
     $5=="OZ285882.1"{$5="chr7"} 
     $5=="OZ285883.1"{$5="chr8"} 
     $5=="OZ285884.1"{$5="chr9"} 
     $5=="OZ285885.1"{$5="chr10"} 
     $5=="OZ285886.1"{$5="chr11"} 
     $5=="OZ285887.1"{$5="chr12"} 
     $5=="OZ285888.1"{$5="chr13"} 
     $5=="OZ285889.1"{$5="chr14"} 
     $5 ~ /^chr[0-9]+$/ {
        split($11, cls, "/")
        if (cls[1]=="Unknown" || cls[1]=="RC" || cls[1]=="DNA" || 
            cls[1]=="LTR" || cls[1]=="SINE" || cls[1]=="LINE") print
     }' iric_unmasked.fna.out > iric_TEs.out

# I. ricinus - exclude all "iRic" and "iScap" unknowns
awk '{
    split($11, cls, "/")
    if (cls[1]=="Unknown" && ($10 ~ /^iRic/ || $10 ~ /^iSca/)) {
        next   # skip this line
    }
    print
}' iric_TEs.out > iric_TEs_filtered.out ##download to R folder

# I. scapularis - rename chromosomes & filter main TE classes
cd /scratch/kcd88651/ticks/ch3/karyoploter/iscap
awk '$5=="NW_024609835.1"{$5="chr1"} 
     $5=="NW_024609846.1"{$5="chr2"} 
     $5=="NW_024609857.1"{$5="chr3"} 
     $5=="NW_024609868.1"{$5="chr4"} 
     $5=="NW_024609879.1"{$5="chr5"} 
     $5=="NW_024609880.1"{$5="chr6"} 
     $5=="NW_024609881.1"{$5="chr7"} 
     $5=="NW_024609883.1"{$5="chr8"} 
     $5=="NW_024609882.1"{$5="chr9"} 
     $5=="NW_024609884.1"{$5="chr10"} 
     $5=="NW_024609836.1"{$5="chr11"} 
     $5=="NW_024609839.1"{$5="chr12"} 
     $5=="NW_024609837.1"{$5="chr13"} 
     $5=="NW_024609838.1"{$5="chr14"} 
     $5 ~ /^chr[0-9]+$/ {
        split($11, cls, "/")
        if (cls[1]=="Unknown" || cls[1]=="RC" || cls[1]=="DNA" || 
            cls[1]=="LTR" || cls[1]=="SINE" || cls[1]=="LINE") print
     }' iscap_unmasked.fna.out > iscap_TEs.out

# I. scapularis - exclude all "iRic" and "iScap" unknowns
awk '{
    split($11, cls, "/")
    if (cls[1]=="Unknown" && ($10 ~ /^iRic/ || $10 ~ /^iSca/)) {
        next   # skip this line
    }
    print
}' iscap_TEs.out > iscap_TEs_filtered.out ##download to R folder

##########
## STEP 7 - karyoploteR
##########
# Visualize and take measure of each assembly's data using these R scripts:
