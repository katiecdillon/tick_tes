##########
### STEP 1: Download assemblies
##########
cd /scratch/kcd88651/ticks/ch3/te_proportions/unmasked_genomes
# I. scapularis
curl -s https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/920/785/GCF_016920785.2_ASM1692078v2/GCF_016920785.2_ASM1692078v2_genomic.fna.gz | gunzip -c > utpal.fna
# I. ricinus
curl -s https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/964/199/275/GCA_964199275.3_IXRI_v3/GCA_964199275.3_IXRI_v3_genomic.fna.gz | gunzip -c > ixri.fna
curl -s https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/043/698/115/GCA_043698115.1_Maya_Ir_2.3.1/GCA_043698115.1_Maya_Ir_2.3.1_genomic.fna.gz | gunzip -c > maya.fna
# D. reticulatus
curl -s https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/057/332/255/GCA_057332255.1_DR_Devon_UK_G_1.2/GCA_057332255.1_DR_Devon_UK_G_1.2_genomic.fna.gz | gunzip -c > devon.fna
curl -s https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/051/549/975/GCA_051549975.1_ASM5154997v1/GCA_051549975.1_ASM5154997v1_genomic.fna.gz | gunzip -c > elska.fna

##########
## STEP 2: Unmask assemblies
##########
ml SeqKit/2.9.0

cd /scratch/kcd88651/ticks/ch3/te_proportions/unmasked_genomes
for i in devon elska ixri maya utpal
do
     seqkit seq -u ${i}.fna > ${i}_unmasked.fna
done

module unload SeqKit/2.9.0

##########
## STEP 3: Run RepeatModeler w/LTRStruct
##########
#!/bin/bash
#SBATCH --job-name=repeatmodeler
#SBATCH --partition=highmem_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=900gb
#SBATCH --export=NONE
#SBATCH --time=7-00:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-user=kcd88651@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL,ARRAY_TASKS
#SBATCH --array=1-5

config=/scratch/kcd88651/ticks/home/RepMod_config_20260911.txt
# Array_ID    Genome  Library
# 1   /scratch/kcd88651/ticks/ch3/te_proportions/unmasked_genomes/devon_unmasked.fna  devon_lib
# 2   /scratch/kcd88651/ticks/ch3/te_proportions/unmasked_genomes/elska_unmasked.fna  elska_lib
# 3   /scratch/kcd88651/ticks/ch3/te_proportions/unmasked_genomes/ixri_unmasked.fna   ixri_lib
# 4   /scratch/kcd88651/ticks/ch3/te_proportions/unmasked_genomes/maya_unmasked.fna   maya_lib
# 5   /scratch/kcd88651/ticks/ch3/te_proportions/unmasked_genomes/utpal_unmasked.fna  utpal_lib

unmasked_genome=$(awk -v id=$SLURM_ARRAY_TASK_ID '$1==id {print $2}' "$config")
library=$(awk -v id=$SLURM_ARRAY_TASK_ID '$1==id {print $3}' "$config")

cd /scratch/kcd88651/ticks/ch3/te_proportions/repeatmodeler_output

ml RepeatModeler/2.0.4-foss-2022a

BuildDatabase -name "$library" "$unmasked_genome"

RepeatModeler -database "$library" -threads 32 -LTRStruct

module unload RepeatModeler/2.0.4-foss-2022a

##########
## STEP 4: Combine RepeatModeler libraries
##########
# Rename individual TE libraries and move them to a new directory
cd /scratch/kcd88651/ticks/ch3/te_proportions/repeatmodeler_output
for i in devon elska ixri maya utpal
do
    mv "${i}_lib-families.fa" /scratch/kcd88651/ticks/ch3/te_proportions/te_libraries/"${i}_auto_te_library.fa"
done

# Modify repeat sequence names to reflect which TE library they came from
cd /scratch/kcd88651/ticks/ch3/te_proportions/te_libraries
for i in devon elska ixri maya utpal
do
    sed "s/^>/>${i}_/" "${i}_auto_te_library.fa" > "${i}_auto_te_library_prefixed.fa"
done

# Combine all the prefixed libraries into one automated TE library
cat devon_auto_te_library_prefixed.fa elska_auto_te_library_prefixed.fa ixri_auto_te_library_prefixed.fa maya_auto_te_library_prefixed.fa utpal_auto_te_library_prefixed.fa > combined_auto_te_library.fa

##########
## STEP 5: Run RepeatMasker w/combined automated TE library
##########
#!/bin/bash
#SBATCH --job-name=RepMask
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
#SBATCH --array=1-5

config=/scratch/kcd88651/ticks/home/ReMa_config_20260911.txt
# Array_ID    Genome  Library
# 1   /scratch/kcd88651/ticks/ch3/te_proportions/unmasked_genomes/devon_unmasked.fna  devon_lib
# 2   /scratch/kcd88651/ticks/ch3/te_proportions/unmasked_genomes/elska_unmasked.fna  elska_lib
# 3   /scratch/kcd88651/ticks/ch3/te_proportions/unmasked_genomes/ixri_unmasked.fna   ixri_lib
# 4   /scratch/kcd88651/ticks/ch3/te_proportions/unmasked_genomes/maya_unmasked.fna   maya_lib
# 5   /scratch/kcd88651/ticks/ch3/te_proportions/unmasked_genomes/utpal_unmasked.fna  utpal_lib
unmasked_genome=$(awk -v id=$SLURM_ARRAY_TASK_ID '$1==id {print $2}' "$config")
RMDIR=/scratch/kcd88651/ticks/ch3/te_proportions/repeatmasker_output/automated

cd "$RMDIR"

ml RepeatMasker/4.2.3-foss-2023a

RepeatMasker -pa 32 \
-e ncbi \
-lib /scratch/kcd88651/ticks/ch3/te_proportions/te_libraries/combined_auto_te_library.fa \
-xsmall \
-gff \
-a \
-dir ./ \
"$unmasked_genome"

module unload RepeatMasker/4.2.3-foss-2023a

##########
## STEP 6: Run RepeatMasker w/combined manual TE library
##########
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
#SBATCH --array=1-5

config=/scratch/kcd88651/ticks/home/ReMa_config_20260911.txt
# Array_ID    Genome  Library
# 1   /scratch/kcd88651/ticks/ch3/te_proportions/unmasked_genomes/devon_unmasked.fna  devon_lib
# 2   /scratch/kcd88651/ticks/ch3/te_proportions/unmasked_genomes/elska_unmasked.fna  elska_lib
# 3   /scratch/kcd88651/ticks/ch3/te_proportions/unmasked_genomes/ixri_unmasked.fna   ixri_lib
# 4   /scratch/kcd88651/ticks/ch3/te_proportions/unmasked_genomes/maya_unmasked.fna   maya_lib
# 5   /scratch/kcd88651/ticks/ch3/te_proportions/unmasked_genomes/utpal_unmasked.fna  utpal_lib
unmasked_genome=$(awk -v id=$SLURM_ARRAY_TASK_ID '$1==id {print $2}' "$config")
RMDIR=/scratch/kcd88651/ticks/ch3/te_proportions/repeatmasker_output/manual

cd "$RMDIR"

ml RepeatMasker/4.2.3-foss-2023a

RepeatMasker -pa 32 \
-e ncbi \
-lib /scratch/kcd88651/ticks/ch3/te_proportions/te_libraries/manual_te_library.fa \
-xsmall \
-gff \
-a \
-dir ./ \
"$unmasked_genome"

module unload RepeatMasker/4.2.3-foss-2023a

##########
## STEP 7: Filter RepeatMasker.out files (MANUAL ONLY!!!!)
##########
## Filter the .out file to remove iRic & iScap unknowns (do this for MANUAL ONLY!!!)
for k in devon elska ixri maya utpal
do
    awk '{
        split($11, cls, "/")
        if (cls[1]=="Unknown" && ($10 ~ /^iRic/ || $10 ~ /^iSca/)) {
            next   # skip this line
        }
        print
    }' ${k}_unmasked.fna.out > ${k}_filtered.out
done

##########
## STEP 8: Run te_plotting2.py (automated & manual)
##########
#!/bin/bash
#SBATCH --job-name=te_plotting
#SBATCH --partition=iob_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=100gb
#SBATCH --export=NONE
#SBATCH --time=1-00:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-user=kcd88651@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL,ARRAY_TASKS
#SBATCH --array=1-10

source ~/.bashrc
conda activate te_plotting
config=/scratch/kcd88651/ticks/home/teplotting2_config_20260929.txt
# Array_ID	genome	library	rm_out	rm_tbl
# 1	devon	manual	/scratch/kcd88651/ticks/ch3/te_proportions/repeatmasker_output/manual/devon_filtered.out	/scratch/kcd88651/ticks/ch3/te_proportions/repeatmasker_output/manual/devon_unmasked.fna.tbl
# 2	elska	manual	/scratch/kcd88651/ticks/ch3/te_proportions/repeatmasker_output/manual/elska_filtered.out	/scratch/kcd88651/ticks/ch3/te_proportions/repeatmasker_output/manual/elska_unmasked.fna.tbl
# 3	ixri	manual	/scratch/kcd88651/ticks/ch3/te_proportions/repeatmasker_output/manual/ixri_filtered.out /scratch/kcd88651/ticks/ch3/te_proportions/repeatmasker_output/manual/ixri_unmasked.fna.tbl
# 4	maya	manual	/scratch/kcd88651/ticks/ch3/te_proportions/repeatmasker_output/manual/maya_filtered.out /scratch/kcd88651/ticks/ch3/te_proportions/repeatmasker_output/manual/maya_unmasked.fna.tbl
# 5	utpal	manual	/scratch/kcd88651/ticks/ch3/te_proportions/repeatmasker_output/manual/utpal_filtered.out    /scratch/kcd88651/ticks/ch3/te_proportions/repeatmasker_output/manual/utpal_unmasked.fna.tbl
# 6	devon	automated	/scratch/kcd88651/ticks/ch3/te_proportions/repeatmasker_output/automated/devon_unmasked.fna.out	/scratch/kcd88651/ticks/ch3/te_proportions/repeatmasker_output/automated/devon_unmasked.fna.tbl
# 7	elska	automated	/scratch/kcd88651/ticks/ch3/te_proportions/repeatmasker_output/automated/elska_unmasked.fna.out	/scratch/kcd88651/ticks/ch3/te_proportions/repeatmasker_output/automated/elska_unmasked.fna.tbl
# 8	ixri	automated	/scratch/kcd88651/ticks/ch3/te_proportions/repeatmasker_output/automated/ixri_unmasked.fna.out	/scratch/kcd88651/ticks/ch3/te_proportions/repeatmasker_output/automated/ixri_unmasked.fna.tbl
# 9	maya	automated	/scratch/kcd88651/ticks/ch3/te_proportions/repeatmasker_output/automated/maya_unmasked.fna.out	/scratch/kcd88651/ticks/ch3/te_proportions/repeatmasker_output/automated/maya_unmasked.fna.tbl
# 10	utpal	automated	/scratch/kcd88651/ticks/ch3/te_proportions/repeatmasker_output/automated/utpal_unmasked.fna.out	/scratch/kcd88651/ticks/ch3/te_proportions/repeatmasker_output/automated/utpal_unmasked.fna.tbl

GENOME=$(awk -v id=$SLURM_ARRAY_TASK_ID '$1==id {print $2}' "$config")
LIBRARY=$(awk -v id=$SLURM_ARRAY_TASK_ID '$1==id {print $3}' "$config")
RM_OUT=$(awk -v id=$SLURM_ARRAY_TASK_ID '$1==id {print $4}' "$config")
RM_TBL=$(awk -v id=$SLURM_ARRAY_TASK_ID '$1==id {print $5}' "$config")

SCRIPT=/scratch/kcd88651/ticks/ch3/te_proportions/te_plotting2/te_plotting2_v2.py
OUTDIR=/scratch/kcd88651/ticks/ch3/te_proportions/te_plotting2/${LIBRARY}

cd "$OUTDIR"

python3 "$SCRIPT" \
    -r "$RM_OUT" \
    -s "$RM_TBL" \
    -op "$GENOME" \
    -c LINE,SINE,LTR,DNA,RC,Unknown

conda deactivate

##########
## STEP 9: Combine te_plotting2.py output (automated & manual)
##########
# combine pie data
cd /scratch/kcd88651/ticks/ch3/te_proportions/te_plotting2
combined="/scratch/kcd88651/ticks/ch3/te_proportions/te_plotting2/all_genomes_pie_data_combined.tsv"
printf "Genome\tTE_approach\tClass\tProportion\n" > "$combined"

for lib in manual automated
do
    case "$lib" in
        manual)    label="manual" ;;
        automated) label="auto" ;;
    esac

    dir="/scratch/kcd88651/ticks/ch3/te_proportions/te_plotting2/${lib}"

    for k in devon elska ixri maya utpal
    do
        pie_file="${dir}/${k}_pie_data.tsv"
        if [ ! -f "$pie_file" ]; then
            echo "WARNING: Missing $pie_file — skipping" >&2
            continue
        fi
        tail -n +2 "$pie_file" | awk -v gen="$k" -v appr="$label" 'BEGIN{OFS="\t"} {print gen, appr, $0}' >> "$combined"
    done
done
# download combined tsv

# combine divergence data
combined="/scratch/kcd88651/ticks/ch3/te_proportions/te_plotting2/all_genomes_stacked_bar_data_combined.tsv"
printf "Genome\tTE_approach\tBin\tDNA\tLINE\tLTR\tRC\tSINE\tUnknown\n" > "$combined"

for lib in manual automated
do
    case "$lib" in
        manual)    label="manual" ;;
        automated) label="auto" ;;
    esac

    dir="/scratch/kcd88651/ticks/ch3/te_proportions/te_plotting2/${lib}"

    for k in devon elska ixri maya utpal
    do
        bar_file="${dir}/${k}_stacked_bar_data.tsv"
        if [ ! -f "$bar_file" ]; then
            echo "WARNING: Missing $bar_file — skipping" >&2
            continue
        fi

        awk -v gen="$k" -v appr="$label" '
        BEGIN {FS=OFS="\t"; split("Bin,DNA,LINE,LTR,RC,SINE,Unknown", target, ",")}
        NR==1 {
            for (i=1; i<=NF; i++) colpos[$i]=i
            for (t=1; t<=7; t++) {
                if (!(target[t] in colpos)) {
                    print "WARNING: column " target[t] " missing in " FILENAME > "/dev/stderr"
                }
            }
            next
        }
        {
            printf "%s\t%s", gen, appr
            for (t=1; t<=7; t++) {
                col = target[t]
                if (col in colpos) {
                    val = $(colpos[col])
                } else {
                    val = 0
                }
                printf "\t%s", val
            }
            printf "\n"
        }' "$bar_file" >> "$combined"
    done
done
# download combined tsv


##########
## STEP 10: Prepare input files for bedtools intersect
##########
## Select for only the major TE classes (AUTOMATED ONLY!!!)
cd /scratch/kcd88651/ticks/ch3/te_proportions/repeatmasker_output/automated
for i in devon elska ixri maya utpal
do
    awk 'NR<=3 {print; next}
        {
            split($11, cf, "/")
            if (cf[1]=="LINE" || cf[1]=="SINE" || cf[1]=="LTR" || cf[1]=="DNA" || cf[1]=="RC" || cf[1]=="Unknown")
                print
        }' ${i}_unmasked.fna.out > ${i}_TEs.out

    python3 /scratch/kcd88651/ticks/ch3/te_proportions/repeatmasker_output/process_repeatmasker_v8.py \
        -i ${i}_TEs.out \
        -ot bed \
        -p ${i}_overlap \
        -ov longer_element
done

## Select for only the major TE classes (MANUAL ONLY!!!)
cd /scratch/kcd88651/ticks/ch3/te_proportions/repeatmasker_output/manual
for i in devon elska ixri maya utpal
do
    awk 'NR<=3 {print; next}
        {
            split($11, cf, "/")
            if (cf[1]=="LINE" || cf[1]=="SINE" || cf[1]=="LTR" || cf[1]=="DNA" || cf[1]=="RC" || cf[1]=="Unknown")
                print
        }' ${i}_filtered.out > ${i}_filtered_TEs.out

    python3 /scratch/kcd88651/ticks/ch3/te_proportions/repeatmasker_output/process_repeatmasker_v8.py \
        -i ${i}_filtered_TEs.out \
        -ot bed \
        -p ${i}_overlap \
        -ov longer_element
    
done

##########
## STEP 11: Run bedtools intersect
##########
#!/bin/bash
#SBATCH --job-name=intersect
#SBATCH --partition=iob_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=500gb
#SBATCH --export=NONE
#SBATCH --time=1-00:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-user=kcd88651@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL,ARRAY_TASKS

ml BEDTools/2.31.1-GCC-13.3.0

for k in devon elska ixri maya utpal
do
    MANUAL=/scratch/kcd88651/ticks/ch3/te_proportions/repeatmasker_output/manual/${k}_overlap.bed
    AUTO=/scratch/kcd88651/ticks/ch3/te_proportions/repeatmasker_output/automated/${k}_overlap.bed
    TBL=/scratch/kcd88651/ticks/ch3/te_proportions/repeatmasker_output/manual/${k}_unmasked.fna.tbl
    OUTDIR=/scratch/kcd88651/ticks/ch3/te_proportions/intersect/${k}

    mkdir -p "$OUTDIR"
    cd "$OUTDIR"

    GENOMESIZE=$(grep -i "total length:" "$TBL" | awk '{print $3}')

    # Split combined Class/Family (column 7) into two separate columns before sorting,
    # so every downstream file already has clean, separate Class and Family fields
    MANUAL_SPLIT="${k}_manual_split.bed"
    AUTO_SPLIT="${k}_auto_split.bed"

    awk -F'\t' 'BEGIN{OFS="\t"} {
        split($7, cf, "/")
        fam = (cf[2]=="") ? cf[1] : cf[2]
        print $1, $2, $3, $4, $5, $6, cf[1], fam
    }' "$MANUAL" > "$MANUAL_SPLIT"

    awk -F'\t' 'BEGIN{OFS="\t"} {
        split($7, cf, "/")
        fam = (cf[2]=="") ? cf[1] : cf[2]
        print $1, $2, $3, $4, $5, $6, cf[1], fam
    }' "$AUTO" > "$AUTO_SPLIT"

    MANUAL_SORTED="${k}_manual_sorted.bed"
    AUTO_SORTED="${k}_auto_sorted.bed"
    sort -k1,1 -k2,2n "$MANUAL_SPLIT" > "$MANUAL_SORTED"
    sort -k1,1 -k2,2n "$AUTO_SPLIT" > "$AUTO_SORTED"

    bedtools intersect -v -a "$MANUAL_SORTED" -b "$AUTO_SORTED" > ${k}_manual_only.bed
    bedtools intersect -v -a "$AUTO_SORTED" -b "$MANUAL_SORTED" > ${k}_auto_only.bed
    bedtools intersect -wo -a "$MANUAL_SORTED" -b "$AUTO_SORTED" > ${k}_manual_vs_auto_overlap.tsv
    bedtools jaccard -a "$MANUAL_SORTED" -b "$AUTO_SORTED" > ${k}_jaccard.txt

    printf "Count\tManual_Class\tAutomated_Class\n" > ${k}_class_agreement_matrix.txt
    awk -F'\t' 'BEGIN{OFS="\t"} {print $7, $15}' ${k}_manual_vs_auto_overlap.tsv | sort | uniq -c | \
        awk 'BEGIN{OFS="\t"} {print $1, $2, $3}' >> ${k}_class_agreement_matrix.txt

    printf "Manual_Class\tAutomated_Class\tOverlap_bp\tPercent_of_Genome\n" > ${k}_class_agreement_bp.txt
    awk -F'\t' -v gsize="$GENOMESIZE" '{sum[$7"\t"$15]+=$17} END {
        for (key in sum) printf "%s\t%d\t%.4f\n", key, sum[key], (sum[key]/gsize)*100
    }' ${k}_manual_vs_auto_overlap.tsv | sort -t$'\t' -k3 -rn >> ${k}_class_agreement_bp.txt
done

# Combine all the class agreement files into one .txt file
cd /scratch/kcd88651/ticks/ch3/te_proportions/intersect
combined="all_genomes_class_agreement_bp.txt"
printf "Genome\tManual_Class\tAutomated_Class\tOverlap_Mbp\tPercent_of_Genome\n" > "$combined"

for k in devon elska ixri maya utpal
do
    file="/scratch/kcd88651/ticks/ch3/te_proportions/intersect/${k}/${k}_class_agreement_bp.txt"
    if [ ! -f "$file" ]; then
        echo "WARNING: Missing $file — skipping" >&2
        continue
    fi
    tail -n +2 "$file" | awk -v gen="$k" 'BEGIN{OFS="\t"} {
        mbp = $3 / 1000000
        printf "%s\t%s\t%s\t%.4f\t%s\n", gen, $1, $2, mbp, $4
    }' >> "$combined"
done # download to R directory