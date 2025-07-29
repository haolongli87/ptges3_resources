#!/bin/bash
# chipseq_pipeline.sh
# Fully integrated ChIP-seq pipeline

set -euo pipefail

# 0. Set Paths and Parameters

REF_PATH="/fh/fast/nelson_p/whan/ChIPseq_Analyzed_Files"
TRIMMOMATIC_ADAPTERS_PATH="$REF_PATH/trim_adapters"
BWA_INDEX_PATH="$REF_PATH/BWA_index/hg38_bwa"
BLACKLIST="$REF_PATH/Blacklist/GRCh38_unified_blacklist.bed"
HOMER_RUN_PATH="$REF_PATH/homer"
output="/hpc/temp/nelson_p/user/whan/HL"

FASTQ_PATH="$output/fastq"
TEMP_PATH="$output/temp"
SCRIPT_PATH="$output/scripts"
TRIM_PATH="$output/trim"
QC_PATH="$output/qc"
SAM_PATH="$output/sam"
BAM_PATH="$output/bam"
BW_PATH="$output/bw"
MACS3_PATH="$output/peaks"
MERGED_PEAK="$MACS3_PATH/merged"
HOMER_PATH="$output/homer"
HEATMAP_PATH="$output/heatmap"
MATRIX_PATH="$HEATMAP_PATH/matrix"

mkdir -p "$FASTQ_PATH" "$TEMP_PATH" "$SCRIPT_PATH" "$TRIM_PATH" "$QC_PATH" "$SAM_PATH" "$BAM_PATH" "$BW_PATH" "$MACS3_PATH" "$MERGED_PEAK" "$HOMER_PATH" "$HEATMAP_PATH" "$MATRIX_PATH"

# 1. QC and Trimming (FastQC, Trimmomatic)

for file_R1 in "$FASTQ_PATH"/*_R1_001.fastq.gz; do
    sample_name=$(basename "$file_R1" | sed 's/_R1_001.fastq.gz//')
    file_R2="$FASTQ_PATH/${sample_name}_R2_001.fastq.gz"
    [ ! -e "$file_R2" ] && echo "Paired file for $file_R1 not found" && continue
    ml FastQC/0.11.9-Java-11
    fastqc -o "$QC_PATH" "$file_R1" "$file_R2" -d "$TEMP_PATH"
    ml Trimmomatic/0.39-Java-11
    java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -threads 4 "$file_R1" "$file_R2" \
        "$TRIM_PATH/${sample_name}_R1_paired.fastq" "$TRIM_PATH/${sample_name}_R1_unpaired.fastq" \
        "$TRIM_PATH/${sample_name}_R2_paired.fastq" "$TRIM_PATH/${sample_name}_R2_unpaired.fastq" \
        ILLUMINACLIP:$TRIMMOMATIC_ADAPTERS_PATH/TruSeq3-PE.fa:2:30:10:8:true SLIDINGWINDOW:4:20
    fastqc -o "$QC_PATH" "$TRIM_PATH/${sample_name}_R1_paired.fastq" "$TRIM_PATH/${sample_name}_R2_paired.fastq" -d "$TEMP_PATH"
done

# 2. Alignment (BWA-MEM)

ml BWA/0.7.17-GCCcore-12.2.0
ml zlib/1.3.1-GCCcore-13.3.0
ml SAMtools/1.19.2-GCC-13.2.0

for file_R1 in "$TRIM_PATH"/*_R1_paired.fastq; do
    sample_name=$(basename "$file_R1" | sed 's/_R1_paired.fastq//')
    file_R2="$TRIM_PATH/${sample_name}_R2_paired.fastq"
    [ ! -e "$file_R2" ] && echo "Paired file for $file_R1 not found" && continue
    bwa mem -t 8 "$BWA_INDEX_PATH" "$file_R1" "$file_R2" > "$SAM_PATH/$sample_name.sam"
    samtools view -b "$SAM_PATH/$sample_name.sam" -o "$BAM_PATH/$sample_name.bam"
    samtools sort "$BAM_PATH/$sample_name.bam" -o "$BAM_PATH/$sample_name.sorted.bam"
    samtools index "$BAM_PATH/$sample_name.sorted.bam"
done

# 3. Merge Replicates

for rep in 1 2; do
    files_R1=("$BAM_PATH"/*_"${rep}"_1*.sorted.bam)
    for file_R1 in "${files_R1[@]}"; do
        filename_R1=$(basename "$file_R1")
        sample_name="${filename_R1%_${rep}_1*.sorted.bam}"
        file_R2=$(ls "$BAM_PATH/${sample_name}_${rep}_2"*.sorted.bam 2>/dev/null || true)
        [ -z "$file_R2" ] && echo "Paired file for $file_R1 not found" && continue
        samtools merge "$BAM_PATH/${sample_name}_${rep}.merged.bam" "$file_R1" "$file_R2"
    done
done

for file_R1 in "$BAM_PATH"/input_*_1_S*.sorted.bam; do
    filename_R1=$(basename "$file_R1")
    sample_name="${filename_R1%_1_S*.sorted.bam}"
    file_R2=$(ls "$BAM_PATH/${sample_name}_2_S"*.sorted.bam 2>/dev/null || true)
    [ -z "$file_R2" ] && echo "Paired file for $file_R1 not found" && continue
    samtools merge "$BAM_PATH/${sample_name}.merged.bam" "$file_R1" "$file_R2"
    samtools index "$BAM_PATH/${sample_name}.merged.bam"
done

# 4. Filtering, Deduplication, BigWig Generation

ml BEDTools/2.31.0-GCC-12.3.0
ml SAMtools/1.17-GCC-12.2.0
ml picard/2.25.1-Java-11
ml deepTools/3.5.4.post1-gfbf-2022b

for bam_file in "$BAM_PATH"/*.merged.bam; do
    sample_name=$(basename "${bam_file%.bam}")
    bedtools intersect -v -abam "$bam_file" -b "$BLACKLIST" \
    | samtools view -b -q 30 -o "$BAM_PATH/${sample_name}.filtered.q30.bam"
    java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
        I="$BAM_PATH/${sample_name}.filtered.q30.bam" \
        O="$BAM_PATH/${sample_name}.filtered.rmdup.bam" \
        M="$TEMP_PATH/${sample_name}.marked_dup_metrics.txt" \
        REMOVE_SEQUENCING_DUPLICATES=true \
        CREATE_INDEX=true \
        TMP_DIR="$TEMP_PATH" \
        VALIDATION_STRINGENCY=SILENT
    bamCoverage \
        -b "$BAM_PATH/${sample_name}.filtered.rmdup.bam" \
        -o "$BW_PATH/${sample_name}.bw" \
        --binSize 10 \
        --numberOfProcessors 8 \
        --normalizeUsing RPGC \
        --effectiveGenomeSize 2913022398 \
        --outFileFormat 'bigwig'
done

# 5. Peak Calling (MACS3)

ml MACS3/3.0.0-foss-2022b

groups=("sgNTC_1" "sgNTC_2" "sgNTC_R1881_1" "sgNTC_R1881_2" "sgPTGES3_1" "sgPTGES3_2" "sgPTGES3_R1881_1" "sgPTGES3_R1881_2")

for group in "${groups[@]}"; do

    sample_files=($(find "$BAM_PATH" -name "*${group}*.filtered.rmdup.bam" | grep -v "input"))

    for sample_file in "${sample_files[@]}"; do

        case "$group" in
            sgNTC_1) input_file=$(find "$BAM_PATH" -name "input_LNCaP_sgNTC_1*.filtered.rmdup.bam" | head -n 1) ;;
            sgNTC_2) input_file=$(find "$BAM_PATH" -name "input_LNCaP_sgNTC_2*.filtered.rmdup.bam" | head -n 1) ;;
            sgNTC_R1881_1) input_file=$(find "$BAM_PATH" -name "input_LNCaP_sgNTC_R1881_1*.filtered.rmdup.bam" | head -n 1) ;;
            sgNTC_R1881_2) input_file=$(find "$BAM_PATH" -name "input_LNCaP_sgNTC_R1881_2*.filtered.rmdup.bam" | head -n 1) ;;
            sgPTGES3_1) input_file=$(find "$BAM_PATH" -name "input_LNCaP_sgPTGES3_1*.filtered.rmdup.bam" | head -n 1) ;;
            sgPTGES3_2) input_file=$(find "$BAM_PATH" -name "input_LNCaP_sgPTGES3_2*.filtered.rmdup.bam" | head -n 1) ;;
            sgPTGES3_R1881_1) input_file=$(find "$BAM_PATH" -name "input_LNCaP_sgPTGES3_R1881_1*.filtered.rmdup.bam" | head -n 1) ;;
            sgPTGES3_R1881_2) input_file=$(find "$BAM_PATH" -name "input_LNCaP_sgPTGES3_R1881_2*.filtered.rmdup.bam" | head -n 1) ;;
            *) echo "Invalid group $group"; continue ;;
        esac

        if [[ -z "$input_file" ]]; then
            echo "Skipping $sample_file due to missing Input file"
            continue
        fi

        macs3 callpeak -t "$sample_file" -c "$input_file" -q 0.01 -f BAMPE -g hs \
            -n "$MACS3_PATH/$(basename "$sample_file").vs.$(basename "$input_file").q0.01" || {
            echo "MACS3 narrow failed for $(basename "$sample_file") vs $(basename "$input_file")"
            exit 1
        }

        macs3 callpeak -t "$sample_file" -c "$input_file" --broad -g hs --broad-cutoff 0.01 \
            -n "$MACS3_PATH/$(basename "$sample_file").vs.$(basename "$input_file").q0.01.BROAD" || {
            echo "MACS3 broad failed for $(basename "$sample_file") vs $(basename "$input_file")"
            exit 1
        }
    done
done

# 6. Peak Merging and Top Peak Selection

ml BEDTools/2.31.0-GCC-12.3.0
cd "$MACS3_PATH"

cat $(ls *.filtered.rmdup.bam.q0.01_peaks.narrowPeak | grep -Ev 'R1881|merged') | \
  sort -k1,1 -k2,2n | \
  bedtools merge -i - -c 5 -o max | \
  awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "peak_"NR, $4}' > $MERGED_PEAK/All_AR_LNCaP_noR1881.filtered.merged_final.bed

cat $(ls *LNCaP_sgNTC*.filtered.rmdup.bam.q0.01_peaks.narrowPeak | grep -Ev 'R1881|merged') | \
  sort -k1,1 -k2,2n | \
  bedtools merge -i - -c 5 -o max | \
  awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "peak_"NR, $4}' > $MERGED_PEAK/AR_LNCaP_sgNTC.filtered.merged_final.bed
sort -k5,5nr $MERGED_PEAK/AR_LNCaP_sgNTC.filtered.merged_final.bed | head -2000 > $MERGED_PEAK/top2000_peaks_AR_LNCaP_sgNTC.bed

cat $(ls *LNCaP_sgPTGES3*.filtered.rmdup.bam.q0.01_peaks.narrowPeak | grep -Ev 'R1881|merged') | \
  sort -k1,1 -k2,2n | \
  bedtools merge -i - -c 5 -o max | \
  awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "peak_"NR, $4}' > $MERGED_PEAK/AR_LNCaP_sgPTGES3.filtered.merged_final.bed
sort -k5,5nr $MERGED_PEAK/AR_LNCaP_sgPTGES3.filtered.merged_final.bed | head -2000 > $MERGED_PEAK/top2000_peaks_AR_LNCaP_sgPTGES3.bed

# 7. Motif Discovery and Annotation (HOMER)

ml Homer/4.11-Perl-5.30.0

for peak in "top2000_peaks_AR_LNCaP_sgNTC.bed" "top2000_peaks_AR_LNCaP_sgPTGES3.bed"; do
    label=$(basename "$peak" .bed)
    mkdir -p "$HOMER_PATH/$label"
    $HOMER_RUN_PATH/bin/findMotifsGenome.pl "$MERGED_PEAK/$peak" hg38 "$HOMER_PATH/$label/motifs" -size 200 -mask
    $HOMER_RUN_PATH/bin/annotatePeaks.pl "$MERGED_PEAK/$peak" hg38 > "$HOMER_PATH/$label/Annotation_${label}.txt"
done

# 8. Signal Matrix and Heatmap (deepTools)

ml deepTools/3.5.4.post1-gfbf-2022b

cd "$BW_PATH"

computeMatrix reference-point -S \
    AR_LNCaP_sgNTC_1.merged.bw \
    AR_LNCaP_sgNTC_2.merged.bw \
    AR_LNCaP_sgPTGES3_1.merged.bw \
    AR_LNCaP_sgPTGES3_2.merged.bw \
    -R "$PEAK_PATH/merged/All_AR_LNCaP_noR1881.filtered.merged_final.bed" \
    -b 2000 -a 2000 \
    --referencePoint center \
    --outFileName "$HEATMAP_PATH/matrix/AR_LNCaP_merged_output.gz" \
    --outFileSortedRegions "$HEATMAP_PATH/matrix/AR_LNCaP_merged_sorted_regions.bed" \
    --samplesLabel sgNTC_1 sgNTC_2 sgPTGES3_1 sgPTGES3_2 \
    --skipZeros --missingDataAsZero \
    -p 36

plotHeatmap -m "$HEATMAP_PATH/matrix/AR_LNCaP_merged_output.gz" -out "$HEATMAP_PATH/AR_LNCaP_merged_heatmap.pdf" \
    --dpi 500 --plotType lines --colorMap Reds --legendLocation "best" --whatToShow "heatmap and colorbar" --yAxisLabel 'All ARBS'

echo "ChIP-seq pipeline finished successfully!"
