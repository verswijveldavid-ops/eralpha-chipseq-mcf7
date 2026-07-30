#!/usr/bin/env bash
set -euo pipefail

# Storage-conscious raw-read workflow for GSE14664.
# The comparison is E2-treated ERalpha ChIP versus mock-treated ERalpha ChIP.
# It is not ChIP versus input chromatin.

project_dir="${1:-$PWD/work}"
bowtie2_index="${BOWTIE2_INDEX:?Set BOWTIE2_INDEX to a compatible hg18/NCBI36 index}"
gtf_path="${GTF_PATH:?Set GTF_PATH to the matching hg18/NCBI36 GTF annotation}"
chrom_prefix="${CHROM_PREFIX:-chr}"
threads="${THREADS:-8}"

mkdir -p "$project_dir"/{raw_fastq,bam,qc,tracks,peaks,annotation,motifs}
cd "$project_dir"

command -v fasterq-dump fastqc bowtie2 samtools macs3 \
  bamCoverage computeMatrix plotHeatmap bedtools >/dev/null

declare -A accession=(
  [mock]="SRR015349"
  [E2]="SRR015350"
)

for sample in mock E2; do
  run="${accession[$sample]}"
  fasterq-dump "$run" --threads "$threads" --outdir raw_fastq
  mv "raw_fastq/${run}.fastq" "raw_fastq/ERalpha_${sample}.fastq"
  fastqc --threads 2 --outdir qc "raw_fastq/ERalpha_${sample}.fastq"

  bowtie2 \
    -x "$bowtie2_index" \
    -U "raw_fastq/ERalpha_${sample}.fastq" \
    -N 1 -L 20 --no-unal -p "$threads" \
    2> "qc/${sample}_bowtie2.log" |
    samtools sort -@ "$threads" -o "bam/ERalpha_${sample}.sorted.bam"

  samtools index "bam/ERalpha_${sample}.sorted.bam"
  samtools flagstat "bam/ERalpha_${sample}.sorted.bam" \
    > "qc/${sample}_flagstat.txt"

  bamCoverage \
    -b "bam/ERalpha_${sample}.sorted.bam" \
    -o "tracks/ERalpha_${sample}.bw" \
    --normalizeUsing CPM \
    --minMappingQuality 30 \
    --binSize 10 \
    --smoothLength 50 \
    -p "$threads"

  rm "raw_fastq/ERalpha_${sample}.fastq"
done

macs3 callpeak \
  -t bam/ERalpha_E2.sorted.bam \
  -c bam/ERalpha_mock.sorted.bam \
  -n ERalpha_E2_vs_mock \
  --outdir peaks \
  -g hs \
  -q 0.05 \
  --call-summits

computeMatrix reference-point \
  -S tracks/ERalpha_E2.bw tracks/ERalpha_mock.bw \
  -R peaks/ERalpha_E2_vs_mock_summits.bed \
  --referencePoint center \
  -b 2000 -a 2000 \
  --binSize 10 \
  --skipZeros \
  -p "$threads" \
  -o tracks/matrix_peaks.gz

plotHeatmap \
  -m tracks/matrix_peaks.gz \
  -o tracks/heatmap_peaks.png \
  --colorMap Reds Greys \
  --samplesLabel "ERalpha E2" "ERalpha mock" \
  --plotTitle "ERalpha signal at E2-enriched summits (+/-2 kb)" \
  --whatToShow "heatmap and colorbar"

# Ensembl release 54 has no separate gene/transcript features. Reconstruct each
# gene boundary from exon rows and derive its strand-aware transcription start.
awk -v prefix="$chrom_prefix" 'BEGIN{OFS="\t"} $3=="exon" {
    match($0,/gene_id "([^"]+)"/,id)
    match($0,/gene_name "([^"]+)"/,name)
    key=id[1]
    chrom=$1
    if(chrom !~ /^chr/ && prefix!="") chrom=prefix chrom
    if(!(key in min_start) || $4<min_start[key]) min_start[key]=$4
    if(!(key in max_end) || $5>max_end[key]) max_end[key]=$5
    gene_chrom[key]=chrom
    gene_strand[key]=$7
    gene_name[key]=name[1]
  }
  END {
    for(key in min_start) {
      if(gene_strand[key]=="+") print gene_chrom[key],min_start[key]-1,min_start[key],gene_name[key],".",gene_strand[key]
      else print gene_chrom[key],max_end[key]-1,max_end[key],gene_name[key],".",gene_strand[key]
    }
  }' "$gtf_path" |
  sort -k1,1 -k2,2n |
  uniq > annotation/hg18_TSS.bed

awk 'BEGIN{OFS="\t"} {
  midpoint=int(($2+$3)/2)
  print $1,midpoint,midpoint+1,$4,$5,$6
}' peaks/ERalpha_E2_vs_mock_peaks.narrowPeak \
  > annotation/peaks_midpoints.bed

bedtools closest \
  -a annotation/peaks_midpoints.bed \
  -b annotation/hg18_TSS.bed \
  -d > annotation/peaks_TSS_distance.tsv

if command -v findMotifsGenome.pl >/dev/null; then
  awk 'BEGIN{OFS="\t"} {
    start=$2-100
    if(start<0) start=0
    print $1,start,$2+100,$4,$5
  }' peaks/ERalpha_E2_vs_mock_summits.bed \
    > motifs/summits_200bp.bed

  findMotifsGenome.pl \
    motifs/summits_200bp.bed \
    hg19 \
    motifs/homer_output \
    -size given \
    -mask \
    -p "$threads"
fi

du -sh "$project_dir"
