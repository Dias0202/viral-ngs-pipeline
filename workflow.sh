#!/bin/bash

set -e

# === VERBOSE COLORS  ===
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m'

echo -e "${YELLOW}[INFO] Starting NGS pipeline...${NC}"

# === DEPENDENCY CHECK ===
check_dep() {
    if ! command -v $1 &> /dev/null; then
        echo -e "${YELLOW}[INFO] Installing $1...${NC}"
        conda install -y -c bioconda $1 || conda install -y -c conda-forge $1
    else
        echo -e "${GREEN}[✔] $1 available.${NC}"
    fi
}

# Main dependencies
check_dep fastqc
check_dep cutadapt
check_dep spades.py
check_dep bwa
check_dep samtools
check_dep ivar
check_dep python

# Python dependencies for plotting
python - <<END
try:
    import pandas, matplotlib
except ImportError:
    import os
    os.system('pip install pandas matplotlib')
END

# === Paths ===
REF="reference/reference.fasta"
PRIMERS="primer.fasta"

# === Directory creation ===
mkdir -p fastqc_raw fastqc_trimmed trimmed spades bams consensus quast coverage_plots logs

# === BWA indexing (if needed) ===
if [ ! -f "${REF}.bwt" ]; then
    echo -e "${YELLOW}[INFO] Indexing reference genome...${NC}"
    bwa index $REF
else
    echo -e "${GREEN}[✔] BWA index already exists.${NC}"
fi

# === Detect samples ===
samples=($(ls fastq/*_R1_001.fastq.gz | sed 's/.*\///; s/_R1_001\.fastq\.gz//' | sort -u))

# === Initialize QUAST summary file ===
echo -e "Sample\tContigs\tN50" > quast_summary.tsv

# === Main loop ===
for sample in "${samples[@]}"; do
    echo -e "${YELLOW}[INFO] Processing sample: $sample${NC}"

    R1="fastq/${sample}_R1_001.fastq.gz"
    R2="fastq/${sample}_R2_001.fastq.gz"

    echo -e "${YELLOW}[1/8] Raw FastQC...${NC}"
    fastqc -o fastqc_raw "$R1" "$R2"

    echo -e "${YELLOW}[2/8] Cutadapt (primer removal)...${NC}"
    cutadapt -g file:"$PRIMERS" -G file:"$PRIMERS" \
        -o trimmed/${sample}_R1_trimmed.fastq.gz \
        -p trimmed/${sample}_R2_trimmed.fastq.gz \
        "$R1" "$R2"

    echo -e "${YELLOW}[3/8] Post-trimming FastQC...${NC}"
    fastqc -o fastqc_trimmed trimmed/${sample}_R1_trimmed.fastq.gz trimmed/${sample}_R2_trimmed.fastq.gz

    echo -e "${YELLOW}[4/8] Assembly with SPAdes...${NC}"
    spades.py --careful \
        -1 trimmed/${sample}_R1_trimmed.fastq.gz \
        -2 trimmed/${sample}_R2_trimmed.fastq.gz \
        -o spades/$sample

    echo -e "${YELLOW}[5/8] Mapping with BWA + BAM generation...${NC}"
    bwa mem -t 4 "$REF" trimmed/${sample}_R1_trimmed.fastq.gz trimmed/${sample}_R2_trimmed.fastq.gz | \
        samtools sort -@ 4 -o bams/${sample}.bam
    samtools index bams/${sample}.bam

    echo -e "${YELLOW}[6/8] Consensus generation with iVar...${NC}"
    samtools mpileup -A -d 0 --reference "$REF" -B -Q 0 bams/${sample}.bam | \
        ivar consensus -p consensus/${sample}_consensus -q 20 -t 0.6 -m 10 -n N

    echo -e "${YELLOW}[7/8] Evaluation with QUAST...${NC}"
	QUAST_PATH=$(find "$CONDA_PREFIX" -name "quast-lg.py" 2>/dev/null | head -n 1)
	if [ -z "$QUAST_PATH" ]; then
        echo -e "${RED}[ERROR] Could not locate quast-lg.py. Make sure the 'quast' environment is activated and QUAST is installed.${NC}"
	else
        export PATH=$(dirname "$QUAST_PATH"):$PATH
        quast-lg.py "spades/$sample/scaffolds.fasta" \
            -r "$REF" \
            -o "quast/$sample" \
            --threads 4 \
            --min-contig 500 \
            --min-alignment 100 \
            &> logs/${sample}_quast.log \
            || echo -e "${YELLOW}[WARNING] QUAST could not align contigs. Check logs in quast/$sample.${NC}"
        
        contigs=$(grep ">" spades/$sample/scaffolds.fasta | wc -l)
        N50=$(grep "N50" quast/$sample/report.txt | awk '{print $2}')
        echo -e "$sample\t$contigs\t$N50" >> quast_summary.tsv
		echo -e "${YELLOW}[QUAST summary available at quast_summary.tsv]${NC}"
	fi

    echo -e "${YELLOW}[8/8] Coverage plot generation...${NC}"
    samtools depth bams/${sample}.bam > coverage_plots/${sample}_coverage.tsv

    python3 - <<EOF
import pandas as pd
import matplotlib.pyplot as plt

try:
    df = pd.read_csv("coverage_plots/${sample}_coverage.tsv", sep="\t", header=None, names=["chrom", "position", "depth"])
    plt.figure(figsize=(12,4))
    plt.plot(df["position"], df["depth"], color='darkblue')
    plt.title("Coverage - $sample")
    plt.xlabel("Position (bp)")
    plt.ylabel("Depth")
    plt.tight_layout()
    plt.savefig("coverage_plots/${sample}_coverage.png", dpi=300)
except Exception as e:
    print(f"[ERROR] Failed to generate plot for {sample}: {e}")
EOF

    echo -e "${GREEN}[✔] Sample $sample processed successfully!${NC}"
done

echo -e "${GREEN}[✔] Workflow completed successfully!${NC}"

