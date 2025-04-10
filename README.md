# Viral Genome Assembly Pipeline

This pipeline performs a full analysis of NGS paired-end sequencing data (Illumina) for viral genome assembly, including quality control, primer trimming, de novo assembly, reference mapping, consensus generation, and assembly evaluation.

---

## 🛠 Requirements

- **Operating system**: Linux
- **Conda**: Miniconda3 recommended

Required tools (via Conda):

- `fastqc`
- `cutadapt`
- `spades`
- `bwa`
- `samtools`
- `ivar`
- `quast`
- `python` (>=3.7) with `pandas` and `matplotlib`

### Create and activate the environment:

```bash
conda create -n viral-ngs -c bioconda -c conda-forge fastqc cutadapt spades bwa samtools ivar quast python=3.10
conda activate viral-ngs
pip install pandas matplotlib
```

---

## 📁 Expected Folder Structure

Before running the pipeline, organize your files as follows:

```
├── fastq/                      # Raw sequencing files (.fastq.gz)
│   ├── sample1_R1_001.fastq.gz
│   └── sample1_R2_001.fastq.gz
├── reference/                 # Reference genome for mapping and QUAST
│   └── reference.fasta
├── primer.fasta               # FASTA file containing primers to be trimmed
├── workflow.sh                # Main script
```
📌 Creating the primer.fasta file
The primer.fasta file must contain all primer sequences used during library preparation, in standard FASTA format. Each primer must start with a header line beginning with >, followed by the primer sequence on the next line.

```
>primer_F1
ACTGTTGGAAGGTTTGTAGG
>primer_R1
TCCATGAGGTTCTTGAGTGT
>primer_F2
GGTGACAGTTGAGAGACACC
>primer_R2
AACCATCAGGTGTTCTCTG
```

---

## ▶️ How to Run

Make the script executable:

```bash
chmod +x workflow.sh
```

Then simply run:

```bash
./workflow.sh
```

Samples are detected automatically based on the paired-end file names in the `fastq/` folder.

---

## 🔄 Pipeline Steps

For each sample, the pipeline executes:

1. **FastQC** – Raw read quality control.
2. **Cutadapt** – Primer removal using sequences from `primer.fasta`.
3. **FastQC** – Post-trimming quality control.
4. **SPAdes** – Error-corrected de novo assembly (`--careful` mode).
5. **BWA + Samtools** – Mapping trimmed reads to reference genome.
6. **iVar** – Generation of consensus genome from mapped reads.
7. **QUAST** – Assembly metrics (contigs, N50, GC content, length).
8. **Matplotlib** – Coverage plot (read depth per position).

---

## 📦 Output Structure

Results are automatically organized into:

```
├── bams/                    # BAM files (mapped and indexed)
├── consensus/              # Consensus genome FASTA files
├── coverage_plots/         # PNG plots + depth tables (TSV)
├── fastqc_raw/             # Raw FastQC reports
├── fastqc_trimmed/         # Trimmed FastQC reports
├── logs/                   # QUAST logs
├── quast/                  # QUAST reports for each sample
├── quast_summary.tsv       # QUAST summary (N50, contigs etc.)
├── spades/                 # SPAdes assemblies per sample
├── trimmed/                # Trimmed reads (after primer removal)
```

---

## ✅ Tips & Recommendations

- Optimized for **short viral genomes** (<30 kb).
- Adjust SPAdes and iVar parameters if you want more/less conservative assemblies.
- If using a different reference genome, replace `reference.fasta` (keep the same name or edit the script), and BWA index will be generated automatically if missing.
- Coverage plots are based on reference alignment, not contig scaffolds.

---

## 📬 
Feel free to fork, open issues, or suggest improvements!

---
