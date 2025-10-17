## Raw Fastq Quality Control

To perform quality control on raw FASTQ files using fastp, follow these steps:

### Step 1: Install fastp

Using conda or mamba to install fastp:

```bash
# conda
conda install -c bioconda fastp
# mamba
mamba install -c bioconda fastp
```

### Step 2: Run fastp for Quality Control

For paired-end reads, run fastp with the following command:

```bash
fastp -i sample_R1.fastq.gz -I sample_R2.fastq.gz \
      -o sample_R1.clean.fastq.gz -O sample_R2.clean.fastq.gz \
      -h sample_fastp.html -j sample_fastp.json \
      --thread 8
```

Parameters:

- `-i`: Input R1 FASTQ file
- `-I`: Input R2 FASTQ file (for paired-end)
- `-o`: Output R1 FASTQ file after QC
- `-O`: Output R2 FASTQ file after QC
- `-h`: HTML report file
- `-j`: JSON report file
- `--thread`: Number of threads to use

### Step 3: Review Quality Control Report

fastp automatically performs several QC steps including:

- Adapter trimming
- Quality filtering
- Length filtering
- Per-read quality pruning
- Base correction for overlapped paired-end reads

Open the HTML report (`sample_fastp.html`) to review detailed quality metrics before and after filtering.

### Step 4: Additional Options (Optional)

You can customize fastp with additional parameters:

```bash
fastp -i sample_R1.fastq.gz -I sample_R2.fastq.gz \
      -o sample_R1.clean.fastq.gz -O sample_R2.clean.fastq.gz \
      -h sample_fastp.html -j sample_fastp.json \
      --qualified_quality_phred 20 \
      --length_required 50 \
      --thread 8
```

- `--qualified_quality_phred`: Minimum quality value (default: 15)
- `--length_required`: Minimum read length to keep (default: 15)

Paper: https://academic.oup.com/bioinformatics/article/34/17/i884/5093234