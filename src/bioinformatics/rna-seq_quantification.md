## Quantification

To quantify gene expression from RNA-seq alignments using featureCounts, follow these steps:

### Step 1: Install featureCounts

featureCounts is part of the Subread package. Install it using conda or mamba:

```bash
# conda
conda install -c bioconda subread
# mamba
mamba install -c bioconda subread
```

### Step 2: Run featureCounts

Use featureCounts to count reads mapped to genomic features (genes) from your BAM file:

```bash
featureCounts -a /path/to/annotations.gtf \
    -o gene_counts.txt \
    -t exon \
    -g gene_id \
    -T 8 \
    -p \
    output_prefix_Aligned.sortedByCoord.out.bam
```

Parameters:

- `-a`: Path to the gene annotation GTF file
- `-o`: Output file name for the count matrix
- `-t`: Feature type to count (default: exon)
- `-g`: Attribute type in GTF file to use as feature ID (default: gene_id)
- `-T`: Number of threads to use
- `-p`: Specify that input data is paired-end

### Step 3: Additional Options

You can customize featureCounts with additional parameters for strand-specific RNA-seq:

```bash
featureCounts -a /path/to/annotations.gtf \
    -o gene_counts.txt \
    -t exon \
    -g gene_id \
    -T 8 \
    -p \
    -s 2 \
    output_prefix_Aligned.sortedByCoord.out.bam
```

- `-s`: Specify strand-specificity (0: unstranded, 1: stranded, 2: reversely stranded)
- `-M`: Count multi-mapping reads
- `-O`: Assign reads to overlapping features
- `--fraction`: Assign fractional counts to features

### Step 4: Output Files

featureCounts generates two main output files:

- `gene_counts.txt`: Tab-delimited file with gene counts per sample
- `gene_counts.txt.summary`: Summary statistics including assigned, unassigned, and ambiguous reads

### Step 5: Review Count Results

The output count file contains columns for:

- Geneid: Gene identifier
- Chr: Chromosome
- Start: Gene start position
- End: Gene end position
- Strand: Gene strand
- Length: Gene length
- Sample columns: Read counts for each BAM file

Review the summary file to check the proportion of successfully assigned reads and identify potential issues with your alignment or annotation.

Paper: https://academic.oup.com/bioinformatics/article/30/7/923/232889