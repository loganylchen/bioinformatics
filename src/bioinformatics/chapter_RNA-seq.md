## 2.1 Raw Fastq Quality Control

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

## 2.2 RNA-seq Mapping

To map paired-end (PE) FASTQ reads to a reference genome using STAR (Spliced Transcripts Alignment to a Reference) (Paper:  https://academic.oup.com/bioinformatics/article/29/1/15/272537?login=false), follow these steps:

### Step 1: Install STAR

Using conda or mamba to install STAR is your best choice:

```bash
# conda
conda install -c bioconda star
# mamba
mamba install -c bioconda star
```

### Step 2: Generate Genome Index

Before mapping, you need to generate a genome index. This step is done once per reference genome. (You should know this step if you encounter similar tasks in the future, but this time, I will build the index for you.)

```bash
STAR --runMode genomeGenerate \
     --genomeDir /path/to/genomeDir \
     --genomeFastaFiles /path/to/genome.fa \
     --sjdbGTFfile /path/to/annotations.gtf \
     --runThreadN 8
     
# Above is the general command, I will also show you how I run the index building
cd /home/shared/reference # move into the reference folder
STAR --runMode genomeGenerate \
		--genomeDir STAR_INDEX \
		--genomeFastaFiles ./Homo_sapiens.GRCh38.dna.primary_assembly.fa \
		--sjdbGTFfile ./Homo_sapiens.GRCh38.104.gtf \
		--runThreadN 8
```

Parameters:

- `--genomeDir`: Directory where the genome index will be stored
- `--genomeFastaFiles`: Path to the reference genome FASTA file
- `--sjdbGTFfile`: Path to the gene annotation GTF file
- `--runThreadN`: Number of threads to use

So the index is in `/home/shared/reference/STAR_INDEX`, you can use it directly in the later steps.

### Step 3: Map Paired-End Reads

Run STAR to align your paired-end FASTQ files to the indexed genome:

```bash
STAR --genomeDir /path/to/genomeDir \
     --readFilesIn sample_R1.fastq.gz sample_R2.fastq.gz \
     --readFilesCommand zcat \
     --outFileNamePrefix output_prefix_ \
     --outSAMtype BAM SortedByCoordinate \
     --runThreadN 8 
     
# For the mapping step, I will just show you the above command, you need to work
# out how the command works for your samples
# You can use STAR --help for help.
```

Parameters:

- `--genomeDir`: Path to the genome index directory
- `--readFilesIn`: Paired-end FASTQ files (R1 and R2)
- `--readFilesCommand zcat`: For compressed .gz files (omit if uncompressed)
- `--outFileNamePrefix`: Prefix for output files
- `--outSAMtype BAM SortedByCoordinate`: Output sorted BAM file
- `--runThreadN`: Number of threads

### Step 4: Output Files

STAR generates several output files:

- `Aligned.sortedByCoord.out.bam`: Sorted BAM file with aligned reads
- `Log.final.out`: Summary mapping statistics
- `Log.out`: Detailed log of the run
- `Log.progress.out`: Progress information
- `SJ.out.tab`: Splice junction information

### Step 5: Index BAM File (Optional)

For downstream analysis, you may need to index the BAM file:

```bash
samtools index output_prefix_Aligned.sortedByCoord.out.bam
```

## 2.3 QC on the alignments

To perform quality control on RNA-seq alignments using QualiMap, follow these steps:

### Step 1: Install QualiMap

Using conda or mamba to install QualiMap:

```bash
# conda
conda install -c bioconda qualimap
# mamba
mamba install -c bioconda qualimap
```

### Step 2: Run QualiMap on BAM File

Use QualiMap's RNA-seq mode to analyze your sorted BAM file:

```bash
qualimap rnaseq \
    -bam output_prefix_Aligned.sortedByCoord.out.bam \
    -gtf /path/to/annotations.gtf \
    -outdir qualimap_results \
    --java-mem-size=4G
```

Parameters:

- `-bam`: Input sorted BAM file from STAR alignment
- `-gtf`: Path to the gene annotation GTF file (same as used in STAR indexing)
- `-outdir`: Output directory for QualiMap results
- `--java-mem-size`: Amount of memory to allocate (adjust based on your system)

### Step 3: Additional Options

You can customize QualiMap with additional parameters:

```bash
qualimap rnaseq \
    -bam output_prefix_Aligned.sortedByCoord.out.bam \
    -gtf /path/to/annotations.gtf \
    -outdir qualimap_results \
    -pe \
    --java-mem-size=4G
```

- `-pe`: Specify paired-end mode
- `-p`: Sequencing protocol (strand-specific or non-strand-specific)

### Step 4: Review Quality Control Report

QualiMap generates an HTML report with comprehensive quality metrics including:

- Coverage profile along genes
- Junction analysis
- Reads genomic origin (exonic, intronic, intergenic)
- Coverage per chromosome
- Insert size distribution (for paired-end reads)
- Mapping quality distribution

Open the HTML report in the output directory (`qualimap_results/qualimapReport.html`) to review detailed quality metrics of your RNA-seq alignments.

Paper: https://academic.oup.com/bioinformatics/article/28/20/2678/206551

## 2.4 Quantification

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