## RNA-seq Mapping

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