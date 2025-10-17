## QC on the alignments

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