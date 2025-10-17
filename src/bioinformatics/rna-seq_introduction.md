# Introduction — RNA-seq analysis

RNA sequencing (RNA‑seq) measures transcriptome-wide RNA abundance and sequence, enabling discovery and quantification of expressed genes, isoforms, and sequence variants. This chapter provides a concise practical introduction to the computational steps, core concepts, and best practices used in typical RNA‑seq studies.

## Goals
- Quantify gene and transcript expression across samples
- Identify differentially expressed genes and pathways
- Detect alternative splicing, fusion transcripts, and RNA editing when relevant
- Produce reproducible, well-documented results ready for interpretation

## Key concepts
- Reads (FASTQ), alignments (BAM), and quantified counts/tpm
- Library type: strandedness, read length, single/paired-end
- Normalization to account for sequencing depth and composition
- Experimental design and covariates (batch, replication)

## Typical computational workflow
1. Project setup: metadata, sample sheet, versioned code and environment (conda/Docker/Singularity)
2. Raw-data QC: fastp, MultiQC
3. Alignment: STAR
4. Quantification: featureCounts, S
5. Normalization and exploratory QC: PCA, sample clustering, library size checks
6. Differential expression: DESeq2, edgeR, limma-voom
7. Functional analysis: GSEA, pathway enrichment, gene-set visualization
8. Optional: isoform analysis, fusion detection, variant calling
9. Reporting: figures, tables, reproducible notebooks/workflows

## Outputs and quality checks
- Per-sample QC reports and alignment metrics
- Count matrix and normalized expression table
- Differential expression results with effect sizes and adjusted p-values
- Diagnostic plots: MA, PCA, heatmaps, volcano plots

## Best practices
- Record sample metadata and processing parameters early
- Use workflow managers (Snakemake, Nextflow) and containerized environments
- Include replication and model covariates in statistical tests
- Validate key findings with orthogonal methods if possible

This chapter will expand each step with commands, recommended tools, example workflows, and interpretation guidance.