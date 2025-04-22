# Genomics Workflow: Read Processing, Assembly, and Comparative Analysis

This project outlines a complete genomics workflow to process raw sequencing data, assemble genomes, and perform comparative genomics analysis. The pipeline involves the use of various command-line tools and bioinformatics packages for quality control, de novo assembly, species identification, and genome completeness assessment.

## 🧬 Overview

The goal of this workflow is to:
- Download and clean raw sequencing reads from the NCBI SRA
- Assemble high-quality genomes using SPAdes
- Filter contigs and compare assemblies with a known reference using FastANI
- Identify species and sequence types using MLST
- Assess genome completeness and contamination with BUSCO

This workflow was executed on a 2-core machine using Conda-managed environments for reproducibility and dependency control.

## 🔧 Steps Included

### 1. **Environment Setup**
- A `genomics3` Conda environment was created with all necessary tools for data preprocessing and analysis.

### 2. **Read Download and Preprocessing**
- SRA reads (SRR27160578–80) were downloaded using `prefetch` and converted to FASTQ format with `fastq-dump`.
- FASTQ files were compressed using `pigz` and quality-checked/cleaned using `fastp`.

### 3. **Genome Assembly**
- Reads were assembled using **SPAdes**, and resulting `contigs.fasta` files were renamed and filtered to retain contigs >500 bp using `seqtk`.

### 4. **Comparative Genomics with FastANI**
- Each filtered assembly was compared to a reference genome (Bordetella parapertussis) using **FastANI**.
- Output files were processed to compute percent alignment and aligned basepairs, and results were combined into a summary table.

### 5. **Species Identification with MLST**
- MLST was run on the reference and assemblies to determine species and sequence types.
- Results were formatted into a clean TSV with

### 6. **Genome Quality Assessment with BUSCO**
- BUSCO was run using the `bacteria_odb12` lineage dataset to estimate genome completeness and contamination.
- Results were parsed and summarized in a TSV file.

## 📂 Output Files

- `fastani.tsv`: Pairwise ANI comparisons between assemblies and the reference.
- `mlst.tsv`: Sequence type information for each assembly.
- `quality.tsv`: BUSCO-based genome completeness and contamination metrics.

## 💻 Dependencies and Tools

| Tool        | Use Case                          |
|-------------|-----------------------------------|
| `sra-tools` | Downloading sequencing reads      |
| `fastp`     | Read quality control              |
| `SPAdes`    | De novo genome assembly           |
| `seqtk`     | Contig filtering                  |
| `FastANI`   | Genome similarity comparison      |
| `mlst`      | Species/strain identification     |
| `BUSCO`     | Genome completeness estimation    |

All tools were installed via **bioconda** or **conda-forge** channels using Conda.

## 📝 Notes
- Execution time for SPAdes assembly is longer on low-core CPUs (~20 mins per sample).
- BUSCO results were manually interpreted due to the lack of explicit contamination metrics.

