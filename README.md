# fungal-ITS-identification-pipeline
This repository contains a bioinformatics pipeline for fungal species identification based on Nanopore sequencing of PCR-amplified ITS regions from clinical fungal isolates.
![image](https://github.com/user-attachments/assets/aad0ee9d-2e40-4e7a-9578-29545650b256)

## Goal
This study was initiated in order to investigate whether Nanopore sequencing of PCR-amplified ITS regions, using the primers ITS1-F_KYO2a and RCA95m, permits accurate fungal species identification from clinical fungal isolates. This approach serves as a model for evaluating the potential of ITS-based diagnostics in settings with degraded DNA, such as FFPE tissues, even though the current study focuses on DNA extracted from high-quality clinical fungal isolates.

---

## Table of Contents
- [Overview](#overview)
- [Materials](#materials)
- [Implementation](#implementation)
  - [Sub-Workflow 1: GermGenie](#sub-workflow-1-germgenie)
  - [Sub-Workflow 2: Consensus-Based Identification](#sub-workflow-2-consensus-based-identification)
  - [Sub-Workflow 3: Mapping & Specificity](#sub-workflow-3-mapping--specificity)
- [Output Structure](#output-structure)
- [Interpretation of Results](#interpretation-of-results)
- [Acknowledgements](#acknowledgements)
- [License](#license)

---

## Overview
This pipeline processes Nanopore long-read data from PCR-amplified ITS regions using three distinct subworkflows:

- Direct ITS classification via GermGenie
- Consensus-based taxonomic identification (ITSx + BLAST)
- Mapping and specificity assessment using Minimap2 and Samtools

Each step supports downstream analysis for high-quality Fungal isolates.

---

## Materials
The following software tools, platforms, and databases were used to build and run the fungal ITS identification pipeline:

### Software Tools
- `Filtlong` – for quality filtering of Nanopore reads  
- `NanoPlot` – for sequencing QC and length distribution plots  
- `EMU` (via GermGenie) – for direct ITS-based taxonomic classification  
- `Flye` – for de novo consensus assembly  
- `wf-amplicon` – for variant-based and de novo consensus workflows  
- `seqtk` – for FASTQ to FASTA conversion  
- `ITSx` – for extraction of ITS1 and ITS2 regions from consensus FASTA  
- `BLAST+` – for fungal species identification against UNITE  
- `Minimap2` – for aligning reads to fungal reference genomes  
- `Samtools` – for BAM processing and calculation of mapped read counts via `samtools view`  
- `Bedtools intersect` – for determining overlap of mapped reads with annotated ITS regions (PCR specificity)

### Reference Databases
- **UNITE fungal ITS database** – used for BLAST-based identification  
- **Custom EMU ITS fungi database** – used for GermGenie classification  
- **Fungal reference genomes + GFF annotations** – for mapping and intersect analyses

### Platforms and Workflow Managers
- **Oxford Nanopore MinION** – sequencing platform used to generate PCR amplicon reads  
- **Asgard** – HPC/cloud platform used to run workflows and heavy analyses  
- **Conda** – for environment and package management  
- **Nextflow** – to execute the `wf-amplicon` pipeline

