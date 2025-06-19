# fungal-ITS-identification-pipeline
This repository contains a bioinformatics pipeline for fungal species identification based on Nanopore sequencing of PCR-amplified ITS regions from clinical fungal isolates.
![image](https://github.com/user-attachments/assets/d78acfce-ede3-4ec1-9483-106b6c8d2090)

## Goal
This study was initiated in order to investigate whether Nanopore sequencing of PCR-amplified ITS regions, using the primers ITS1-F_KYO2a and RCA95m, permits accurate fungal species identification from clinical fungal isolates. This approach serves as a model for evaluating the potential of ITS-based diagnostics in settings with degraded DNA, such as FFPE tissues, even though the current study focuses on DNA extracted from high-quality clinical fungal isolates.

---

## Table of Contents
- [Overview](#overview)
- [Materials](#materials)
- [Implementation](#implementation)
  - [Preprocessing: Read Filtering & Quality Control](#preprocessing-read-filtering--quality-control)
  - [Sub-Workflow 1: Mapping & Specificity](#sub-workflow-1-mapping--specificity)
  - [Sub-Workflow 2: GermGenie](#sub-workflow-2-germgenie)
  - [Sub-Workflow 3: Consensus-Based/convertfasta-based Identification](#sub-workflow-3-consensus-based-identification)
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
- `BLAST+` – for fungal species identification against UNITE and primer blast against reference genome. 
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

---

## Preprocessing: Read Filtering & Quality Control

Before running any identification workflow, raw Nanopore reads are filtered and quality-checked using `NanoPlot` and `Filtlong`.  
All project data is stored in the `project_data/` directory, containing raw FASTQ files from 10 PCR-amplified clinical fungal isolates.

All commands are executed from the following working directory: /mnt/studentfiles/2025/2025MBI06.


### Installation

Create a dedicated conda environment for QC tools:
```bash
conda create -n QC_env nanoplot filtlong -c bioconda -c conda-forge
conda activate QC_env
```

### Visualize Raw Reads with NanoPlot

```bash
#!/bin/bash

# Create an output directory for NanoPlot results
mkdir -p nanoplot_samples

# Loop through all FASTQ files in the input folder
for file in project_data/*.fastq; do
    # Get the filename without the extension
    name=$(basename "$file" .fastq)

    # Create a subfolder for each sample
    mkdir -p "nanoplot_samples/$name"

    # Run NanoPlot and output results into the subfolder
    NanoPlot --fastq "$file" --outdir "nanoplot_samples/$name" --loglength --N50
done
```

### Filter Reads with Filtlong

To remove low-quality and short Nanopore reads, we used `Filtlong` to retain only high-quality sequences.  
Each `.fastq` file (one per clinical isolate) was filtered with a minimum read length of 1000 bp, keeping only the top 70% of reads.

Filtered files were saved to the directory: `filtlong_samples_70/`.

```bash
filtlong --min_length 1000 --keep_percent 70 project_data/2425-008_barcodeXX.fastq > filtlong_samples_70/barcodeXX_filtered.fastq
```
Replace XX with the barcode number (01 to 10).
For example:
```bash
filtlong --min_length 1000 --keep_percent 70 project_data/2425-008_barcode01.fastq > filtlong_samples_70/barcode01_filtered.fastq
```
> **Note:**  
> We tested filtering at different thresholds (`--keep_percent 70`, `60`, and `50`).  
> Based on read quality and yield, we continued all downstream analyses using the 70% filtered datasets.

### Visualize Filtered Reads with NanoPlot

After filtering, each `.fastq` file is visualized using `NanoPlot`.  
Results are saved in a separate folder for each barcode inside the same output directory.

```bash
#!/bin/bash

# Set the folder where the filtered FASTQ files are stored
map="filtlong_samples_70"

# Go through all filtered FASTQ files
for file in "$map"/barcode*_filtered.fastq; do
    # Get the barcode name (like barcode01)
    barcode=$(basename "$file" | cut -d'_' -f1)

    # Create a folder to save the plots
    mkdir -p "$map/$barcode"

    # Run NanoPlot on the file and save output
    NanoPlot --fastq "$file" --outdir "$map/$barcode" --plots hex dot --loglength --N50
done
```

--- 

## Sub-Workflow 1: Mapping & Specificity

This workflow evaluates how specifically each sample maps to the correct fungal species and to the annotated ITS region of that species.

It uses:
- `minimap2` to align reads to fungal reference genomes
- `samtools` to convert and process the alignment files
- `samtools view` + `blasted primers` to calculate the percentage % specific mapped

All commands are executed from the following working directory: /mnt/studentfiles/2025/2025MBI06.


### Reference Genomes

Fungal reference genome files were downloaded from NCBI and pre-unzipped before running this step. Each clinical isolate (barcode01–barcode10) was assigned a matching `.fna` reference genome, stored in the `All_fna/` folder.

Examples:
- `barcode01 -> c_albicans_genomic.fna`
- `barcode02 -> c_neofromans_genomic.fna`
- ...
- `barcode10 -> m_pachydermatis_genomic.fna`

### Installation

Install all required tools in a dedicated Conda environment:

```bash
conda create -n mapping_env minimap2 samtools bedtools blastn -c bioconda -c conda-forge
conda activate mapping_env
```

### Mapping Reads to Reference Genomes

Create a folder for the mapping results:

```bash
mkdir -p Minimap_samples
```
Then map each sample to its corresponding reference genome. Example for barcode01:
```bash
minimap2 -ax map-ont All_fna/c_albicans_genomic.fna project_data/2425-008_barcode01.fastq > Minimap_samples/barcode01_c_albicans.sam
```

### Convert and Sort Alignment Files with Samtools

Each `.sam` file from the `minimap_samples/` folder was converted to `.bam`, sorted, and indexed.  
The sorted files were stored in a subfolder called `sorted_bam/`, located inside `minimap_samples/`.

```bash
# create output folder
mkdir -p minimap_samples/sorted_bam

# Convert SAM to BAM
samtools view -bS minimap_samples/barcode01_c_albicans_genomic.sam > minimap_samples/barcode01_c_albicans_genomic.bam

# Sort the BAM file
samtools sort minimap_samples/barcode01_c_albicans_genomic.bam -o minimap_samples/sorted_bam/barcode01.sorted.bam

# Index the sorted BAM
samtools index minimap_samples/sorted_bam/barcode01.sorted.bam
```

Repeat these three commands for each barcode (barcode01 to barcode10), updating the filenames accordingly.


### Primer-to-Genome BLAST Search

To verify where primers bind in the fungal genomes, each primer was aligned against each `.fna` reference genome using `blastn`.
This step helps determine expected amplicon positions for alignment validation and extraction.

Prepare Primer FASTA File

Each primer (forward or reverse) must be saved as a `.fasta` file with a meaningful name. Example:
```fasta
>Fw1
TAGAGGAAGTAAAAGTCGTAA
>Rv5
...
```

### Run PRIMER BLAST Against All Genomes
```bash
#!/bin/bash

# Path to the primer and genome directory
PRIMER="/mnt/StudentFiles/2025/2025MBI06/all_fna/primerF1.fasta"
DB_DIR="/mnt/StudentFiles/2025/2025MBI06/all_fna"
OUT_DIR="./primerblastnew"

mkdir -p "$OUT_DIR"

# Loop through each genome file and run BLAST
for DB in "$DB_DIR"/*.fna; do
    BASENAME=$(basename "$DB" .fna)

    # Create BLAST database
    makeblastdb -in "$DB" -dbtype nucl

    # Run blastn with word_size 22 for short primer hits
    blastn -query "$PRIMER" -db "$DB" -word_size 22 -out "$OUT_DIR/${BASENAME}_blast.txt" -outfmt 6
done
```

### Calculate Mapping Specificity Percentages

To evaluate how specifically reads map to their expected fungal ITS regions, we calculated the percentage of mapped reads that fall within the primer-defined amplicon coordinates.

This was done using samtools view in combination with primer coordinates (from BLAST) stored in .bed files. Read counts were obtained by converting BAM output to FASTQ and dividing the number of lines by 4.

The percentage is calculated using the formula:
( number of reads in ITS region / total number of mapped reads ) × 100.

```bash
# Example for Barcode 01 (Candida albicans)
# Get the number of reads mapped to the expected ITS region
samtools view -b -F 4 minimap_samples/sorted_bam/barcode01_c_albicans_genomic.sorted.bam NC_032096.1:1892965-1895597 \
| samtools fastq | awk 'END {print NR/4}'

# Get the total number of mapped reads
samtools view -b -F 4 minimap_samples/sorted_bam/barcode01_c_albicans_genomic.sorted.bam \
| samtools fastq | awk 'END {print NR/4}'
```

--- 

## Sub-Workflow 2: GermGenie

This workflow performs direct taxonomic classification of fungal ITS reads using [GermGenie](https://github.com/czbiohub/GermGenie), which internally uses the EMU classifier.

Before using EMU, download the prebuilt fungal ITS database provided by the UNITE Community:

**Database source:**  
https://doi.plutof.ut.ee/doi/10.15156/BIO/1280049

### Download and Extract the Database

```bash
# Extract the prebuilt UNITE database (downloaded as .tgz and .tar)
tar -xvzf sh_general_release.tgz
tar -xvzf sh_general_release.tar
```

### Create a dedicated Conda environment for EMU and activate the environment
```bash
conda create -n emu_env -c bioconda -c conda-forge emu
conda activate emu_env
```

### Configure EMU Environment Variables
```bash
# Set EMU database directory and database name
export EMU_DATABASE_DIR=~/emu_dbs
export EMU_PREBUILT_DB='unite-fungi'

# Create and move into the database directory
mkdir -p ${EMU_DATABASE_DIR}
cd ${EMU_DATABASE_DIR}
```



### Creating new folder and run EMU on a Sample
```bash
# Create new foler
mkdir -p emu_results
# Example: classify reads from barcode01
emu abundance /filtlong_samples_70/barcode01_filtered.fastq \
  --db $EMU_DATABASE_DIR \
  --output-dir ./emu_results/barcode01
```
Repeat this step for each filtered sample.

## Sub-Workflow 3: Consensus-Based Identification

This sub-workflow focuses on generating consensus sequences from the Nanopore amplicons, either via *de novo* assembly or reference-guided variant calling. The goal is to reconstruct the fungal ITS regions for more accurate identification.

It includes three main strategies:
- `Flye` assembly of filtered reads (de novo)
- `wf-amplicon` using both *de novo* and *variant calling* modes
- Conversion of filtered `.fastq` files to `.fasta` for ITSx

All commands are executed from the working directory:  
`/mnt/studentfiles/2025/2025MBI06`

### Flye installation
Install Flye

```bash
conda create -n flye_env -c bioconda -c conda-forge flye
conda activate flye_env
```

### Run Flye for All Barcodes with filtered reads (70%)
Assemble filtered reads (70%) from Filtlong

```bash
#!/bin/bash

for i in $(seq -w 1 10); do
    flye \
      --nano-raw "filtlong_samples_70/barcode${i}_filtered.fastq" \
      --out-dir "./ITSflye/barcode${i}" \
      --threads 8
done
```

### wf-amplicon (de novo and reference-based)
wf-amplicon was used to reconstruct amplicons either:
- Without a reference (de novo mode)
- Using an ITS-region reference (variant calling mode)

De Novo Mode Example:
```bash
nextflow run epi2me-labs/wf-amplicon \
  --fastq project_data/2425-008_barcode01.fastq \
  --sample barcode01 \
  --out_dir wfamplicon_denovo/barcode01 \
  -profile standard
```
Variant Calling Mode (Requires ITS Reference) example:
Reference FASTA files were extracted from the expected ITS regions of each fungal genome using BLAST coordinates. These reference files were stored in the chromosome/ directory.
Variant calling example:
```bash
nextflow run epi2me-labs/wf-amplicon \
  --fastq project_data/2425-008_barcode01.fastq \
  --sample barcode01 \
  --reference chromosome/albicans.fasta \
  --out_dir wfampliconref_new/barcode01 \
  -profile standard
```

### Convert FASTQ to FASTA
Before running ITSx, the filtered .fastq reads of the clinical isolates must be converted to .fasta, since ITSx only supports FASTA input.

### Install seqtk
```bash
conda create -n seqtk_env -c bioconda -c conda-forge seqtk
conda activate seqtk_env
```

### Conversion Script
```bash
#!/bin/bash

mkdir -p ITSconvertfasta

for i in $(seq -w 1 10); do
  seqtk seq -A "filtlong_samples_70/barcode${i}_filtered.fastq" > \
  "ITSconvertfasta/barcode${i}_filtered.fasta"
done
```

### 4. ITS Region Extraction with ITSx

ITSx was used to extract the Internal Transcribed Spacer (ITS) regions from consensus sequences obtained via:

- `Flye` (de novo)
- `wf-amplicon` (reference-based)
- `wf-amplicon` (de novo)

Each set of consensus FASTA sequences was processed separately. For each workflow, a new output directory was created.

Example: ITSx on wf-amplicon (reference-based)

```bash
#!/bin/bash

INPUT_DIR="/mnt/StudentFiles/2025/2025MBI06/wfampliconref_new"
OUTPUT_DIR="/mnt/StudentFiles/2025/2025MBI06/ITSwfamplicon"

mkdir -p "$OUTPUT_DIR"

for folder in "$INPUT_DIR"/barcode*/; do
  sample=$(basename "$folder")
  fasta_file="$folder/reference_sanitized_seqIDs.fasta"

  mkdir -p "$OUTPUT_DIR/$sample"

  ITSx \
    -i "$fasta_file" \
    -o "$OUTPUT_DIR/$sample/$sample" \
    -t F \
    --cpu 4
done
```

> Output: One directory per sample containing ITSx results (e.g., ITSx output, .ITS1.fasta, .ITS2.fasta, .log.txt, etc.)
