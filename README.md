# Fungal-ITS-identification-pipeline
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
  - [Sub-Workflow 3: Consensus-Based/convertfasta-based Identification](#sub-workflow-3-consensus-basedconvertfasta-based-identification)
- [Output Structure](#output-structure)
- [Interpretation of Results](#interpretation-of-results)
- [Acknowledgements](#acknowledgements)
- [License](#license)

---

## Overview
This pipeline processes Nanopore long-read data from PCR-amplified ITS regions using three distinct subworkflows:

- Direct ITS classification via GermGenie
- Consensus-based taxonomic identification (ITSx + BLAST)
- Convert fasta for taxonomic identification (ITSx + BLAST)
- Mapping and specificity assessment using Minimap2, Samtools + bedtools

Each step supports downstream analysis for high-quality Fungal isolates.

---

## Materials
The following software tools, platforms, and databases were used to build and run the fungal ITS identification pipeline:

### Input Data

This project analyzed sequencing data from 10 clinically relevant fungal isolates. These isolates originated from pure cultures obtained from microbiology laboratories.

- **Samples**: 10 fungal isolates from clinical origin
- **Target region**: Internal Transcribed Spacer (ITS1–5.8S–ITS2)
- **Primers used**:  
  - Forward: `ITS1-F_KYO2a`  
  - Reverse: `RCA95m`
- **Sequencing platform**: Oxford Nanopore MinION
- **Read type**: Single-end (FASTQ format)
- **Read characteristics**:
  - Mean read length: ~2,625 bp

#### Sample Overview

| FASTQ File                  | Barcode | Species                     |
|----------------------------|---------|-----------------------------|
| 2425-008_barcode01.fastq   | 01      | *Candida albicans*          |
| 2425-008_barcode02.fastq   | 02      | *Cryptococcus neoformans*   |
| 2425-008_barcode03.fastq   | 03      | *Lichtheimia ramosa*        |
| 2425-008_barcode04.fastq   | 04      | *Trichosporon asahii*       |
| 2425-008_barcode05.fastq   | 05      | *Candida tropicalis*        |
| 2425-008_barcode06.fastq   | 06      | *Trichophyton indotineae*   |
| 2425-008_barcode07.fastq   | 07      | *Aspergillus fumigatus*     |
| 2425-008_barcode08.fastq   | 08      | *Candida parapsilosis*      |
| 2425-008_barcode09.fastq   | 09      | *Aspergillus flavus*        |
| 2425-008_barcode10.fastq   | 10      | *Malassezia pachydermatis*  |


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
- `Bedtools` – used to:
  - convert `.bam` to `.bed` format with `bedtools bamtobed`  
  - merge overlapping mapped regions with `bedtools merge`  
  - identify and summarize genomic regions where reads aligned

### Reference Databases
- **UNITE fungal ITS database** – used for BLAST-based identification
  Downloaded from: UNITE
  File name: UNITE_public_19.02.2025.fasta.gz (131 MB)
  Download link: https://doi.plutof.ut.ee/doi/10.15156/BIO/3301227
  
- **Custom EMU ITS fungi database** – used for GermGenie classification
  Downloaded from: EMU
  File name: 	sh_general_release_10.05.2021.tgz (131 MB)
  Download link: https://doi.plutof.ut.ee/doi/10.15156/BIO/1280049
  
- **Fungal reference genomes** – for mapping
  NCBI reference sequences

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
conda create -n QC_env nanoplot filtlong bedtools -c bioconda -c conda-forge
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
conda create -n mapping_env minimap2 samtools blastn -c bioconda -c conda-forge
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

### BEDTOOLS: Convert BAM to BED and Merge Mapped Regions

After sorting the aligned reads (from `minimap2`) into `.sorted.bam` files, we used `bedtools` to extract the genomic coordinates of the aligned regions.

This two-step process was applied to all clinical isolate samples:

1. **Convert** `.bam` to `.bed` using `bedtools bamtobed`
2. **Merge** overlapping/adjacent regions using `bedtools merge`

This provides a compact overview of where reads are mapping in the reference genome 
### BEDTOOLS: Convert BAM to BED
Example: Candida albicans

```bash
mkdir -n bedtools_bamtobed

bedtools bamtobed -i minimap_samples/sorted_bam/barcode01_c_albicans_genomic.sorted.bam \
> bedtools_bamtobed/albicans_reads.bed
```
> Repeat this step for each barcode.

### BEDTOOLS: Merge Mapped Read Regions
Example: Candida albicans
```bash
mkdir -n bedtools_merge

bedtools merge -i bedtools_bamtobed/albicans_reads.bed \
> bedtools_merge/albicans_merged.bed
```
> Repeat this step for each barcode.

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

---

## Sub-Workflow 3: Consensus-Based/convertfasta-based Identification

This sub-workflow focuses on generating consensus sequences from Nanopore amplicons, using either *de novo* assembly or reference-guided variant calling. The goal is to reconstruct accurate ITS regions for downstream tools such as ITSx and BLAST, which require FASTA input.

It includes three main strategies:
- `Flye` assembly of filtered reads (de novo) for ITSx -> BLASTn
- `wf-amplicon` using both *de novo* and *variant calling* modes for ITSx -> BLASTn
- Conversion of filtered `.fastq` files to `.fasta` for ITSx -> BLASTn

All commands are executed from the working directory:  
`/mnt/studentfiles/2025/2025MBI06`

### Flye
Flye was used to perform de novo assembly of Nanopore amplicon reads for each clinical isolate (barcode01 to barcode10), without using a reference genome.

#### Flye installation
```bash
conda create -n flye_env -c bioconda -c conda-forge flye
conda activate flye_env
```

#### Run Flye for All Barcodes with filtered reads (70%)
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

### Wf-amplicon (de novo and reference-based)
wf-amplicon was used to reconstruct amplicons either:
- Without a reference (de novo mode)
- Using an ITS-region reference (variant calling mode)

#### wf-amplicon installation
```bash
conda create -n wfamplicon_env -c bioconda -c conda-forge nextflow
conda activate wfamplicon_env
```

#### De Novo Mode Example:
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

#### Install seqtk
```bash
conda create -n seqtk_env -c bioconda -c conda-forge seqtk
conda activate seqtk_env
```

#### Conversion Script
```bash
#!/bin/bash

mkdir -p ITSconvertfasta

for i in $(seq -w 1 10); do
  seqtk seq -A "filtlong_samples_70/barcode${i}_filtered.fastq" > \
  "ITSconvertfasta/barcode${i}_filtered.fasta"
done
```

### ITS Region Extraction with ITSx

ITSx was used to extract the Internal Transcribed Spacer (ITS) regions from consensus FASTA sequences. These FASTA files were obtained from four different workflows:
- Flye assemblies (assembly.fasta)
- wf-amplicon (reference-guided) output (reference_sanitized_seqIDs.fasta)
- wf-amplicon (de novo) output (consensus.fasta)
- seqtk converted filtered reads (barcode##_filtered.fasta)

Each dataset was placed in a separate input folder, and a corresponding output folder was created for the extracted ITS regions.

#### ITSx installation
```bash
conda create -n itsx_env -c bioconda -c conda-forge ITSx
conda activate istx_env
```

#### Example: ITSx on wf-amplicon (reference-based)

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
> To process the other datasets, only the INPUT_DIR, OUTPUT_DIR and fasta_file variable need to be adjusted accordingly.
> Repeat this process for the other consensus sources by adjusting the INPUT_DIR:
> - ITSflye -> ITSflye_output
> - wfamplicon_denovo -> ITSwfamplicon_denovo

### Run BLASTn on Extracted ITS Regions

After extracting ITS regions using ITSx, we performed BLASTn searches to identify fungal species.

We used the UNITE + INSD fungal ITS database, available at:
https://doi.plutof.ut.ee/doi/10.15156/BIO/3301227

The BLASTn search was applied to consensus ITS1 and ITS2 sequences extracted from the following sources:
- ITSflye (from Flye-based assemblies)
- ITSwfampliconref (from reference-guided wf-amplicon)
- ITSwfamplicon_de_novo (from wf-amplicon de novo assemblies)
- ITSconvertfasta (filtered reads converted to FASTA)

Each directory contains FASTA files named sample.ITS1.fasta and sample.ITS2.fasta, which were used as input for the BLASTn searches.

#### blastn installation
```bash
conda create -n blast_env -c bioconda -c conda-forge blastn
conda activate blast_env
```

#### Running blastn on the sources
```bash
#!/bin/bash
ITS_DIR="/mnt/StudentFiles/2025/2025MBI06/ITSwfamplicon"
DB="/mnt/StudentFiles/2025/2025MBI06/ITS_database/unite_ITS_db"
# Loop over each barcode folder and BLAST ITS1 and ITS2 sequences
for folder in "$ITS_DIR"/barcode*/; do
    sample=$(basename "$folder")

    ITS1="$folder/${sample}.ITS1.fasta"
    ITS2="$folder/${sample}.ITS2.fasta"

    # BLAST ITS1 region
    blastn -query "$ITS1" -db "$DB" -out "$folder/${sample}.ITS1.blast.txt" \
      -evalue 1e-5 -outfmt 6 -num_threads 4

    # BLAST ITS2 region
    blastn -query "$ITS2" -db "$DB" -out "$folder/${sample}.ITS2.blast.txt" \
      -evalue 1e-5 -outfmt 6 -num_threads 4
done
```

--- 

## Output Structure
Each sub-workflow produces its own organized set of results, stored in dedicated directories per sample (barcode01 to barcode10). Below is an overview of the key output folders:

### Preprocessing

- `nanoplot_samples/`  
  Contains NanoPlot outputs per sample. For each barcode, this includes a `NanoStats.txt` file with:
  - Mean read length and quality
  - Read length N50
  - Quality distribution
  - Top longest and highest quality reads

- `filtlong_samples_70/`  
  Filtered FASTQ reads retaining the best 70% per sample (used in all downstream workflows).
  

### Sub-Workflow 1: Mapping & Specificity

- `minimap_samples/`  
  Contains `.sam` files generated by `minimap2` aligning each barcode to its matching reference genome.

- `minimap_samples/sorted_bam/`  
  Contains sorted and indexed `.bam` files created with `samtools`. These are used to calculate mapping percentages and specificity.

- `a_aprimerblastnew/`  
  Contains `.txt` files with ITS amplicon coordinates derived from BLAST of primers. These regions are used to extract region-specific mappings with `samtools view`.


### Sub-Workflow 2: GermGenie

- `emu_results/`  
  Taxonomic classification results from `EMU`, using the prebuilt UNITE ITS fungal database.  
  Includes `abundance.tsv` per barcode with classified species and their relative abundance.


### Sub-Workflow 3: Consensus-Based/convertfasta-Based Identification

- `ITSflye/`  
  De novo consensus assemblies per barcode generated with `Flye`, used as input for ITSx.

- `wfampliconref_new/`  
  Consensus sequences generated by `wf-amplicon` using reference-based variant calling, used as input for ITSx.

- `wfamplicon_de_novo/`  
  De novo consensus sequences generated by `wf-amplicon` without a reference, used as input for ITSx.

- `ITSconvertfasta/`  
  FASTA-converted filtered reads using `seqtk`, used as input for ITSx.

- `ITSwfamplicon/`, `ITSflye_output/`, `ITSde_novo_output/`, `ITSconvert_output/`  
  ITSx results for ITS1 and ITS2 regions per input type. Each folder contains:
  - `*.ITS1.fasta`  
  - `*.ITS2.fasta`

- `*.ITS1.blast.txt`, `*.ITS2.blast.txt`  
  BLASTn results from querying extracted ITS regions (from ITSx) against the UNITE+INSD fungal ITS database.


--- 

## Interpretation of Results

The output files of each sub-workflow provide distinct types of biological and technical insight into the fungal ITS sequencing analysis. Below is a detailed explanation of how to interpret the most important result types.

### Quality Control (NanoPlot)

Each sample includes a `NanoStats.txt` file with summary statistics. Example fields:

- **Mean read length**: Average size of reads (e.g. `3,942.8` bp).
- **Mean read quality**: Average Phred quality score (e.g. `19.3`). Scores >10 indicate high basecall confidence.
- **N50**: Length at which 50% of total bases are in reads of this size or longer.
- **Quality cutoffs (Q15, Q20, Q25, etc.)**: Percentages of reads exceeding each quality threshold.

**Interpretation:**  
High-quality samples typically show:
- Mean read length >2000 bp
- Q ≥ 15
Use this to confirm if samples are of sufficient quality for downstream ITS analysis.


### Taxonomic Classification (EMU via GermGenie)
The `emu` output includes a tab-separated table with taxonomic ranks and **relative abundance per barcode**:

| tax_id | abundance | genus      | species          |
|--------|-----------|------------|------------------|
| 53     | 0.7783    | Candida    | Candida albicans |
| 3127   | 0.1188    | Candida    | -                |

**Interpretation:**
- The **dominant taxon** (abundance >0.5)
- Secondary taxa may indicate low-level contamination, database overlap, or off-target amplification.
- The output validates whether the intended fungal species was successfully amplified and sequenced.

**Example Interpretation**:
If `Candida albicans` appears with an abundance of `0.77`, this indicates dominant presence. Secondary hits (e.g. <0.05) may reflect aspecific amplification.

EMU uses a prebuilt UNITE-based fungal ITS database, so accuracy depends on database coverage.



### ITS Region Extraction (ITSx)

ITSx extracts the ITS1 and ITS2 regions from consensus FASTA sequences.

For each sample, the tool outputs:

- `*.ITS1.fasta` – Extracted ITS1 region
- `*.ITS2.fasta` – Extracted ITS2 region
- `.summary.txt`, `.log` – Summary and processing logs

Presence of both ITS1 and ITS2 suggests complete and usable ITS region extraction for BLAST identification.


### BLAST Results (ITS1 & ITS2)

BLAST results (`*.blast.txt`) are formatted in tabular form (outfmt 6) and include:

| Column | Description                |
|--------|----------------------------|
| 1      | Query sequence ID          |
| 2      | Matched subject ID         |
| 3      | % Identity                 |
| 4      | Alignment length           |
| 11     | E-value                    |
| 12     | Bit score                  |

**Interpretation Guidelines**:

- **% Identity** ≥ 97–100% indicates strong species-level match.
- **Low E-value** (e.g. `7.22e-42`) confirms statistical significance.
- **Taxon names** in hit descriptions guide fungal species assignment (e.g. *Candida parapsilosis*).



### Mapping Specificity (Minimap2 + Samtools + bedtools)

Mapping specificity is calculated as:
( Number of reads mapped to ITS region / Total mapped reads ) × 100

This is done with:
```bash
# Reads mapped to ITS region (coordinates in BED)
samtools view -b -F 4 -L region.bed sorted.bam | samtools fastq | awk 'END {print NR/4}'

# Total mapped reads
samtools view -b -F 4 sorted.bam | samtools fastq | awk 'END {print NR/4}'
```

A high percentage (e.g. >80%) suggests accurate and specific amplification of the fungal ITS region.
Low values may indicate poor primer performance, low-quality reads, or wrong reference genomes.

> The .sorted.bam files can be visually inspected using IGV (Integrative Genomics Viewer). Load both the sorted BAM and its index .bai file alongside the reference genome to verify whether reads cluster within  the expected ITS coordinates.

BEDTools analysis for mapping region summary:

We used bedtools bamtobed to convert BAM to BED, and then bedtools merge to summarize all mapped regions:
Interpretation:
- The merged BED file reveals the total extent of coverage.
- This can be compared to ITS region coordinates to assess overlap and detect non-specific mapping.

---

## Acknowledgements

This project makes use of several open-source tools and community resources. We gratefully acknowledge the developers and contributors of the following tools and databases:

- [`Filtlong`](https://github.com/rrwick/Filtlong) – for quality filtering of Nanopore reads  
- [`NanoPlot`](https://github.com/wdecoster/NanoPlot) – for quality control and read length statistics  
- [`EMU`](https://github.com/treangenlab/emu) – for ITS-based taxonomic classification  
- [`Flye`](https://github.com/fenderglass/Flye) – for de novo genome assembly from long reads  
- [`wf-amplicon`](https://github.com/epi2me-labs/wf-amplicon) – for amplicon consensus generation and variant calling  
- [`seqtk`](https://github.com/lh3/seqtk) – for conversion between FASTQ and FASTA formats  
- [`ITSx`](https://anaconda.org/bioconda/itsx) – for ITS1 and ITS2 region extraction  
- [`BLASTn`](https://github.com/JacobLondon/Blastn) – for sequence alignment and species identification  
- [`minimap2`](https://github.com/lh3/minimap2) – for mapping reads to fungal reference genomes  
- [`samtools`](https://github.com/samtools/samtools) – for BAM processing and alignment statistics  
- [`BEDTools`](https://github.com/arq5x/bedtools2) – for converting and merging to get the locations
- [`gunzip`](https://www.gnu.org/software/gzip/) – for decompressing .gz archive files  
- [`IGV (Integrative Genomics Viewer)`](https://github.com/igvteam/igv) – for visualizing read alignments and genome features  

**Additional resources:**
- [UNITE Database](https://unite.ut.ee/) – for curated fungal ITS reference sequences  
- [NCBI RefSeq](https://www.ncbi.nlm.nih.gov/refseq/) – for fungal genome assemblies and annotations

---

## License

This project is licensed under the [GNU General Public License v3.0 (GPLv3)](https://www.gnu.org/licenses/gpl-3.0.en.html).
