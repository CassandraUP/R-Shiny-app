## Nitrogen Symbiosis RNA-seq Heatmap App
This repository provides an R-based interactive application for initially designed for **exploring nitrogen-fixing symbiosis gene expression**. The app was developed for ***Encephalartos natalensis*** and related cyanobacterial host plants, with comparative insights from legumes and other well-studied nitrogen-fixing symbiotic systems. This app can be modified to also support broader comparative transcriptomics analyses across plant species.

### About
The app allows researchers to upload RNA-seq expression datasets and interactively generate **publication-quality plots and heatmaps**. It was initially designed to help uncover patterns in nitrogen symbiosis, focusing on:
- ***E. natalensis*** (our focal species)
- **Cyanobacteria-hosting plants** (four representative hosts)
- **Legumes and well-studied symbiotic plants** for comparative analysis

**Key Features**
- **Interactive User Interface** – built with Shiny; no coding required.
- **Single-gene exploration** – plot expression profiles across tissues for any *E. natalensis* gene ID.
- **Heatmap visualization** – cluster multiple genes to compare tissue-specific expression patterns.
- **Orthogroup integration** – link genes with orthologous groups for comparative studies.
- **Customizable outputs** – adjust clustering, scaling, and figure dimensions.
- **Export options** – download plots (PNG/PDF) and result tables for downstream analysis or publication.

This app bridges exploratory data analysis with publication-ready visualization, making it a useful resource for plant symbiosis research and beyond.

### Required Input Files
The Nitrogen Symbiosis RNA-seq Heatmap App requires the following input files to run properly:

#### 1. Expression Matrix
- **Format:** Tab-delimited .txt or .tsv file.
- **Structure:**
  - Rows = Gene IDs.
  - Columns = samples (e.g., tissues, experimental conditions).

```
GeneIDs  +N_R1    -N_S3    +N_L2
Gene1    12.4     5.2      3.1
Gene2    2.1      0.0      15.2
Gene3    8.3      4.5      3.4
```

#### 2. Sample Metadata (Optional, but recommended)
*Links each column in the expression matrix to biological metadata (e.g., tissue type, treatment, species).*
- **Format:** Tab-delimited .txt or .tsv.
- **Structure:**
  - First column = Sample name (must match headers in expression matrix).
  - Additional columns = Experimental details (e.g., tissue, condition, replicate).

```
SampleID  Tissue  Condition  Species
+N_R1     Root    Control    E. natalensis
-N_S3     Stem    Deficient  E. natalensis
+N_L2     Leaf    Control    E. natalensis
```

#### 3. Orthogroup Mapping File (Optional)
*Allows cross-species comparison by mapping genes to orthogroups.*
- **Format:** Tab-delimited .txt or .tsv.
- **Structure:**
  - Columns = Orthogroup ID, Gene ID, Species.

```
Orthogroup   GeneID      Species
OG0001       Enat001     E.natalensis
OG0012       Glyma.123   G. max
OG0001       Medtr.456   M. truncatula
```

### Usage

#### Running the App
1. Clone or download this repository.
2. Open R (or RStudio) in the repository directory.
3. Run the script:

```
source("nitro_app_v0.03.5.R")
```

4. The Nitrogen Symbiosis UI will launch automatically in your browser.

*All dependencies are installed automatically on first run.*
