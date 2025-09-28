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

### Usage

## Running the App
1. Clone or download this repository.
2. Open R (or RStudio) in the repository directory.
3. Run the script:

```
source("nitro_app_v0.03.5.R")
```

4. The Nitrogen Symbiosis UI will launch automatically in your browser.

*All dependencies are installed automatically on first run.*
