# Differential Expression Analysis of Metabolic Genes in ER+ Breast Cancer

## Overview

This repository contains R code for paired differential expression analysis of ketone body and one-carbon metabolism genes in ER+ breast cancer tumors compared to matched adjacent normal tissue.

**Dataset:** GSE58135 (NCBI Gene Expression Omnibus)  
**Analysis:** DESeq2 paired design with strict 1:1 tumor-normal matching  
**Focus Genes:** 9 genes involved in ketone body metabolism (BDH1, BDH2, OXCT1, OXCT2, ACAT1, ACAT2) and one-carbon metabolism (SHMT1, SHMT2, BCAT1)

## Key Findings

- **OXCT1** (log2FC = -1.41, adjusted p = 5.75×10⁻⁷) - Significantly downregulated in tumors
- **SHMT1** (log2FC = -1.31, adjusted p = 1.62×10⁻⁵) - Significantly downregulated in tumors  
- **ACAT1** (log2FC = -1.07, adjusted p = 1.07×10⁻⁴) - Significantly downregulated in tumors

These results suggest selective suppression of ketone body catabolism and cytoplasmic one-carbon metabolism in ER+ breast cancer.

## Requirements

### Software
- **R version:** 4.5.1 (or higher)
- **RStudio:** Recommended but not required

### R Packages

#### CRAN packages:
- `data.table` (1.17.8)
- `ggplot2` (4.0.0)
- `pheatmap` (1.0.12)
- `reshape2` (1.4.4)

#### Bioconductor packages:
- `GEOquery` (2.76.0)
- `DESeq2` (1.48.1)
- `org.Hs.eg.db` (3.21.0)

## Installation

### 1. Install R and RStudio
Download and install R from [CRAN](https://cran.r-project.org/) and RStudio from [Posit](https://posit.co/download/rstudio-desktop/).

### 2. Install Required Packages
The script will automatically install required packages. Alternatively, run:

```r
# Install BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE)) 
    install.packages("BiocManager")

# Install Bioconductor packages
BiocManager::install(c("GEOquery", "DESeq2", "org.Hs.eg.db"))

# Install CRAN packages
install.packages(c("data.table", "ggplot2", "pheatmap", "reshape2"))
```

## Usage

### Quick Start

1. Clone this repository:
```bash
git clone https://github.com/mweichhaus/Metabolomics-and-Ketones.git
cd Metabolomics-and-Ketones
```

2. Open R or RStudio and run the script:
```r
source("GSE58135_ER_Breast_Cancer_DESeq2_Analysis.R")
```

### What the Script Does

The analysis pipeline consists of 12 main steps:

1. **Download Data** - Fetches NCBI-generated raw counts from GEO (GSE58135)
2. **Fetch Metadata** - Downloads sample annotation from GEO
3. **Identify ER+ Samples** - Filters for ER+ tumors and adjacent normal tissue
4. **Match Tumor-Normal Pairs** - Extracts patient IDs from sample names
5. **Filter Count Table** - Subsets to paired samples only
6. **Create Sample Metadata** - Builds colData for DESeq2
7. **Apply Strict 1:1 Matching** - Keeps only patients with exactly one tumor and one normal sample
8. **Map Gene Symbols** - Converts gene symbols to Entrez IDs
9. **Run DESeq2** - Performs paired differential expression analysis (`~patient + condition`)
10. **Print Results** - Displays differential expression statistics
11. **Prepare Visualization Data** - VST-normalizes counts and formats data
12. **Generate Plots** - Creates 4 publication-ready figures

### Expected Runtime

- **Total runtime:** ~5-10 minutes (depending on internet speed)
- **Data download:** ~2-3 minutes
- **DESeq2 analysis:** ~2-3 minutes
- **Visualization:** ~1 minute

## Output

### Console Output
- Sample counts and filtering statistics
- Patient pairing verification
- Differential expression results table
- Progress messages throughout analysis

### Differential Expression Results Table
```
     symbol  baseMean log2FoldChange     lfcSE       pvalue         padj
     ACAT1  1718.53    -1.074         0.258    3.09e-05    1.07e-04
     OXCT1   753.99    -1.406         0.264    1.06e-07    5.75e-07
     SHMT1  3747.66    -1.309         0.284    3.94e-06    1.62e-05
     ...
```

### Visualizations

**Figure 1 (p1):** Paired line plot showing individual patient trajectories from normal to tumor tissue  
**Figure 2 (p2):** Boxplots with jittered points showing expression distributions  
**Figure 3 (p3):** Heatmap of VST-normalized expression (rows = genes, columns = samples)  
**Figure 4 (p4):** Log2 fold change barplot with error bars and significance indicators  

**Significance indicators:**
- `*` = adjusted p < 0.05
- `**` = adjusted p < 0.01
- `***` = adjusted p < 0.001

### Optional: Save Results
Uncomment the final lines of the script to save:
- `DESeq2_results_ER_breast_cancer.csv` - Results table
- `Figure1_paired_expression.pdf` - Paired line plot
- `Figure2_fold_change.pdf` - Fold change barplot

## Dataset Details

**GEO Accession:** [GSE58135](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE58135)  
**Title:** Breast Cancer RNA-seq  
**Platform:** Illumina HiSeq 2000  
**Samples Used:** 20 ER+ tumor samples and 20 matched adjacent normal samples (40 total)  
**Genome Build:** GRCh38.p13  
**Gene Annotation:** NCBI Homo sapiens Annotation Release 109.20190905

### Citation for Original Dataset
Varley KE, Myers RM. *Breast Cancer RNA-seq.* Gene Expression Omnibus (GEO) GSE58135. Published June 11, 2014.

## Methods Summary

Raw RNA-seq count data were obtained from NCBI-generated count matrices (GRCh38.p13). Strict 1:1 matching was applied to retain only patients with exactly one tumor and one matched normal sample (n=20 pairs). Genes with <10 counts in <2 samples were filtered. Differential expression analysis used DESeq2 with a paired design (`~patient + condition`) to control for patient-specific effects. Statistical significance was determined using Benjamini-Hochberg adjusted p-values.

## File Structure

```
.
├── README.md                                          # This file
├── GSE58135_ER_Breast_Cancer_DESeq2_Analysis.R       # Main analysis script
└── LICENSE                                            # License file (optional)
```

## Troubleshooting

### Common Issues

**Problem:** Package installation fails  
**Solution:** Make sure you have an active internet connection and try installing packages individually

**Problem:** "Cannot open URL" error during data download  
**Solution:** Check internet connection; GEO servers may be temporarily unavailable. Wait and retry.

**Problem:** "All samples have 0 counts" error  
**Solution:** This usually means column name mismatch. Ensure you're using the correct version of the script with fixed patient matching.

**Problem:** Memory issues  
**Solution:** Close other applications; the analysis requires ~2-4 GB RAM

## Contributing

Contributions are welcome! Please open an issue or submit a pull request.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Contact

For questions or issues, please open an issue on GitHub or contact [your email].

## Acknowledgments

- NCBI GEO for providing public access to the GSE58135 dataset
- Varley KE and Myers RM for generating and sharing the original data
- Bioconductor community for developing DESeq2 and related packages

## Citation

If you use this code in your research, please cite:

```
Weichhaus, M. (2026). Differential Expression Analysis of Metabolic Genes in ER+ Breast Cancer. 
GitHub repository: https://github.com/mweichhaus/Metabolomics-and-Ketones.git
```

And cite the original dataset:
```
Varley KE, Myers RM. Breast Cancer RNA-seq. Gene Expression Omnibus (GEO) GSE58135. 
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE58135
```