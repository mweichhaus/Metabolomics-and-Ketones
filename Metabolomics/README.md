# GC×GC-MS Metabolomics Analysis Pipeline

## Overview

This Python script performs comprehensive metabolomics analysis on GC×GC-MS data from breast cancer cell lines (MCF-7 and T47D) under different glucose and beta-hydroxybutyrate (BHB) treatment conditions.

## Treatment Groups

- **High Glucose (HiGlc)**: Standard DMEM medium
- **Low Glucose (LoGlc)**: 5% glucose (225 mg/L)
- **Low Glucose + BHB**: 5% glucose + 10 mM beta-hydroxybutyrate

## Analysis Methods

1. **Data Preprocessing**
   - Exclusion of solvents, contaminants, and derivatization byproducts
   - Metabolite identification and standard nomenclature mapping
   - Averaging of duplicate metabolite entries

2. **Fisher Ratio Analysis**
   - Calculates F-statistic: F = MS_between / MS_within
   - Compares to critical F-value at α = 0.05
   - Identifies metabolites with significant treatment effects

3. **Visualization**
   - Hierarchical clustering heatmaps (Z-score normalized)
   - PCA plots with 95% confidence ellipses
   - Bar plots with SEM error bars and significance markers

## Installation

```bash
# Create virtual environment (recommended)
python -m venv metabolomics_env
source metabolomics_env/bin/activate  # Linux/Mac
# or: metabolomics_env\Scripts\activate  # Windows

# Install dependencies
pip install -r requirements.txt
```

## Dependencies

- Python >= 3.10
- pandas >= 2.0
- numpy >= 1.24
- scipy >= 1.11
- matplotlib >= 3.7
- seaborn >= 0.12
- scikit-learn >= 1.3
- openpyxl >= 3.1

## Usage

```bash
python GCMS_Metabolomics_Analysis.py --mcf7 <mcf7_file.xlsx> --t47d <t47d_file.xlsx> --output <output_dir>
```

### Arguments

| Argument | Short | Description |
|----------|-------|-------------|
| `--mcf7` | `-m` | Path to MCF-7 Excel data file (required) |
| `--t47d` | `-t` | Path to T47D Excel data file (required) |
| `--output` | `-o` | Output directory (default: `output`) |

### Example

```bash
python GCMS_Metabolomics_Analysis.py \
  --mcf7 2021May_MCF7_All_sets_wout_HexaneMSTFA_Session_1_Discovery.xlsx \
  --t47d 2021May_T47D_All_sets_wout_HexaneMSTFA_Session_1_Discovery.xlsx \
  --output results/
```

## Output Files

| File | Description |
|------|-------------|
| `MCF7_Fisher_Ratio.csv` | Fisher Ratio statistics for MCF-7 |
| `T47D_Fisher_Ratio.csv` | Fisher Ratio statistics for T47D |
| `MCF7_statistics.csv` | Descriptive statistics (means, SEMs) for MCF-7 |
| `T47D_statistics.csv` | Descriptive statistics (means, SEMs) for T47D |
| `MCF7_heatmap.png` | Hierarchical clustering heatmap for MCF-7 |
| `T47D_heatmap.png` | Hierarchical clustering heatmap for T47D |
| `MCF7_PCA.png` | PCA scores and loadings for MCF-7 |
| `T47D_PCA.png` | PCA scores and loadings for T47D |
| `MCF7_barplot.png` | Bar plot with significance markers for MCF-7 |
| `T47D_barplot.png` | Bar plot with significance markers for T47D |
| `Supplementary_Table.png` | Combined statistics table for publication |

## Significance Markers

- `*` : p < 0.05 (statistically significant)
- `†` : p < 0.10 (trend, approaching significance)

## Data Format Requirements

The input Excel files should contain:
- Sheet named `MCF7 Compounds Table` or `T47D Compounds Table`
- Row 1: Treatment conditions
- Row 2: Set/batch numbers
- Row 6+: Compound data
- Column 3: Compound names
- Column 10+: Peak intensity values

## Configuration

Key parameters can be modified at the top of the script:

```python
# Experimental sets to include (default: 3, 4, 5)
SETS_TO_KEEP = [3, 4, 5]

# Compounds to exclude
EXCLUDE_LIST = ['Toluene', 'n-Hexane', ...]

# Metabolite name mapping
METABOLITE_MAP = {
    'L-Isoleucine, TMS derivative': 'Isoleucine',
    ...
}
```

## Citation

If you use this analysis pipeline, please cite:

> [Your publication details here]

## License

This code is provided for academic research purposes.

## Contact

Dr. Michael's Laboratory  
Chaminade University of Honolulu  
Laboratory of Molecular Cancer Research

---
Generated: January 2026
