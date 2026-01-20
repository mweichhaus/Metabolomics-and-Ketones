#!/usr/bin/env python3
"""
GC×GC-MS Metabolomics Analysis Pipeline
========================================

Analysis of breast cancer cell metabolomics data comparing MCF-7 and T47D cells
under different glucose and beta-hydroxybutyrate (BHB) treatment conditions.

Treatment Groups:
    - High Glucose (HiGlc): Standard DMEM medium
    - Low Glucose (LoGlc): 5% glucose (225 mg/L)
    - Low Glucose + BHB: 5% glucose + 10 mM beta-hydroxybutyrate

Analysis Methods:
    - Fisher Ratio analysis for biomarker identification
    - Principal Component Analysis (PCA)
    - Hierarchical clustering heatmaps
    - Bar plots with statistical significance markers

Author: Generated for Dr. Michael's Lab, Chaminade University of Honolulu
Date: January 2026
Version: 1.0

Dependencies:
    - Python >= 3.10
    - pandas >= 2.0
    - numpy >= 1.24
    - scipy >= 1.11
    - matplotlib >= 3.7
    - seaborn >= 0.12
    - scikit-learn >= 1.3
    - openpyxl >= 3.1 (for Excel file reading)

Usage:
    python GCMS_Metabolomics_Analysis.py --mcf7 <mcf7_file.xlsx> --t47d <t47d_file.xlsx> --output <output_dir>

"""

import argparse
import os
import sys
from pathlib import Path
from typing import Tuple, List, Dict, Optional

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from matplotlib.patches import Ellipse
from matplotlib.lines import Line2D
import warnings

warnings.filterwarnings('ignore')

# =============================================================================
# CONFIGURATION
# =============================================================================

# Metabolite name mapping: Raw GC-MS names -> Standard nomenclature
METABOLITE_MAP = {
    'L-Isoleucine, TMS derivative': 'Isoleucine',
    'Glycine, 3TMS derivative': 'Glycine',
    'DL-Phenylalanine, TMS derivative': 'Phenylalanine',
    'α-Hydroxyisobutyric acid, 2TMS derivative': 'α-Hydroxyisobutyrate',
    'Acetic acid, TMS derivative': 'Acetate',
    'Pyroglutamic acid, TMS derivative': 'Pyroglutamate',
    'Lactic Acid, 2TMS derivative': 'Lactate',
    'Propanedioic acid, 2TMS derivative': 'Malonate',
    'Ethanolamine, 3TMS derivative': 'Ethanolamine',
    'Glycerol, 3TMS derivative': 'Glycerol',
    'β-D-(-)-Ribopyranose, 4TMS derivative': 'Ribose',
}

# Compounds to exclude (solvents, contaminants, derivatization byproducts)
EXCLUDE_LIST = [
    'Toluene', 'n-Hexane', 'Hexane', 'Pentane', 'Cyclopentane', '1-Propanol',
    'Nonadecane', 'Heptacosane', 'N,N-Dimethyltrifluoroacetamide', 
    'Butylated Hydroxytoluene', '2,4,4,6-Tetramethyl', '2,4-Di-tert-butylphenol',
    'Methylphosphonic acid', '3-Buten-1-ol', 'Ethanol, 2-(trimethylsilyl)',
    '1,1-Dimethylethanol', 'Butyrolactone', 'Butanoic acid, 2-methylbutyl',
    'Acetaldehyde, O-methyloxime', 'chloro'
]

# Color scheme for treatment groups
COLORS = {
    'HiGlc': '#2ecc71',  # Green
    'LoGlc': '#e74c3c',  # Red
    'BHB': '#3498db'     # Blue
}

# Treatment labels for publication
TREATMENT_LABELS = {
    'HiGlc': 'High Glucose',
    'LoGlc': 'Low Glucose',
    'BHB': 'Low Glucose + BHB'
}

# Sets to include in analysis (excluding sets with batch effects)
SETS_TO_KEEP = [3, 4, 5]


# =============================================================================
# DATA PROCESSING FUNCTIONS
# =============================================================================

def should_exclude(compound_name: str) -> bool:
    """
    Check if a compound should be excluded from analysis.
    
    Parameters
    ----------
    compound_name : str
        Name of the compound from GC-MS data
        
    Returns
    -------
    bool
        True if compound should be excluded
    """
    compound_lower = str(compound_name).lower()
    for excl in EXCLUDE_LIST:
        if excl.lower() in compound_lower:
            return True
    return False


def get_biological_name(raw_name: str) -> Optional[str]:
    """
    Convert raw GC-MS compound name to standard biological nomenclature.
    
    Parameters
    ----------
    raw_name : str
        Raw compound name from GC-MS analysis
        
    Returns
    -------
    Optional[str]
        Standardized metabolite name, or None if not a target metabolite
    """
    name_str = str(raw_name).strip()
    for key, value in METABOLITE_MAP.items():
        if key in name_str:
            return value
    return None


def load_and_process_data(
    filepath: str,
    sheet_name: str,
    sets_to_keep: List[int] = SETS_TO_KEEP
) -> Tuple[np.ndarray, List[str], List[str], List[int]]:
    """
    Load GC-MS data from Excel file and process for analysis.
    
    Parameters
    ----------
    filepath : str
        Path to Excel file containing GC-MS data
    sheet_name : str
        Name of sheet containing compound data
    sets_to_keep : List[int]
        Experimental sets to include (default: [3, 4, 5])
        
    Returns
    -------
    Tuple containing:
        - data_matrix : np.ndarray
            Matrix of metabolite intensities (metabolites × samples)
        - metabolite_names : List[str]
            List of metabolite names
        - treatments : List[str]
            Treatment condition for each sample
        - sets : List[int]
            Experimental set number for each sample
    """
    # Read Excel file
    df = pd.read_excel(filepath, sheet_name=sheet_name)
    
    # Extract metadata and data
    treatments_raw = df.iloc[0, 9:].tolist()
    sets_raw = df.iloc[1, 9:].tolist()
    compounds = df.iloc[5:, 2].values
    data = df.iloc[5:, 9:].values
    
    # Filter samples by set and treatment
    sample_mask = []
    sample_treatments = []
    sample_sets = []
    
    for i, (trt, s) in enumerate(zip(treatments_raw, sets_raw)):
        try:
            set_num = int(s)
        except (ValueError, TypeError):
            sample_mask.append(False)
            continue
            
        if set_num not in sets_to_keep:
            sample_mask.append(False)
            continue
        
        # Parse treatment condition
        trt_str = str(trt)
        if 'A:' in trt_str or '(+) Glucose' in trt_str:
            trt_short = 'HiGlc'
        elif 'B:' in trt_str or ('(-) Glucose' in trt_str and 'BHB' not in trt_str):
            trt_short = 'LoGlc'
        elif 'C:' in trt_str or 'BHB' in trt_str:
            trt_short = 'BHB'
        else:
            sample_mask.append(False)
            continue
        
        sample_mask.append(True)
        sample_treatments.append(trt_short)
        sample_sets.append(set_num)
    
    sample_mask = np.array(sample_mask)
    
    # Process metabolites
    metabolite_data = {}
    
    for i, compound in enumerate(compounds):
        if should_exclude(compound):
            continue
            
        bio_name = get_biological_name(compound)
        if bio_name is None:
            continue
        
        row_data = data[i, sample_mask].astype(float)
        
        # Aggregate duplicate metabolite entries
        if bio_name not in metabolite_data:
            metabolite_data[bio_name] = []
        metabolite_data[bio_name].append(row_data)
    
    # Average duplicate entries
    final_names = []
    final_data = []
    
    for name, data_list in metabolite_data.items():
        stacked = np.vstack(data_list)
        averaged = np.nanmean(stacked, axis=0)
        final_names.append(name)
        final_data.append(averaged)
    
    return np.array(final_data), final_names, sample_treatments, sample_sets


# =============================================================================
# STATISTICAL ANALYSIS FUNCTIONS
# =============================================================================

def compute_fisher_ratio(
    data: np.ndarray,
    names: List[str],
    treatments: List[str],
    sets: List[int]
) -> pd.DataFrame:
    """
    Compute Fisher Ratio statistics for each metabolite.
    
    The Fisher Ratio (F-statistic) is calculated as:
        F = MS_between / MS_within
    
    where MS_between is the mean square between groups and MS_within is the
    mean square within groups.
    
    Parameters
    ----------
    data : np.ndarray
        Data matrix (metabolites × samples)
    names : List[str]
        Metabolite names
    treatments : List[str]
        Treatment condition for each sample
    sets : List[int]
        Experimental set for each sample
        
    Returns
    -------
    pd.DataFrame
        DataFrame with Fisher Ratio statistics for each metabolite
    """
    results = []
    
    for i, name in enumerate(names):
        row = data[i]
        
        # Group samples by treatment
        groups = {'HiGlc': [], 'LoGlc': [], 'BHB': []}
        
        for val, trt in zip(row, treatments):
            if not np.isnan(val) and val > 0:
                groups[trt].append(val)
        
        # Skip if insufficient data
        if any(len(g) < 3 for g in groups.values()):
            continue
        
        # Calculate within-group statistics
        within_variances = []
        group_means = {}
        group_ns = {}
        
        for grp_name, vals in groups.items():
            var = np.var(vals, ddof=1)
            within_variances.append(var)
            group_means[grp_name] = np.mean(vals)
            group_ns[grp_name] = len(vals)
        
        # Calculate overall statistics
        all_vals = []
        for vals in groups.values():
            all_vals.extend(vals)
        
        grand_mean = np.mean(all_vals)
        n_total = len(all_vals)
        k = 3  # Number of groups
        
        # Between-group mean square
        ss_between = sum(
            group_ns[g] * (group_means[g] - grand_mean)**2 
            for g in groups
        )
        df_between = k - 1
        ms_between = ss_between / df_between
        
        # Within-group mean square
        ss_within = sum(
            (group_ns[g] - 1) * within_variances[i] 
            for i, g in enumerate(groups)
        )
        df_within = n_total - k
        ms_within = ss_within / df_within if df_within > 0 else 0
        
        # Fisher Ratio (F-statistic)
        fisher_ratio = ms_between / ms_within if ms_within > 0 else np.nan
        
        # Critical F-value and p-value
        f_crit = stats.f.ppf(0.95, df_between, df_within)
        p_value = 1 - stats.f.cdf(fisher_ratio, df_between, df_within) \
                  if not np.isnan(fisher_ratio) else 1.0
        
        significant = fisher_ratio > f_crit if not np.isnan(fisher_ratio) else False
        
        results.append({
            'Compound': name,
            'n_HiGlc': group_ns['HiGlc'],
            'n_LoGlc': group_ns['LoGlc'],
            'n_BHB': group_ns['BHB'],
            'n_total': n_total,
            'Mean_HiGlc': group_means['HiGlc'],
            'Mean_LoGlc': group_means['LoGlc'],
            'Mean_BHB': group_means['BHB'],
            'Var_HiGlc': within_variances[0],
            'Var_LoGlc': within_variances[1],
            'Var_BHB': within_variances[2],
            'MS_Between': ms_between,
            'MS_Within': ms_within,
            'Fisher_Ratio': fisher_ratio,
            'df_between': df_between,
            'df_within': df_within,
            'F_crit_0.05': f_crit,
            'p_value': p_value,
            'Significant': significant
        })
    
    return pd.DataFrame(results).sort_values('Fisher_Ratio', ascending=False)


def compute_descriptive_statistics(
    data: np.ndarray,
    names: List[str],
    treatments: List[str],
    sets: List[int],
    fisher_df: pd.DataFrame
) -> pd.DataFrame:
    """
    Compute descriptive statistics (means, SEMs) for publication.
    
    Parameters
    ----------
    data : np.ndarray
        Data matrix (metabolites × samples)
    names : List[str]
        Metabolite names
    treatments : List[str]
        Treatment conditions
    sets : List[int]
        Experimental sets
    fisher_df : pd.DataFrame
        Fisher Ratio results
        
    Returns
    -------
    pd.DataFrame
        DataFrame with means, SEMs, and p-values
    """
    results = []
    
    for i, name in enumerate(names):
        row = data[i]
        
        # Average technical replicates to get biological replicates
        bio_reps = {'HiGlc': [], 'LoGlc': [], 'BHB': []}
        
        for s in sorted(set(sets)):
            for trt in ['HiGlc', 'LoGlc', 'BHB']:
                idx = [j for j, (t, st) in enumerate(zip(treatments, sets)) 
                       if t == trt and st == s]
                if len(idx) > 0:
                    bio_reps[trt].append(np.nanmean([row[j] for j in idx]))
        
        # Calculate means and SEMs
        means = {trt: np.mean(vals) for trt, vals in bio_reps.items()}
        sems = {trt: stats.sem(vals) if len(vals) > 1 else 0 
                for trt, vals in bio_reps.items()}
        
        # Get Fisher Ratio p-value
        fisher_row = fisher_df[fisher_df['Compound'] == name]
        fisher_p = fisher_row['p_value'].values[0] if len(fisher_row) > 0 else 1.0
        
        results.append({
            'name': name,
            'mean_HiGlc': means['HiGlc'],
            'mean_LoGlc': means['LoGlc'],
            'mean_BHB': means['BHB'],
            'sem_HiGlc': sems['HiGlc'],
            'sem_LoGlc': sems['LoGlc'],
            'sem_BHB': sems['BHB'],
            'fisher_p': fisher_p,
            'n': len(bio_reps['HiGlc'])
        })
    
    return pd.DataFrame(results).sort_values('fisher_p')


# =============================================================================
# VISUALIZATION FUNCTIONS
# =============================================================================

def create_heatmap(
    data: np.ndarray,
    names: List[str],
    treatments: List[str],
    sets: List[int],
    cell_line: str,
    output_path: str
) -> None:
    """
    Create publication-ready heatmap with hierarchical clustering.
    
    Parameters
    ----------
    data : np.ndarray
        Data matrix (metabolites × samples)
    names : List[str]
        Metabolite names
    treatments : List[str]
        Treatment conditions
    sets : List[int]
        Experimental sets
    cell_line : str
        Cell line name for title
    output_path : str
        Path to save figure
    """
    # Average technical replicates
    unique_bio_reps = sorted(
        set(zip(treatments, sets)),
        key=lambda x: (0 if x[0]=='HiGlc' else 1 if x[0]=='LoGlc' else 2, x[1])
    )
    
    avg_data = []
    group_labels = []
    group_colors = []
    
    for trt, s in unique_bio_reps:
        idx = [i for i, (t, st) in enumerate(zip(treatments, sets)) 
               if t == trt and st == s]
        avg_data.append(np.nanmean(data[:, idx], axis=1))
        
        rep_num = [s for s in sorted(set(sets))].index(s) % 3 + 1
        group_labels.append(f"{TREATMENT_LABELS[trt][:4]}-{rep_num}")
        group_colors.append(COLORS[trt])
    
    avg_matrix = np.column_stack(avg_data)
    
    # Log2 transform and Z-score normalize
    avg_log = np.log2(np.where(avg_matrix > 0, avg_matrix, np.nan))
    avg_zscore = (avg_log - np.nanmean(avg_log, axis=1, keepdims=True)) / \
                 np.nanstd(avg_log, axis=1, keepdims=True)
    avg_zscore = np.nan_to_num(avg_zscore, 0)
    
    # Create DataFrame
    df = pd.DataFrame(avg_zscore, index=names, columns=group_labels)
    
    # Create figure
    fig_height = max(5, len(names) * 0.5)
    g = sns.clustermap(
        df,
        cmap='RdBu_r', center=0, vmin=-2.5, vmax=2.5,
        col_colors=[group_colors],
        row_cluster=True, col_cluster=False,
        dendrogram_ratio=(0.15, 0.02),
        cbar_pos=(0.02, 0.8, 0.03, 0.15),
        figsize=(10, fig_height),
        xticklabels=True, yticklabels=True,
        linewidths=0.5
    )
    
    g.ax_heatmap.set_xlabel('')
    g.ax_heatmap.set_ylabel('')
    g.ax_heatmap.tick_params(axis='x', labelsize=9, rotation=45)
    g.ax_heatmap.tick_params(axis='y', labelsize=10)
    
    # Add treatment dividers
    n_per_group = 3
    g.ax_heatmap.axvline(x=n_per_group, color='black', linewidth=2)
    g.ax_heatmap.axvline(x=n_per_group*2, color='black', linewidth=2)
    
    g.ax_cbar.set_ylabel('Z-score', fontsize=10)
    g.fig.suptitle(f'{cell_line}', fontsize=14, fontweight='bold', y=1.02)
    
    # Legend
    legend_elements = [
        Line2D([0], [0], marker='s', color='w', markerfacecolor=COLORS['HiGlc'],
               markersize=12, label='High Glucose'),
        Line2D([0], [0], marker='s', color='w', markerfacecolor=COLORS['LoGlc'],
               markersize=12, label='Low Glucose'),
        Line2D([0], [0], marker='s', color='w', markerfacecolor=COLORS['BHB'],
               markersize=12, label='Low Glucose + BHB'),
    ]
    g.ax_heatmap.legend(
        handles=legend_elements, loc='upper left',
        bbox_to_anchor=(1.15, 1.0), fontsize=9, frameon=True
    )
    
    g.fig.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()


def create_pca_plot(
    data: np.ndarray,
    names: List[str],
    treatments: List[str],
    sets: List[int],
    cell_line: str,
    output_path: str
) -> None:
    """
    Create publication-ready PCA plot with scores and loadings.
    
    Parameters
    ----------
    data : np.ndarray
        Data matrix (metabolites × samples)
    names : List[str]
        Metabolite names
    treatments : List[str]
        Treatment conditions
    sets : List[int]
        Experimental sets
    cell_line : str
        Cell line name for title
    output_path : str
        Path to save figure
    """
    # Average technical replicates
    unique_groups = sorted(
        set(zip(treatments, sets)),
        key=lambda x: (0 if x[0]=='HiGlc' else 1 if x[0]=='LoGlc' else 2, x[1])
    )
    
    avg_data = []
    group_info = []
    
    for trt, s in unique_groups:
        idx = [i for i, (t, st) in enumerate(zip(treatments, sets)) 
               if t == trt and st == s]
        avg_data.append(np.nanmean(data[:, idx], axis=1))
        group_info.append({'trt': trt, 'set': s})
    
    X = np.log2(np.column_stack(avg_data).T)
    X = np.nan_to_num(X, nan=np.nanmean(X))
    
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    pca = PCA(n_components=min(len(X), len(names)))
    X_pca = pca.fit_transform(X_scaled)
    exp_var = pca.explained_variance_ratio_ * 100
    
    # Create figure
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # Scores plot
    for i, info in enumerate(group_info):
        ax1.scatter(
            X_pca[i, 0], X_pca[i, 1],
            c=COLORS[info['trt']], marker='o',
            s=120, edgecolors='black', linewidths=1, zorder=3
        )
    
    # Add confidence ellipses
    for trt in ['HiGlc', 'LoGlc', 'BHB']:
        idx = [i for i, info in enumerate(group_info) if info['trt'] == trt]
        if len(idx) >= 2:
            pts = X_pca[idx, :2]
            center = pts.mean(axis=0)
            
            cov = np.cov(pts.T)
            eigenvalues, eigenvectors = np.linalg.eigh(cov)
            angle = np.degrees(np.arctan2(eigenvectors[1, 1], eigenvectors[0, 1]))
            width, height = 2 * np.sqrt(eigenvalues) * 1.96  # 95% CI
            
            ellipse = Ellipse(
                center, width, height, angle=angle,
                facecolor=COLORS[trt], alpha=0.15,
                edgecolor=COLORS[trt], linewidth=2, zorder=1
            )
            ax1.add_patch(ellipse)
    
    ax1.axhline(y=0, color='gray', linestyle='--', alpha=0.4, linewidth=0.8)
    ax1.axvline(x=0, color='gray', linestyle='--', alpha=0.4, linewidth=0.8)
    ax1.set_xlabel(f'PC1 ({exp_var[0]:.1f}%)', fontsize=12)
    ax1.set_ylabel(f'PC2 ({exp_var[1]:.1f}%)', fontsize=12)
    ax1.set_title(f'{cell_line} - PCA Scores', fontsize=13, fontweight='bold')
    
    # Legend
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor=COLORS['HiGlc'],
               markersize=10, markeredgecolor='black', label='High Glucose'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor=COLORS['LoGlc'],
               markersize=10, markeredgecolor='black', label='Low Glucose'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor=COLORS['BHB'],
               markersize=10, markeredgecolor='black', label='Low Glucose + BHB'),
    ]
    ax1.legend(handles=legend_elements, loc='best', fontsize=9, frameon=True)
    
    # Loadings plot
    loadings = pca.components_[:2].T
    max_loading = np.abs(loadings).max()
    
    for i, name in enumerate(names):
        ax2.arrow(
            0, 0, loadings[i, 0], loadings[i, 1],
            head_width=0.03, head_length=0.02,
            fc='#34495e', ec='#34495e', alpha=0.8, linewidth=1.5
        )
        ax2.text(
            loadings[i, 0]*1.12, loadings[i, 1]*1.12,
            name, fontsize=9, ha='center', va='center', fontweight='medium'
        )
    
    ax2.axhline(y=0, color='gray', linestyle='--', alpha=0.4, linewidth=0.8)
    ax2.axvline(x=0, color='gray', linestyle='--', alpha=0.4, linewidth=0.8)
    ax2.set_xlabel(f'PC1 ({exp_var[0]:.1f}%)', fontsize=12)
    ax2.set_ylabel(f'PC2 ({exp_var[1]:.1f}%)', fontsize=12)
    ax2.set_title(f'{cell_line} - PCA Loadings', fontsize=13, fontweight='bold')
    
    # Set axis limits
    max_range = max(np.abs(X_pca[:, :2]).max(), max_loading) * 1.3
    ax1.set_xlim(-max_range * 1.2, max_range * 1.2)
    ax1.set_ylim(-max_range * 1.2, max_range * 1.2)
    ax2.set_xlim(-max_loading * 1.4, max_loading * 1.4)
    ax2.set_ylim(-max_loading * 1.4, max_loading * 1.4)
    
    plt.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()


def create_barplot(
    stats_df: pd.DataFrame,
    cell_line: str,
    output_path: str
) -> None:
    """
    Create publication-ready bar plot with error bars and significance markers.
    
    Parameters
    ----------
    stats_df : pd.DataFrame
        Descriptive statistics DataFrame
    cell_line : str
        Cell line name for title
    output_path : str
        Path to save figure
    """
    n_mets = len(stats_df)
    fig, ax = plt.subplots(figsize=(12, 6))
    
    x = np.arange(n_mets)
    width = 0.25
    
    # Plot bars
    ax.bar(
        x - width, stats_df['mean_HiGlc']/1000, width,
        yerr=stats_df['sem_HiGlc']/1000, capsize=3,
        label='High Glucose', color=COLORS['HiGlc'],
        edgecolor='black', linewidth=0.8, error_kw={'linewidth': 1}
    )
    
    ax.bar(
        x, stats_df['mean_LoGlc']/1000, width,
        yerr=stats_df['sem_LoGlc']/1000, capsize=3,
        label='Low Glucose', color=COLORS['LoGlc'],
        edgecolor='black', linewidth=0.8, error_kw={'linewidth': 1}
    )
    
    ax.bar(
        x + width, stats_df['mean_BHB']/1000, width,
        yerr=stats_df['sem_BHB']/1000, capsize=3,
        label='Low Glucose + BHB', color=COLORS['BHB'],
        edgecolor='black', linewidth=0.8, error_kw={'linewidth': 1}
    )
    
    # Add significance markers
    for i, (idx, row) in enumerate(stats_df.iterrows()):
        max_val = max(row['mean_HiGlc'], row['mean_LoGlc'], row['mean_BHB'])/1000
        max_sem = max(row['sem_HiGlc'], row['sem_LoGlc'], row['sem_BHB'])/1000
        y_pos = max_val + max_sem + 0.1
        
        if row['fisher_p'] < 0.05:
            marker = '*'
        elif row['fisher_p'] < 0.10:
            marker = '†'
        else:
            marker = ''
        
        if marker:
            ax.text(i, y_pos, marker, ha='center', va='bottom', 
                   fontsize=16, fontweight='bold')
    
    ax.set_xticks(x)
    ax.set_xticklabels(stats_df['name'], rotation=45, ha='right', fontsize=11)
    ax.set_ylabel('Peak Intensity (×10³)', fontsize=12)
    ax.set_title(f'{cell_line}', fontsize=14, fontweight='bold')
    ax.legend(loc='upper right', fontsize=10, frameon=True)
    ax.set_ylim(bottom=0)
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    plt.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()


def create_supplementary_table(
    mcf7_stats: pd.DataFrame,
    t47d_stats: pd.DataFrame,
    output_path: str
) -> None:
    """
    Create supplementary statistics table as image.
    
    Parameters
    ----------
    mcf7_stats : pd.DataFrame
        MCF-7 statistics
    t47d_stats : pd.DataFrame
        T47D statistics
    output_path : str
        Path to save figure
    """
    fig, ax = plt.subplots(figsize=(14, 8))
    ax.axis('off')
    
    def format_row(row, cell_line):
        sig = '*' if row['fisher_p'] < 0.05 else '†' if row['fisher_p'] < 0.10 else ''
        return [
            cell_line,
            row['name'],
            f"{row['mean_HiGlc']:.0f} ± {row['sem_HiGlc']:.0f}",
            f"{row['mean_LoGlc']:.0f} ± {row['sem_LoGlc']:.0f}",
            f"{row['mean_BHB']:.0f} ± {row['sem_BHB']:.0f}",
            f"{row['fisher_p']:.4f}",
            sig
        ]
    
    table_data = []
    for _, row in mcf7_stats.iterrows():
        table_data.append(format_row(row, 'MCF-7'))
    table_data.append(['', '', '', '', '', '', ''])
    for _, row in t47d_stats.iterrows():
        table_data.append(format_row(row, 'T47D'))
    
    col_labels = [
        'Cell Line', 'Metabolite', 'High Glucose\n(Mean ± SEM)',
        'Low Glucose\n(Mean ± SEM)', 'Low Glc + BHB\n(Mean ± SEM)',
        'p-value', ''
    ]
    
    table = ax.table(
        cellText=table_data, colLabels=col_labels, loc='center',
        cellLoc='center', colColours=['#e8f4f8']*7
    )
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1.2, 1.8)
    
    # Style header
    for i in range(7):
        table[(0, i)].set_fontsize(10)
        table[(0, i)].set_text_props(weight='bold')
    
    # Highlight significant rows
    for i, row_data in enumerate(table_data, 1):
        if row_data[6] in ['†', '*']:
            for j in range(7):
                table[(i, j)].set_facecolor('#fff3cd')
    
    plt.title(
        'Supplementary Table: Metabolite Quantification\n'
        '(Mean ± SEM; † p<0.10, * p<0.05, one-way ANOVA)',
        fontsize=12, fontweight='bold', pad=20
    )
    
    fig.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()


# =============================================================================
# MAIN ANALYSIS PIPELINE
# =============================================================================

def run_analysis(
    mcf7_filepath: str,
    t47d_filepath: str,
    output_dir: str
) -> None:
    """
    Run complete metabolomics analysis pipeline.
    
    Parameters
    ----------
    mcf7_filepath : str
        Path to MCF-7 data file
    t47d_filepath : str
        Path to T47D data file
    output_dir : str
        Directory for output files
    """
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    print("="*70)
    print("GC×GC-MS METABOLOMICS ANALYSIS PIPELINE")
    print("="*70)
    
    # Load and process data
    print("\n[1/6] Loading and processing data...")
    
    mcf7_data, mcf7_names, mcf7_trt, mcf7_sets = load_and_process_data(
        mcf7_filepath, 'MCF7 Compounds Table'
    )
    print(f"  MCF-7: {len(mcf7_names)} metabolites, {len(mcf7_trt)} samples")
    
    t47d_data, t47d_names, t47d_trt, t47d_sets = load_and_process_data(
        t47d_filepath, 'T47D Compounds Table'
    )
    print(f"  T47D: {len(t47d_names)} metabolites, {len(t47d_trt)} samples")
    
    # Compute Fisher Ratios
    print("\n[2/6] Computing Fisher Ratio statistics...")
    
    mcf7_fisher = compute_fisher_ratio(mcf7_data, mcf7_names, mcf7_trt, mcf7_sets)
    t47d_fisher = compute_fisher_ratio(t47d_data, t47d_names, t47d_trt, t47d_sets)
    
    mcf7_fisher.to_csv(os.path.join(output_dir, 'MCF7_Fisher_Ratio.csv'), index=False)
    t47d_fisher.to_csv(os.path.join(output_dir, 'T47D_Fisher_Ratio.csv'), index=False)
    
    print(f"  MCF-7 significant (p<0.10): {(mcf7_fisher['p_value'] < 0.10).sum()}")
    print(f"  T47D significant (p<0.10): {(t47d_fisher['p_value'] < 0.10).sum()}")
    
    # Compute descriptive statistics
    print("\n[3/6] Computing descriptive statistics...")
    
    mcf7_stats = compute_descriptive_statistics(
        mcf7_data, mcf7_names, mcf7_trt, mcf7_sets, mcf7_fisher
    )
    t47d_stats = compute_descriptive_statistics(
        t47d_data, t47d_names, t47d_trt, t47d_sets, t47d_fisher
    )
    
    mcf7_stats.to_csv(os.path.join(output_dir, 'MCF7_statistics.csv'), index=False)
    t47d_stats.to_csv(os.path.join(output_dir, 'T47D_statistics.csv'), index=False)
    
    # Create visualizations
    print("\n[4/6] Creating heatmaps...")
    
    create_heatmap(
        mcf7_data, mcf7_names, mcf7_trt, mcf7_sets,
        'MCF-7', os.path.join(output_dir, 'MCF7_heatmap.png')
    )
    create_heatmap(
        t47d_data, t47d_names, t47d_trt, t47d_sets,
        'T47D', os.path.join(output_dir, 'T47D_heatmap.png')
    )
    
    print("\n[5/6] Creating PCA plots...")
    
    create_pca_plot(
        mcf7_data, mcf7_names, mcf7_trt, mcf7_sets,
        'MCF-7', os.path.join(output_dir, 'MCF7_PCA.png')
    )
    create_pca_plot(
        t47d_data, t47d_names, t47d_trt, t47d_sets,
        'T47D', os.path.join(output_dir, 'T47D_PCA.png')
    )
    
    print("\n[6/6] Creating bar plots and supplementary table...")
    
    create_barplot(mcf7_stats, 'MCF-7', os.path.join(output_dir, 'MCF7_barplot.png'))
    create_barplot(t47d_stats, 'T47D', os.path.join(output_dir, 'T47D_barplot.png'))
    
    create_supplementary_table(
        mcf7_stats, t47d_stats,
        os.path.join(output_dir, 'Supplementary_Table.png')
    )
    
    # Print summary
    print("\n" + "="*70)
    print("ANALYSIS COMPLETE")
    print("="*70)
    
    print("\nKey Findings:")
    print("-"*40)
    
    print("\nMCF-7 (ranked by significance):")
    for _, row in mcf7_fisher.head(3).iterrows():
        sig = '†' if row['p_value'] < 0.10 else ''
        print(f"  {row['Compound']:<22} F={row['Fisher_Ratio']:.3f}  p={row['p_value']:.4f} {sig}")
    
    print("\nT47D (ranked by significance):")
    for _, row in t47d_fisher.head(3).iterrows():
        sig = '†' if row['p_value'] < 0.10 else ''
        print(f"  {row['Compound']:<22} F={row['Fisher_Ratio']:.3f}  p={row['p_value']:.4f} {sig}")
    
    print(f"\nOutput files saved to: {output_dir}")
    print("\nFiles generated:")
    for f in sorted(os.listdir(output_dir)):
        print(f"  - {f}")


# =============================================================================
# COMMAND LINE INTERFACE
# =============================================================================

def main():
    """Command line entry point."""
    parser = argparse.ArgumentParser(
        description='GC×GC-MS Metabolomics Analysis Pipeline',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python GCMS_Metabolomics_Analysis.py --mcf7 MCF7_data.xlsx --t47d T47D_data.xlsx --output results/
  python GCMS_Metabolomics_Analysis.py -m MCF7.xlsx -t T47D.xlsx -o output/
        """
    )
    
    parser.add_argument(
        '-m', '--mcf7',
        required=True,
        help='Path to MCF-7 Excel data file'
    )
    parser.add_argument(
        '-t', '--t47d',
        required=True,
        help='Path to T47D Excel data file'
    )
    parser.add_argument(
        '-o', '--output',
        default='output',
        help='Output directory (default: output)'
    )
    
    args = parser.parse_args()
    
    # Validate input files
    if not os.path.exists(args.mcf7):
        print(f"Error: MCF-7 file not found: {args.mcf7}")
        sys.exit(1)
    if not os.path.exists(args.t47d):
        print(f"Error: T47D file not found: {args.t47d}")
        sys.exit(1)
    
    # Run analysis
    run_analysis(args.mcf7, args.t47d, args.output)


if __name__ == '__main__':
    main()
