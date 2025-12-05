#!/usr/bin/env python3
"""
Generate Supplementary Figure S1: TPR/TNR distributions for GARFIELD-NGS scores
Shows true variants (TPR, red) and false variants (TNR, green) identified by score
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc
import seaborn as sns

def calculate_tpr_tnr_at_threshold(y_true, y_scores, threshold):
    """Calculate TPR and TNR at a specific threshold"""
    predictions = y_scores >= threshold
    
    # True Positives: correctly identified true variants
    tp = np.sum((predictions == 1) & (y_true == 1))
    # False Negatives: true variants called false
    fn = np.sum((predictions == 0) & (y_true == 1))
    # True Negatives: correctly identified false variants
    tn = np.sum((predictions == 0) & (y_true == 0))
    # False Positives: false variants called true
    fp = np.sum((predictions == 1) & (y_true == 0))
    
    tpr = tp / (tp + fn) if (tp + fn) > 0 else 0  # Sensitivity / Recall
    tnr = tn / (tn + fp) if (tn + fp) > 0 else 0  # Specificity
    
    return tpr, tnr

def find_threshold_for_tnr(y_true, y_scores, target_tnr=0.90):
    """Find threshold that gives approximately target TNR"""
    thresholds = np.linspace(y_scores.min(), y_scores.max(), 1000)
    
    best_threshold = None
    best_diff = float('inf')
    best_tpr = 0
    
    for thresh in thresholds:
        tpr, tnr = calculate_tpr_tnr_at_threshold(y_true, y_scores, thresh)
        diff = abs(tnr - target_tnr)
        
        if diff < best_diff:
            best_diff = diff
            best_threshold = thresh
            best_tpr = tpr
    
    return best_threshold, best_tpr

def find_max_accuracy_threshold(y_true, y_scores):
    """Find threshold that maximizes accuracy"""
    thresholds = np.linspace(y_scores.min(), y_scores.max(), 1000)
    
    best_threshold = None
    best_accuracy = 0
    best_tpr = 0
    best_tnr = 0
    
    for thresh in thresholds:
        predictions = y_scores >= thresh
        accuracy = np.mean(predictions == y_true)
        
        if accuracy > best_accuracy:
            best_accuracy = accuracy
            best_threshold = thresh
            best_tpr, best_tnr = calculate_tpr_tnr_at_threshold(y_true, y_scores, thresh)
    
    return best_threshold, best_tpr, best_tnr, best_accuracy

def plot_tpr_tnr_distribution(y_true, y_scores, variant_type, ax):
    """
    Plot TPR/TNR distribution for a single variant type
    
    Args:
        y_true: Array of true labels (1=true variant, 0=false variant)
        y_scores: Array of GARFIELD-NGS scores
        variant_type: String like "Illumina INS/DELs"
        ax: Matplotlib axis object
    """
    
    # Separate scores by true/false variants
    true_variant_scores = y_scores[y_true == 1]
    false_variant_scores = y_scores[y_true == 0]
    
    # Find max accuracy threshold
    max_acc_thresh, tpr_max, tnr_max, accuracy = find_max_accuracy_threshold(y_true, y_scores)
    
    # Find TPR at 90% and 95% TNR
    thresh_90, tpr_90 = find_threshold_for_tnr(y_true, y_scores, 0.90)
    thresh_95, tpr_95 = find_threshold_for_tnr(y_true, y_scores, 0.95)
    
    # Plot histograms
    ax.hist(true_variant_scores, bins=50, alpha=0.6, color='red', 
            label=f'True variants (n={len(true_variant_scores)})', density=True)
    ax.hist(false_variant_scores, bins=50, alpha=0.6, color='green', 
            label=f'False variants (n={len(false_variant_scores)})', density=True)
    
    # Add vertical line for max accuracy threshold
    ax.axvline(max_acc_thresh, color='black', linestyle='--', linewidth=2,
               label=f'Max accuracy threshold')
    
    # Add title and labels
    ax.set_xlabel('GARFIELD-NGS Score', fontsize=12)
    ax.set_ylabel('Density', fontsize=12)
    ax.set_title(variant_type, fontsize=14, fontweight='bold')
    ax.legend(loc='upper right', fontsize=9)
    
    # Add text box with statistics
    stats_text = (
        f'Max Accuracy Threshold: {max_acc_thresh:.3f}\n'
        f'TPR: {tpr_max:.1%} (true retained)\n'
        f'TNR: {tnr_max:.1%} (false filtered)\n'
        f'\n'
        f'At 90% TNR: TPR = {tpr_90:.1%}\n'
        f'At 95% TNR: TPR = {tpr_95:.1%}'
    )
    
    ax.text(0.02, 0.98, stats_text,
            transform=ax.transAxes,
            fontsize=9,
            verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    ax.grid(True, alpha=0.3)
    
    return {
        'max_acc_threshold': max_acc_thresh,
        'tpr_at_max_acc': tpr_max,
        'tnr_at_max_acc': tnr_max,
        'accuracy': accuracy,
        'tpr_at_90_tnr': tpr_90,
        'tpr_at_95_tnr': tpr_95
    }

def generate_figure_s1(results_dict, output_file='figure_S1.png'):
    """
    Generate complete Figure S1 with 4 subplots
    
    Args:
        results_dict: Dictionary with keys:
            - 'illumina_indels': (y_true, y_scores)
            - 'illumina_snvs': (y_true, y_scores)
            - 'ion_indels': (y_true, y_scores)
            - 'ion_snvs': (y_true, y_scores)
        output_file: Path to save figure
    """
    
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle('Supplementary Figure S1: TPR/TNR Distributions for GARFIELD-NGS Score',
                 fontsize=16, fontweight='bold', y=0.995)
    
    # Subplot titles
    variant_types = [
        ('illumina_indels', 'a) Illumina INS/DELs'),
        ('illumina_snvs', 'b) Illumina SNVs'),
        ('ion_indels', 'c) ION INS/DELs'),
        ('ion_snvs', 'd) ION SNVs')
    ]
    
    stats_summary = {}
    
    for idx, (key, title) in enumerate(variant_types):
        row = idx // 2
        col = idx % 2
        
        if key in results_dict:
            y_true, y_scores = results_dict[key]
            stats = plot_tpr_tnr_distribution(y_true, y_scores, title, axes[row, col])
            stats_summary[key] = stats
        else:
            axes[row, col].text(0.5, 0.5, f'No data for {title}',
                               ha='center', va='center', fontsize=14)
            axes[row, col].set_title(title, fontsize=14, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Figure saved to: {output_file}")
    
    return stats_summary

def main():
    """
    Example usage - replace with your actual data
    """
    
    # Example: Load your predictions and true labels
    # You'll need to replace this with your actual data loading
    
    # Format: DataFrame with columns ['y_true', 'y_score', 'variant_type', 'platform']
    # y_true: 1 = true variant (in GIAB), 0 = false variant
    # y_score: GARFIELD-NGS prediction score (0-1)
    
    # Example data structure:
    # df = pd.read_csv('garfield_predictions.csv')
    # illumina_indels = df[(df['platform'] == 'Illumina') & (df['variant_type'] == 'INDEL')]
    # illumina_snvs = df[(df['platform'] == 'Illumina') & (df['variant_type'] == 'SNV')]
    # etc.
    
    # For demonstration, let's create synthetic data
    print("Generating synthetic data for demonstration...")
    np.random.seed(42)
    
    def generate_synthetic_data(n_true=1000, n_false=1000):
        # True variants: higher scores (mean=0.8)
        true_scores = np.random.beta(8, 2, n_true)
        # False variants: lower scores (mean=0.3)
        false_scores = np.random.beta(2, 5, n_false)
        
        y_true = np.concatenate([np.ones(n_true), np.zeros(n_false)])
        y_scores = np.concatenate([true_scores, false_scores])
        
        return y_true, y_scores
    
    results_dict = {
        'illumina_indels': generate_synthetic_data(1000, 1500),
        'illumina_snvs': generate_synthetic_data(3000, 2000),
        'ion_indels': generate_synthetic_data(800, 1200),
        'ion_snvs': generate_synthetic_data(2500, 1800)
    }
    
    # Generate figure
    stats = generate_figure_s1(results_dict, 'figure_S1_example.png')
    
    # Print summary statistics
    print("\n" + "="*60)
    print("SUMMARY STATISTICS")
    print("="*60)
    for variant_type, stat in stats.items():
        print(f"\n{variant_type.upper()}:")
        print(f"  Max Accuracy Threshold: {stat['max_acc_threshold']:.3f}")
        print(f"  TPR at max accuracy: {stat['tpr_at_max_acc']:.1%}")
        print(f"  TNR at max accuracy: {stat['tnr_at_max_acc']:.1%}")
        print(f"  TPR at 90% TNR: {stat['tpr_at_90_tnr']:.1%}")
        print(f"  TPR at 95% TNR: {stat['tpr_at_95_tnr']:.1%}")

if __name__ == '__main__':
    main()