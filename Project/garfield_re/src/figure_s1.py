#!/usr/bin/env python3
"""
Generate EXACT Figure S1 from the paper - ION INDEL (Panel c)
TPR/TNR curves vs GARFIELD-NGS score threshold
"""

import h2o
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from load_data import GarfieldDataLoader

# --- CONFIGURATION ---
MODEL_PATH = "../models/ion_indel_paper_exact" 
FIGURE_OUTPUT = "Figure_S1_IonIndel_PaperStyle.png"

def calculate_tpr_tnr_curves(y_true, y_scores):
    # (Same function as before, kept for brevity)
    thresholds = np.linspace(0, 1, 1000)
    tpr_values = []
    tnr_values = []
    
    for threshold in thresholds:
        predictions = (y_scores >= threshold).astype(int)
        tp = np.sum((predictions == 1) & (y_true == 1))
        fn = np.sum((predictions == 0) & (y_true == 1))
        tn = np.sum((predictions == 0) & (y_true == 0))
        fp = np.sum((predictions == 1) & (y_true == 0))
        
        tpr = tp / (tp + fn) if (tp + fn) > 0 else 0
        tnr = tn / (tn + fp) if (tn + fp) > 0 else 0
        
        tpr_values.append(tpr)
        tnr_values.append(tnr)
    
    return thresholds, np.array(tpr_values), np.array(tnr_values)

def find_max_accuracy_threshold(y_true, y_scores):
    # (Same function as before)
    thresholds = np.linspace(0, 1, 1000)
    best_accuracy = 0
    best_threshold = 0.5
    for threshold in thresholds:
        predictions = (y_scores >= threshold).astype(int)
        accuracy = np.mean(predictions == y_true)
        if accuracy > best_accuracy:
            best_accuracy = accuracy
            best_threshold = threshold
    return best_threshold, best_accuracy

def find_tpr_at_tnr(thresholds, tpr_values, tnr_values, target_tnr):
    idx = np.argmin(np.abs(tnr_values - target_tnr))
    return tpr_values[idx]

def plot_figure_s1_paper_style():
    """Generate paper-style Figure S1 for ION INDEL"""
    
    # 1. Initialize H2O and Load Model
    h2o.init(max_mem_size="4G")
    h2o.no_progress()
    
    print(f"Loading model from: {MODEL_PATH}")
    try:
        model = h2o.load_model(MODEL_PATH)
    except:
        print(f"ERROR: Could not load model at {MODEL_PATH}")
        return

    # 2. Load Test Data - CORRECTED FOR INDEL
    print("Loading Test Data (ION INDEL)...")
    loader = GarfieldDataLoader()
    datasets = loader.load_ion_indel() # <--- CHANGED TO INDEL
    
    # --- THE FIX: Manually define features without PB/PBP ---
    # Ion Indels typically have these features but NOT PB/PBP
    INDEL_FEATURES = [f for f in loader.ION_FEATURES if f not in ['PB', 'PBP']]
    
    X_test, y_test = loader.preprocess_data(datasets['test'], INDEL_FEATURES)
    
    # Convert to H2OFrame
    test_h2o = h2o.H2OFrame(X_test)
    
    # 3. Get Predictions
    print("Running predictions...")
    preds = model.predict(test_h2o)
    pred_df = preds.as_data_frame()
    scores = pred_df['p1'].values
    y_true = y_test.values
    
    # 4. Calculate TPR/TNR curves
    print("Calculating TPR/TNR curves...")
    thresholds, tpr_values, tnr_values = calculate_tpr_tnr_curves(y_true, scores)
    
    # 5. Find important points
    max_acc_threshold, max_acc = find_max_accuracy_threshold(y_true, scores)
    tpr_at_90_tnr = find_tpr_at_tnr(thresholds, tpr_values, tnr_values, 0.90)
    tpr_at_95_tnr = find_tpr_at_tnr(thresholds, tpr_values, tnr_values, 0.95)
    
    idx_max = np.argmin(np.abs(thresholds - max_acc_threshold))
    tpr_at_max = tpr_values[idx_max]
    tnr_at_max = tnr_values[idx_max]
    
    print(f"\nResults (Ion Indel):")
    print(f"  Max Accuracy Threshold: {max_acc_threshold:.4f}")
    print(f"  TPR at max accuracy: {tpr_at_max:.4f}")
    print(f"  TNR at max accuracy: {tnr_at_max:.4f}")
    
    # 6. Create the plot (paper style)
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Plot TPR curve (red)
    ax.plot(thresholds, tpr_values, 'r-', linewidth=2, label='True Variants (TPR)')
    # Plot TNR curve (green)
    ax.plot(thresholds, tnr_values, 'g-', linewidth=2, label='False Variants (TNR)')
    
    # Add vertical dashed line for max accuracy threshold
    ax.axvline(x=max_acc_threshold, color='black', linestyle='--', linewidth=2,
               label=f'Max accuracy threshold')
    
    # Annotate TPR/TNR
    ax.annotate(f'{tpr_at_max:.4f}', xy=(max_acc_threshold, tpr_at_max),
                xytext=(max_acc_threshold - 0.15, tpr_at_max + 0.05),
                fontsize=12, color='red', weight='bold',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='red'))
    
    ax.annotate(f'{tnr_at_max:.4f}', xy=(max_acc_threshold, tnr_at_max),
                xytext=(max_acc_threshold - 0.15, tnr_at_max - 0.1),
                fontsize=12, color='green', weight='bold',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='green'))
    
    # Add text box with stats
    stats_text = (f"       TPR\nTNR 90  {tpr_at_90_tnr:.4f}\nTNR 95  {tpr_at_95_tnr:.4f}")
    ax.text(0.75, 0.15, stats_text, transform=ax.transAxes, fontsize=11,
            verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    # Labels and styling
    ax.set_xlabel('Score', fontsize=14, weight='bold')
    ax.set_ylabel('TPR/TNR', fontsize=14, weight='bold')
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.05)
    ax.legend(loc='upper left', fontsize=11)
    ax.grid(True, alpha=0.3)
    
    # Title "c" for Indels
    ax.set_title('c', fontsize=18, weight='bold', loc='left')
    
    plt.tight_layout()
    plt.savefig(FIGURE_OUTPUT, dpi=300, bbox_inches='tight')
    print(f"\n✓ Saved to {FIGURE_OUTPUT}")
    
    h2o.cluster().shutdown()

if __name__ == "__main__":
    plot_figure_s1_paper_style()