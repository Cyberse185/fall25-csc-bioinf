"""Generate Supplementary Figure S7: ROC Curves on Training and Validation Sets

ROC curves of final GARFIELD-NGS models on training and validation datasets.
Performance evaluated on:
- Illumina (a) training set and (c) validation set
- ION (b) training set and (d) validation set
Orange lines = INS/DEL variants, Blue lines = SNV variants
"""

import h2o
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc
from pathlib import Path
from load_data import GarfieldDataLoader

def get_predictions(model, X, y, features):
    """Get model predictions and true labels"""
    # Prepare H2O frame
    data = X.copy()
    data['target'] = y.values
    h2o_frame = h2o.H2OFrame(data)
    h2o_frame['target'] = h2o_frame['target'].asfactor()
    
    # Get predictions
    preds = model.predict(h2o_frame)
    pred_df = preds.as_data_frame()
    
    # Extract probability of class 1 (true variant)
    y_score = pred_df['p1'].values
    y_true = y.values
    
    return y_true, y_score

def plot_roc_curve(ax, y_true, y_score, label, color):
    """Plot ROC curve and return AUC"""
    fpr, tpr, _ = roc_curve(y_true, y_score)
    roc_auc = auc(fpr, tpr)
    
    ax.plot(fpr, tpr, color=color, lw=2, 
            label=f'{label} (AUC = {roc_auc:.4f})')
    
    return roc_auc

def create_figure_s7():
    """Generate complete Figure S7 with all 4 subplots"""
    
    print("="*70)
    print("GENERATING FIGURE S7: ROC CURVES")
    print("="*70)
    
    # Initialize H2O
    h2o.init(max_mem_size="4G")
    h2o.no_progress()
    
    # Load data
    loader = GarfieldDataLoader()
    
    print("\n1. Loading datasets...")
    ilm_snv_data = loader.load_illumina_snp()
    ilm_indel_data = loader.load_illumina_indel()
    ion_snv_data = loader.load_ion_snp()
    ion_indel_data = loader.load_ion_indel()
    
    print("\n2. Loading trained models...")
    model_dir = Path("../models")
    ilm_snv_model = h2o.load_model(str(model_dir / "illumina_snv_paper_exact"))
    ilm_indel_model = h2o.load_model(str(model_dir / "illumina_indel_paper_exact"))
    ion_snv_model = h2o.load_model(str(model_dir / "ion_snv_paper_exact"))
    ion_indel_model = h2o.load_model(str(model_dir / "ion_indel_paper_exact"))
    print("  ✓ All models loaded")
    
    print("\n3. Preprocessing data...")
    # Illumina SNV
    ilm_snv_train_X, ilm_snv_train_y = loader.preprocess_data(
        ilm_snv_data['train'], loader.ILLUMINA_FEATURES)
    ilm_snv_val_X, ilm_snv_val_y = loader.preprocess_data(
        ilm_snv_data['validation'], loader.ILLUMINA_FEATURES)
    
    # Illumina INDEL
    ilm_indel_train_X, ilm_indel_train_y = loader.preprocess_data(
        ilm_indel_data['train'], loader.ILLUMINA_FEATURES)
    ilm_indel_val_X, ilm_indel_val_y = loader.preprocess_data(
        ilm_indel_data['validation'], loader.ILLUMINA_FEATURES)
    
    # ION SNV
    ion_snv_train_X, ion_snv_train_y = loader.preprocess_data(
        ion_snv_data['train'], loader.ION_FEATURES)
    ion_snv_val_X, ion_snv_val_y = loader.preprocess_data(
        ion_snv_data['validation'], loader.ION_FEATURES)
    
    # ION INDEL
    ion_indel_train_X, ion_indel_train_y = loader.preprocess_data(
        ion_indel_data['train'], loader.ION_INDEL_FEATURES)
    ion_indel_val_X, ion_indel_val_y = loader.preprocess_data(
        ion_indel_data['validation'], loader.ION_INDEL_FEATURES)
    
    print("\n4. Generating predictions...")
    
    # Illumina Training
    print("  - Illumina training set...")
    ilm_snv_train_true, ilm_snv_train_score = get_predictions(
        ilm_snv_model, ilm_snv_train_X, ilm_snv_train_y, loader.ILLUMINA_FEATURES)
    ilm_indel_train_true, ilm_indel_train_score = get_predictions(
        ilm_indel_model, ilm_indel_train_X, ilm_indel_train_y, loader.ILLUMINA_FEATURES)
    
    # Illumina Validation
    print("  - Illumina validation set...")
    ilm_snv_val_true, ilm_snv_val_score = get_predictions(
        ilm_snv_model, ilm_snv_val_X, ilm_snv_val_y, loader.ILLUMINA_FEATURES)
    ilm_indel_val_true, ilm_indel_val_score = get_predictions(
        ilm_indel_model, ilm_indel_val_X, ilm_indel_val_y, loader.ILLUMINA_FEATURES)
    
    # ION Training
    print("  - ION training set...")
    ion_snv_train_true, ion_snv_train_score = get_predictions(
        ion_snv_model, ion_snv_train_X, ion_snv_train_y, loader.ION_FEATURES)
    ion_indel_train_true, ion_indel_train_score = get_predictions(
        ion_indel_model, ion_indel_train_X, ion_indel_train_y, loader.ION_INDEL_FEATURES)
    
    # ION Validation
    print("  - ION validation set...")
    ion_snv_val_true, ion_snv_val_score = get_predictions(
        ion_snv_model, ion_snv_val_X, ion_snv_val_y, loader.ION_FEATURES)
    ion_indel_val_true, ion_indel_val_score = get_predictions(
        ion_indel_model, ion_indel_val_X, ion_indel_val_y, loader.ION_INDEL_FEATURES)
    
    print("\n5. Creating Figure S7...")
    
    # Create 2x2 subplot
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Color scheme from paper
    color_indel = '#FF8C00'  # Orange for INS/DEL
    color_snv = '#4169E1'    # Blue for SNV
    
    # (a) Illumina Training Set
    ax = axes[0, 0]
    auc_ilm_indel_train = plot_roc_curve(ax, ilm_indel_train_true, ilm_indel_train_score, 
                                          'INS/DELs', color_indel)
    auc_ilm_snv_train = plot_roc_curve(ax, ilm_snv_train_true, ilm_snv_train_score, 
                                        'SNVs', color_snv)
    ax.plot([0, 1], [0, 1], 'k--', lw=1, alpha=0.3)
    ax.set_xlabel('False positive rate', fontsize=11)
    ax.set_ylabel('True positive rate', fontsize=11)
    ax.set_title('(a) Illumina Training Set', fontsize=12, fontweight='bold')
    ax.legend(loc='lower right', fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_xlim([0, 1])
    ax.set_ylim([0, 1])
    
    # (b) ION Training Set
    ax = axes[0, 1]
    auc_ion_indel_train = plot_roc_curve(ax, ion_indel_train_true, ion_indel_train_score, 
                                          'INS/DELs', color_indel)
    auc_ion_snv_train = plot_roc_curve(ax, ion_snv_train_true, ion_snv_train_score, 
                                        'SNVs', color_snv)
    ax.plot([0, 1], [0, 1], 'k--', lw=1, alpha=0.3)
    ax.set_xlabel('False positive rate', fontsize=11)
    ax.set_ylabel('True positive rate', fontsize=11)
    ax.set_title('(b) ION Training Set', fontsize=12, fontweight='bold')
    ax.legend(loc='lower right', fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_xlim([0, 1])
    ax.set_ylim([0, 1])
    
    # (c) Illumina Validation Set
    ax = axes[1, 0]
    auc_ilm_indel_val = plot_roc_curve(ax, ilm_indel_val_true, ilm_indel_val_score, 
                                        'INS/DELs', color_indel)
    auc_ilm_snv_val = plot_roc_curve(ax, ilm_snv_val_true, ilm_snv_val_score, 
                                      'SNVs', color_snv)
    ax.plot([0, 1], [0, 1], 'k--', lw=1, alpha=0.3)
    ax.set_xlabel('False positive rate', fontsize=11)
    ax.set_ylabel('True positive rate', fontsize=11)
    ax.set_title('(c) Illumina Validation Set', fontsize=12, fontweight='bold')
    ax.legend(loc='lower right', fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_xlim([0, 1])
    ax.set_ylim([0, 1])
    
    # (d) ION Validation Set
    ax = axes[1, 1]
    auc_ion_indel_val = plot_roc_curve(ax, ion_indel_val_true, ion_indel_val_score, 
                                        'INS/DELs', color_indel)
    auc_ion_snv_val = plot_roc_curve(ax, ion_snv_val_true, ion_snv_val_score, 
                                      'SNVs', color_snv)
    ax.plot([0, 1], [0, 1], 'k--', lw=1, alpha=0.3)
    ax.set_xlabel('False positive rate', fontsize=11)
    ax.set_ylabel('True positive rate', fontsize=11)
    ax.set_title('(d) ION Validation Set', fontsize=12, fontweight='bold')
    ax.legend(loc='lower right', fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_xlim([0, 1])
    ax.set_ylim([0, 1])
    
    plt.tight_layout()
    
    # Save figure
    output_file = 'figure_s7_roc_curves.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\n✓ Figure saved: {output_file}")
    
    # Print summary
    print("\n" + "="*70)
    print("AUROC SUMMARY")
    print("="*70)
    print(f"\nIllumina Training Set:")
    print(f"  INS/DELs: {auc_ilm_indel_train:.4f} (Paper: 0.9367)")
    print(f"  SNVs:     {auc_ilm_snv_train:.4f} (Paper: 0.7981)")
    
    print(f"\nION Training Set:")
    print(f"  INS/DELs: {auc_ion_indel_train:.4f} (Paper: 0.9708)")
    print(f"  SNVs:     {auc_ion_snv_train:.4f} (Paper: 0.9667)")
    
    print(f"\nIllumina Validation Set:")
    print(f"  INS/DELs: {auc_ilm_indel_val:.4f}")
    print(f"  SNVs:     {auc_ilm_snv_val:.4f}")
    
    print(f"\nION Validation Set:")
    print(f"  INS/DELs: {auc_ion_indel_val:.4f}")
    print(f"  SNVs:     {auc_ion_snv_val:.4f}")
    print("="*70)
    
    # Cleanup
    h2o.cluster().shutdown()
    
    return fig

if __name__ == "__main__":
    create_figure_s7()
    plt.show()