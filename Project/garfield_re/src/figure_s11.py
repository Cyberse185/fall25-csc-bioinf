"""Generate Supplementary Figure S11: Feature Importance in GARFIELD-NGS Models

Shows scaled importance of each feature in prediction models for:
(a) Illumina INS/DELs
(b) ION INS/DELs
(c) Illumina SNVs
(d) ION SNVs
"""

import h2o
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

def get_feature_importance(model):
    """Extract and scale feature importance from H2O model"""
    # Get variable importance
    varimp = model.varimp(use_pandas=True)
    
    # Scale importance to 0-1 range
    if varimp is not None and len(varimp) > 0:
        features = varimp['variable'].values
        importance = varimp['scaled_importance'].values
        return features, importance
    else:
        return None, None

def plot_feature_importance(ax, features, importance, title):
    """Plot feature importance bar chart"""
    # Sort by importance
    sorted_idx = np.argsort(importance)[::-1]
    features_sorted = features[sorted_idx]
    importance_sorted = importance[sorted_idx]
    
    # Create bar plot
    x_pos = np.arange(len(features_sorted))
    ax.bar(x_pos, importance_sorted, color='#4169E1', width=0.7)
    
    # Formatting
    ax.set_xticks(x_pos)
    ax.set_xticklabels(features_sorted, rotation=45, ha='right', fontsize=9)
    ax.set_ylabel('variable scaled importance', fontsize=11)
    ax.set_ylim([0, 1.05])
    ax.set_title(title, fontsize=12, fontweight='bold', loc='left')
    ax.grid(axis='y', alpha=0.3)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

def create_figure_s11():
    """Generate complete Figure S11 with all 4 subplots"""
    
    print("="*70)
    print("GENERATING FIGURE S11: FEATURE IMPORTANCE")
    print("="*70)
    
    # Initialize H2O
    h2o.init(max_mem_size="4G")
    h2o.no_progress()
    
    print("\n1. Loading trained models...")
    model_dir = Path("../models")
    
    ilm_indel_model = h2o.load_model(str(model_dir / "illumina_indel_paper_exact"))
    ion_indel_model = h2o.load_model(str(model_dir / "ion_indel_paper_exact"))
    ilm_snv_model = h2o.load_model(str(model_dir / "illumina_snv_paper_exact"))
    ion_snv_model = h2o.load_model(str(model_dir / "ion_snv_paper_exact"))
    
    print("  ✓ All models loaded")
    
    print("\n2. Extracting feature importance...")
    
    # Get feature importance from each model
    ilm_indel_feat, ilm_indel_imp = get_feature_importance(ilm_indel_model)
    ion_indel_feat, ion_indel_imp = get_feature_importance(ion_indel_model)
    ilm_snv_feat, ilm_snv_imp = get_feature_importance(ilm_snv_model)
    ion_snv_feat, ion_snv_imp = get_feature_importance(ion_snv_model)
    
    print("  ✓ Feature importance extracted")
    
    print("\n3. Creating Figure S11...")
    
    # Create 2x2 subplot
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # (a) Illumina INS/DELs
    if ilm_indel_feat is not None:
        plot_feature_importance(axes[0, 0], ilm_indel_feat, ilm_indel_imp,
                               'a')
        print(f"  ✓ Illumina INS/DELs: Top feature = {ilm_indel_feat[0]}")
    
    # (b) ION INS/DELs
    if ion_indel_feat is not None:
        plot_feature_importance(axes[0, 1], ion_indel_feat, ion_indel_imp,
                               'b')
        print(f"  ✓ ION INS/DELs: Top feature = {ion_indel_feat[0]}")
    
    # (c) Illumina SNVs
    if ilm_snv_feat is not None:
        plot_feature_importance(axes[1, 0], ilm_snv_feat, ilm_snv_imp,
                               'c')
        print(f"  ✓ Illumina SNVs: Top feature = {ilm_snv_feat[0]}")
    
    # (d) ION SNVs
    if ion_snv_feat is not None:
        plot_feature_importance(axes[1, 1], ion_snv_feat, ion_snv_imp,
                               'd')
        print(f"  ✓ ION SNVs: Top feature = {ion_snv_feat[0]}")
    
    plt.tight_layout()
    
    # Save figure
    output_file = 'figure_s11_feature_importance.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\n✓ Figure saved: {output_file}")
    
    # Print detailed feature rankings
    print("\n" + "="*70)
    print("FEATURE IMPORTANCE RANKINGS")
    print("="*70)
    
    if ilm_indel_feat is not None:
        print("\nIllumina INS/DELs:")
        for i, (feat, imp) in enumerate(zip(ilm_indel_feat[:5], ilm_indel_imp[:5]), 1):
            print(f"  {i}. {feat}: {imp:.4f}")
    
    if ion_indel_feat is not None:
        print("\nION INS/DELs:")
        for i, (feat, imp) in enumerate(zip(ion_indel_feat[:5], ion_indel_imp[:5]), 1):
            print(f"  {i}. {feat}: {imp:.4f}")
    
    if ilm_snv_feat is not None:
        print("\nIllumina SNVs:")
        for i, (feat, imp) in enumerate(zip(ilm_snv_feat[:5], ilm_snv_imp[:5]), 1):
            print(f"  {i}. {feat}: {imp:.4f}")
    
    if ion_snv_feat is not None:
        print("\nION SNVs:")
        for i, (feat, imp) in enumerate(zip(ion_snv_feat[:5], ion_snv_imp[:5]), 1):
            print(f"  {i}. {feat}: {imp:.4f}")
    
    print("="*70)
    
    # Cleanup
    h2o.cluster().shutdown()
    
    return fig

if __name__ == "__main__":
    create_figure_s11()
    plt.show()