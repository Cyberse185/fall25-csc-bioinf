import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc
from load_data import GarfieldDataLoader

def plot_rocs_for_dataset(df, features, title, filename):
    """
    Generates a 'Spaghetti Plot' ROC curve for a specific dataset.
    It loops through every feature in the list and draws a line for it.
    """
    plt.figure(figsize=(12, 10))
    
    # Create True/False labels (1 for True, 0 for False)
    # The dataset uses 'T' for True and 'F' for False in the 'Class' column
    y_true = (df['Class'] == 'T').astype(int)
    
    print(f"\n--- Generating {title} ---")
    print(f"Processing {len(features)} features...")
    
    # Generate distinct colors for the lines
    colors = plt.cm.tab20(np.linspace(0, 1, len(features)))
    
    for feature, color in zip(features, colors):
        if feature not in df.columns:
            print(f"Warning: Feature '{feature}' not found in dataframe. Skipping.")
            continue

        # Get the raw numbers for this feature
        # fillna(0) ensures the code doesn't crash on missing data
        y_score = df[feature].fillna(0)
        
        # Calculate ROC
        fpr, tpr, _ = roc_curve(y_true, y_score)
        roc_auc = auc(fpr, tpr)
        
        # --- THE CRITICAL "FLIP" LOGIC ---
        # Some features (like Strand Bias 'FS') are "Bad" when high.
        # This means high FS = False Mutation.
        # A standard ROC curve would dip below the diagonal (AUC < 0.5).
        # To match the paper's style, we flip these so they show positive predictive power.
        if roc_auc < 0.5:
            y_score = -y_score
            fpr, tpr, _ = roc_curve(y_true, y_score)
            roc_auc = auc(fpr, tpr)
            
        # Plot the line
        plt.plot(fpr, tpr, color=color, lw=2, alpha=0.8,
                 label=f'{feature} (AUC = {roc_auc:.3f})')

    # Formatting the Plot (to look like the paper)
    plt.plot([0, 1], [0, 1], 'k--', lw=2, label='Random (0.500)')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate', fontsize=12)
    plt.ylabel('True Positive Rate', fontsize=12)
    plt.title(title, fontsize=14, fontweight='bold')
    
    # Put legend outside the plot if there are too many features
    plt.legend(loc="center left", bbox_to_anchor=(1, 0.5), fontsize=10)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    
    print(f"Saving {filename}...")
    plt.savefig(filename, dpi=300)
    # plt.show() # Uncomment if you want to see them pop up
    plt.close()

def main():
    # 1. Initialize Loader
    loader = GarfieldDataLoader()
    
    # --- PLOT A: ILLUMINA INDEL ---
    # Note: We use the generic ILLUMINA_FEATURES list for both SNV and INDEL 
    # unless you have a specific INDEL list defined.
    print("Loading Illumina INDEL...")
    data = loader.load_illumina_indel()
    plot_rocs_for_dataset(
        data['train'], 
        loader.ILLUMINA_FEATURES, 
        "S16(a) Illumina INS/DELs Features", 
        "Figure_S16a_Illumina_INDEL.png"
    )

    # --- PLOT B: ION INDEL ---
    print("Loading ION INDEL...")
    data = loader.load_ion_indel()
    plot_rocs_for_dataset(
        data['train'], 
        loader.ION_INDEL_FEATURES,  # Uses the special list
        "S16(b) ION INS/DELs Features", 
        "Figure_S16b_Ion_INDEL.png"
    )

    # --- PLOT C: ILLUMINA SNV ---
    print("Loading Illumina SNV...")
    data = loader.load_illumina_snp()
    plot_rocs_for_dataset(
        data['train'], 
        loader.ILLUMINA_FEATURES, 
        "S16(c) Illumina SNVs Features", 
        "Figure_S16c_Illumina_SNV.png"
    )

    # --- PLOT D: ION SNV ---
    print("Loading ION SNV...")
    data = loader.load_ion_snp()
    plot_rocs_for_dataset(
        data['train'], 
        loader.ION_FEATURES, 
        "S16(d) ION SNVs Features", 
        "Figure_S16d_Ion_SNV.png"
    )

    print("\nAll 4 Figures Generated Successfully!")

if __name__ == "__main__":
    main()