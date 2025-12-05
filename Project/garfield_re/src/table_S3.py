"""Reproduce Supplementary Table S3 - GARFIELD-NGS Performance at Different Thresholds"""

import h2o
import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.metrics import (accuracy_score, recall_score, precision_score, 
                             matthews_corrcoef, confusion_matrix, roc_curve, auc)
from load_data import GarfieldDataLoader

def calculate_metrics(y_true, y_pred, y_proba):
    """Calculate all metrics for Table S3"""
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()
    
    acc = accuracy_score(y_true, y_pred)
    tpr = tp / (tp + fn) if (tp + fn) > 0 else 0  # Sensitivity
    fdr = fp / (fp + tp) if (fp + tp) > 0 else 0  # False Discovery Rate
    tnr = tn / (tn + fp) if (tn + fp) > 0 else 0  # Specificity
    ppv = tp / (tp + fp) if (tp + fp) > 0 else 0  # Precision
    npv = tn / (tn + fn) if (tn + fn) > 0 else 0
    mcc = matthews_corrcoef(y_true, y_pred)
    
    return {
        'ACC': acc,
        'TPR': tpr,
        'FDR': fdr,
        'TNR': tnr,
        'PPV': ppv,
        'NPV': npv,
        'MCC': mcc
    }

def find_threshold_for_tpr(y_true, y_proba, target_tpr):
    """Find threshold that achieves target TPR"""
    fpr, tpr, thresholds = roc_curve(y_true, y_proba)
    
    # Find threshold closest to target TPR
    idx = np.argmin(np.abs(tpr - target_tpr))
    return thresholds[idx]

def find_max_accuracy_threshold(y_true, y_proba):
    """Find threshold that maximizes accuracy"""
    fpr, tpr, thresholds = roc_curve(y_true, y_proba)
    
    best_acc = 0
    best_threshold = 0.5
    
    for threshold in thresholds:
        y_pred = (y_proba >= threshold).astype(int)
        acc = accuracy_score(y_true, y_pred)
        if acc > best_acc:
            best_acc = acc
            best_threshold = threshold
    
    return best_threshold

def find_max_f1_threshold(y_true, y_proba):
    """Find threshold that maximizes F1 score"""
    fpr, tpr, thresholds = roc_curve(y_true, y_proba)
    
    best_f1 = 0
    best_threshold = 0.5
    
    for threshold in thresholds:
        y_pred = (y_proba >= threshold).astype(int)
        precision = precision_score(y_true, y_pred, zero_division=0)
        recall = recall_score(y_true, y_pred, zero_division=0)
        f1 = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0
        
        if f1 > best_f1:
            best_f1 = f1
            best_threshold = threshold
    
    return best_threshold, best_f1

def evaluate_model(model_path, test_data, features, variant_type):
    """Evaluate a single model on test set"""
    print(f"\n{'='*70}")
    print(f"Evaluating: {variant_type}")
    print(f"{'='*70}")
    
    # Load model
    print(f"Loading model from: {model_path}")
    model = h2o.load_model(model_path)
    
    # Prepare test data
    X_test = test_data[features].copy()
    y_test = (test_data['Class'] == 'T').astype(int)
    
    # Handle missing values
    X_test = X_test.fillna(X_test.median())
    X_test = X_test.replace([np.inf, -np.inf], np.nan)
    X_test = X_test.fillna(X_test.median())
    
    # Convert to H2O
    test_h2o_data = X_test.copy()
    test_h2o_data['target'] = y_test.values
    test_h2o = h2o.H2OFrame(test_h2o_data)
    test_h2o['target'] = test_h2o['target'].asfactor()
    
    # Get predictions
    predictions = model.predict(test_h2o)
    y_proba = predictions['p1'].as_data_frame().values.flatten()
    y_true = y_test.values
    
    print(f"Test set size: {len(y_true)}")
    print(f"Class balance: {pd.Series(y_true).value_counts().to_dict()}")
    
    # Calculate AUC
    test_auc = auc(*roc_curve(y_true, y_proba)[:2])
    print(f"Test AUC: {test_auc:.4f}")
    
    # Find thresholds
    print("\nFinding optimal thresholds...")
    threshold_tpr_99 = find_threshold_for_tpr(y_true, y_proba, 0.99)
    threshold_tpr_95 = find_threshold_for_tpr(y_true, y_proba, 0.95)
    threshold_max_acc = find_max_accuracy_threshold(y_true, y_proba)
    threshold_max_f1, max_f1_score = find_max_f1_threshold(y_true, y_proba)
    
    # Calculate metrics at each threshold
    results = []
    
    for criterion, threshold in [
        ('TPR 0.99', threshold_tpr_99),
        ('TPR 0.95', threshold_tpr_95),
        ('Maximum accuracy', threshold_max_acc),
        (f'Maximum f1 ({max_f1_score:.3f})', threshold_max_f1)
    ]:
        y_pred = (y_proba >= threshold).astype(int)
        metrics = calculate_metrics(y_true, y_pred, y_proba)
        metrics['criterion'] = criterion
        metrics['threshold'] = threshold
        results.append(metrics)
    
    # Create results DataFrame
    df = pd.DataFrame(results)
    df = df[['criterion', 'ACC', 'TPR', 'FDR', 'TNR', 'PPV', 'NPV', 'MCC', 'threshold']]
    
    print("\n" + "="*70)
    print(f"Results for {variant_type}:")
    print("="*70)
    print(df.to_string(index=False))
    
    return df, test_auc

def main():
    print("="*70)
    print("REPRODUCING SUPPLEMENTARY TABLE S3")
    print("GARFIELD-NGS Performance at Different Thresholds")
    print("="*70)
    
    # Initialize H2O
    h2o.init(max_mem_size="4G")
    h2o.no_progress()
    
    # Load data
    loader = GarfieldDataLoader()
    
    # Model configurations
    models = [
        {
            'name': 'Illumina INS/DELs',
            'path': '/Users/tomislavstojcic/fall25-csc-bioinf/Project/models/illumina_indel_paper_exact',
            'test_file': 'ILM_INDEL_Test.txt',
            'features': loader.ILLUMINA_FEATURES
        },
        {
            'name': 'Illumina SNVs',
            'path': '/Users/tomislavstojcic/fall25-csc-bioinf/Project/models/illumina_snv_paper_exact',
            'test_file': 'ILM_SNP_Test.txt',
            'features': loader.ILLUMINA_FEATURES
        },
        {
            'name': 'ION INS/DELs',
            'path': '/Users/tomislavstojcic/fall25-csc-bioinf/Project/models/ion_indel_paper_exact',
            'test_file': 'ION_INDEL_Test.txt',
            'features': loader.ION_INDEL_FEATURES
        },
        {
            'name': 'ION SNVs',
            'path': '/Users/tomislavstojcic/fall25-csc-bioinf/Project/models/ion_snv_paper_exact',
            'test_file': 'ION_SNP_Test.txt',
            'features': loader.ION_FEATURES
        }
    ]
    
    # Evaluate each model
    all_results = {}
    
    for model_config in models:
        # Load test data
        test_data = loader.load_dataset(model_config['test_file'])
        
        # Evaluate
        results_df, test_auc = evaluate_model(
            model_config['path'],
            test_data,
            model_config['features'],
            model_config['name']
        )
        
        all_results[model_config['name']] = {
            'results': results_df,
            'auc': test_auc
        }
    
    # Save combined results
    print("\n" + "="*70)
    print("SAVING RESULTS")
    print("="*70)
    
    output_dir = Path("results")
    output_dir.mkdir(exist_ok=True)
    
    # Save individual results
    for variant_type, data in all_results.items():
        safe_name = variant_type.replace('/', '_').replace(' ', '_')
        output_file = output_dir / f"table_s3_{safe_name}.csv"
        data['results'].to_csv(output_file, index=False)
        print(f"Saved: {output_file}")
    
    # Create combined table
    print("\n" + "="*70)
    print("FINAL TABLE S3 - GARFIELD-NGS PERFORMANCE")
    print("="*70)
    
    for variant_type, data in all_results.items():
        print(f"\n{variant_type} (AUC: {data['auc']:.4f})")
        print("-" * 70)
        print(data['results'].to_string(index=False))
    
    # Save combined results
    combined_file = output_dir / "table_s3_combined.txt"
    with open(combined_file, 'w') as f:
        f.write("SUPPLEMENTARY TABLE S3 - GARFIELD-NGS PERFORMANCE\n")
        f.write("="*70 + "\n\n")
        for variant_type, data in all_results.items():
            f.write(f"\n{variant_type} (AUC: {data['auc']:.4f})\n")
            f.write("-" * 70 + "\n")
            f.write(data['results'].to_string(index=False))
            f.write("\n\n")
    
    print(f"\nCombined results saved to: {combined_file}")
    
    h2o.cluster().shutdown()
    print("\n✓ Table S3 reproduction complete!")

if __name__ == "__main__":
    main()