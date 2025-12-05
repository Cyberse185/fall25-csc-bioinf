"""Train Illumina INDEL - Exact Paper Implementation"""

import h2o
from h2o.estimators import H2ODeepLearningEstimator, H2OAutoEncoderEstimator
import pandas as pd
import numpy as np
import random
from pathlib import Path
from load_data import GarfieldDataLoader
from model_configs import get_model_config

def balance_training_data(X, y, target_false_ratio=0.20):
    """
    Balance training data by randomly removing true calls to ensure
    at least 20% false variants (as per paper methodology).
    
    Assumes: 1 = True variant, 0 = False variant
    """
    print("\n   Checking class balance...")
    counts = y.value_counts()
    print(f"   Original distribution: {counts.to_dict()}")
    
    # Calculate current false variant ratio
    n_false = (y == 0).sum()
    n_true = (y == 1).sum()
    current_false_ratio = n_false / len(y)
    
    print(f"   Current false variant ratio: {current_false_ratio:.2%}")
    
    if current_false_ratio >= target_false_ratio:
        print(f"   ✓ Already meets paper requirement (>={target_false_ratio:.0%} false variants)")
        return X, y
    
    # Need to downsample true variants
    # Formula: n_false / (n_false + n_true_keep) = target_ratio
    # Solving: n_true_keep = n_false * (1/target_ratio - 1)
    n_true_keep = int(n_false * (1/target_false_ratio - 1))
    n_true_remove = n_true - n_true_keep
    
    print(f"   Downsampling true variants: {n_true} → {n_true_keep} (removing {n_true_remove})")
    
    # Get indices
    false_idx = y[y == 0].index
    true_idx = y[y == 1].index
    
    # Randomly sample true variants to keep
    np.random.seed(42)
    true_keep_idx = np.random.choice(true_idx, n_true_keep, replace=False)
    
    # Combine indices
    keep_idx = np.concatenate([false_idx, true_keep_idx])
    np.random.shuffle(keep_idx)
    
    # Create balanced dataset
    X_balanced = X.loc[keep_idx].reset_index(drop=True)
    y_balanced = y.loc[keep_idx].reset_index(drop=True)
    
    # Verify
    final_false_ratio = (y_balanced == 0).sum() / len(y_balanced)
    print(f"   Final distribution: {y_balanced.value_counts().to_dict()}")
    print(f"   Final false variant ratio: {final_false_ratio:.2%}")
    print(f"   ✓ Balanced to paper specifications")
    
    return X_balanced, y_balanced

def train_illumina_indel_paper_exact():
    # Set seeds for reproducibility
    random.seed(42)
    np.random.seed(42)
    
    print("="*70)
    print("TRAINING: ILLUMINA INDEL - EXACT PAPER IMPLEMENTATION")
    print("="*70)
    print("Model: Rectifier activation, [10, 70, 60, 10, 10] architecture")
    print("="*70)
    
    h2o.init(max_mem_size="4G")
    h2o.no_progress()
    
    loader = GarfieldDataLoader()
    datasets = loader.load_illumina_indel()
    
    print("\n1. Preprocessing...")
    X_pretrain, y_pretrain = loader.preprocess_data(datasets['pretrain'], loader.ILLUMINA_FEATURES)
    X_train, y_train = loader.preprocess_data(datasets['train'], loader.ILLUMINA_FEATURES)
    X_val, y_val = loader.preprocess_data(datasets['validation'], loader.ILLUMINA_FEATURES)
    X_test, y_test = loader.preprocess_data(datasets['test'], loader.ILLUMINA_FEATURES)
    
    # Paper: "pre-training and training sets were balanced by randomly removing 
    # true calls so that they contain at least 20% of false variants"
    print("\n2. Balancing datasets (Paper Method: Manual Downsampling)...")
    X_pretrain, y_pretrain = balance_training_data(X_pretrain, y_pretrain, target_false_ratio=0.20)
    X_train, y_train = balance_training_data(X_train, y_train, target_false_ratio=0.20)
    
    # Validation and test sets are NOT balanced (as per paper)
    print(f"\n   Validation set: {len(X_val)} samples (NOT balanced)")
    print(f"   Test set: {len(X_test)} samples (NOT balanced)")
    
    def prepare_h2o_frame(X, y, name):
        data = X.copy()
        data['target'] = y.values
        h2o_frame = h2o.H2OFrame(data, destination_frame=name)
        h2o_frame['target'] = h2o_frame['target'].asfactor()
        return h2o_frame
    
    print("\n3. Converting to H2O...")
    pretrain_h2o = prepare_h2o_frame(X_pretrain, y_pretrain, "pretrain")
    train_h2o = prepare_h2o_frame(X_train, y_train, "train")
    val_h2o = prepare_h2o_frame(X_val, y_val, "validation")
    test_h2o = prepare_h2o_frame(X_test, y_test, "test")
    
    # Get Illumina INDEL hyperparameters from Supplementary Table S2
    config = get_model_config('illumina_indel')
    
    print("\n4. Autoencoder pre-training...")
    print(f"   Architecture: {config['hidden']}")
    print(f"   Activation: {config['activation']}")
    autoencoder = H2OAutoEncoderEstimator(
        hidden=config['hidden'],  # [10, 70, 60, 10, 10]
        activation=config['activation'],  # "Rectifier"
        l1=config['l1'],  # 0.00015
        l2=config['l2'],  # 0.01311
        rho=config['rho'],  # 0.961
        epsilon=config['epsilon'],  # 1.00E-009
        epochs=1000,  # Paper: 1,000 epochs
        stopping_rounds=5,  # Paper: 5 stopping rounds
        stopping_tolerance=1e-3,  # Paper: 10E-03
        stopping_metric="MSE",
        seed=42,
        reproducible = True
    )
    
    autoencoder.train(x=loader.ILLUMINA_FEATURES, training_frame=pretrain_h2o)
    print("  ✓ Pre-training complete!")
    
    print("\n5. Supervised training...")
    model = H2ODeepLearningEstimator(
        hidden=config['hidden'],  # [10, 70, 60, 10, 10]
        activation=config['activation'],  # "Rectifier"
        l1=config['l1'],  # 0.00015
        l2=config['l2'],  # 0.01311
        rho=config['rho'],  # 0.961
        epsilon=1e-9,       
        epochs=1000,  # Paper: 1,000 epochs
        mini_batch_size=1,  # Paper: mini-batch size = 1
        stopping_rounds=5,  # Paper: 5 stopping rounds
        stopping_tolerance=1e-3,  # Paper: 10E-03
        stopping_metric='logloss',  # Paper: log-loss
        seed=42,
        score_each_iteration=False,
        variable_importances=True,
        balance_classes=False,  # Paper: Manual downsampling, NOT H2O balancing
        pretrained_autoencoder=autoencoder,
        reproducible = True
    )
    
    print("  Training with paper parameters...")
    model.train(x=loader.ILLUMINA_FEATURES, y='target',
               training_frame=train_h2o,
               validation_frame=val_h2o)
    
    print(f"  ✓ Training AUC: {model.auc(train=True):.4f}")
    print(f"  ✓ Validation AUC: {model.auc(valid=True):.4f}")
    
    print("\n6. Evaluation...")
    print("--- DIAGNOSTIC START ---")
    preds = model.predict(test_h2o)
    print(preds.head(10))
    perf = model.model_performance(test_h2o)
    test_auc = perf.auc()
    
    print(f"\n  Test AUC: {test_auc:.4f}")
    print(f"  Paper AUC: 0.9269")  # From Supplementary Table S10
    print(f"  Difference: {abs(test_auc - 0.9269):.4f} ({abs(test_auc - 0.9269)/0.9269*100:.1f}%)")
    
    model_dir = Path("../models")
    model_dir.mkdir(exist_ok=True)
    model_path = h2o.save_model(model, path=str(model_dir), 
                                filename="illumina_indel_paper_exact", force=True)
    print(f"\n7. Saved: {model_path}")
    
    print("\n" + "="*70)
    print(f"ILLUMINA INDEL (Paper Implementation): {test_auc:.4f} (target: 0.9269)")
    print("="*70)
    
    h2o.cluster().shutdown()
    return test_auc

if __name__ == "__main__":
    train_illumina_indel_paper_exact()