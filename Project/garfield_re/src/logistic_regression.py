import h2o
from h2o.estimators import H2OGeneralizedLinearEstimator
import pandas as pd
from load_data import GarfieldDataLoader
from pathlib import Path


def train_logistic_regression():
    """Train Logistic Regression models to match paper's baseline"""
    
    print("="*70)
    print("GARFIELD-NGS: LOGISTIC REGRESSION REPLICATION")
    print("="*70)
    
    h2o.init(max_mem_size="4G")
    h2o.no_progress()
    
    loader = GarfieldDataLoader()
    results = []
    
    # Config with Paper's LR scores
    models_config = [
        {
            'name': 'Illumina SNV',
            'loader_func': loader.load_illumina_snp,
            'features': loader.ILLUMINA_FEATURES,
            'paper_lr': 0.8092
        },
        {
            'name': 'Illumina INDEL',
            'loader_func': loader.load_illumina_indel,
            'features': loader.ILLUMINA_FEATURES,
            'paper_lr': 0.9293
        },
        {
            'name': 'Ion SNV',
            'loader_func': loader.load_ion_snp,
            'features': loader.ION_FEATURES,
            'paper_lr': 0.9768
        },
        {
            'name': 'Ion INDEL',
            'loader_func': loader.load_ion_indel,
            'features': loader.ION_INDEL_FEATURES,
            'paper_lr': 0.9246
        }
    ]
    
    for config in models_config:
        print(f"\nProcessing: {config['name']}")
        print("-" * 30)
        
        # 1. Load Data
        datasets = config['loader_func']()
        X_train, y_train = loader.preprocess_data(datasets['train'], config['features'])
        X_val, y_val = loader.preprocess_data(datasets['validation'], config['features'])
        X_test, y_test = loader.preprocess_data(datasets['test'], config['features'])
        
        # 2. Prepare H2O Frames
        def make_frame(X, y, name):
            data = X.copy()
            data['target'] = y.values
            frame = h2o.H2OFrame(data, destination_frame=name)
            frame['target'] = frame['target'].asfactor()
            return frame
        
        safe_name = config['name'].replace(' ', '_')
        train_h2o = make_frame(X_train, y_train, f"{safe_name}_train")
        val_h2o = make_frame(X_val, y_val, f"{safe_name}_val")
        test_h2o = make_frame(X_test, y_test, f"{safe_name}_test")
        
        # 3. Train Logistic Regression
        print("  > Training H2O Logistic Regression...")
        
        # 'family="binomial"' is the standard setting for binary classification (Logistic Regression)
        # We use seed=42 for reproducibility
        lr_model = H2OGeneralizedLinearEstimator(
            family='binomial',
            seed=42
        )
        
        lr_model.train(
            x=config['features'], 
            y='target',
            training_frame=train_h2o,
            validation_frame=val_h2o
        )
        
        # 4. Evaluate
        perf = lr_model.model_performance(test_h2o)
        lr_auc = perf.auc()
        
        print(f"  > Done.")
        print(f"  > Our LR AUC:   {lr_auc:.4f}")
        print(f"  > Paper LR AUC: {config['paper_lr']:.4f}")
        
        print(f"  > Done.")
        print(f"  > Our LR AUC:   {lr_auc:.4f}")
        print(f"  > Paper LR AUC: {config['paper_lr']:.4f}")
        
        # SAVE MODEL (ADD THIS)
        model_dir = Path("../models")
        model_dir.mkdir(exist_ok=True)
        lr_path = h2o.save_model(lr_model, path=str(model_dir), 
                                 filename=f"lr_{safe_name.lower()}", force=True)
        print(f"  > Saved model: lr_{safe_name.lower()}")
                
        results.append({
            'Dataset': config['name'],
            'Our_LR_AUC': lr_auc,
            'Paper_LR_AUC': config['paper_lr'],
            'Difference': lr_auc - config['paper_lr']
        })

    # Summary Output
    print("\n" + "="*70)
    print("FINAL SUMMARY: LOGISTIC REGRESSION REPLICATION")
    print("="*70)
    
    df = pd.DataFrame(results)
    print(df.to_string(index=False, float_format=lambda x: "{:.4f}".format(x)))
    
    # Save to CURRENT directory
    df.to_csv('logistic_regression_replication.csv', index=False)
    print(f"\nSaved results to logistic_regression_replication.csv")

    h2o.cluster().shutdown()

if __name__ == "__main__":
    train_logistic_regression()