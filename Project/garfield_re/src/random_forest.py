import h2o
from h2o.estimators import H2ORandomForestEstimator
import pandas as pd
from load_data import GarfieldDataLoader
from pathlib import Path


def train_random_forests():
    """Train Random Forest models to match paper's baseline"""
    
    print("="*70)
    print("GARFIELD-NGS: RANDOM FOREST REPLICATION")
    print("="*70)
    
    h2o.init(max_mem_size="4G")
    h2o.no_progress()
    
    loader = GarfieldDataLoader()
    results = []
    
    # Config filtered to just necessary RF data
    models_config = [
        {
            'name': 'Illumina SNV',
            'loader_func': loader.load_illumina_snp,
            'features': loader.ILLUMINA_FEATURES,
            'paper_rf': 0.7989
        },
        {
            'name': 'Illumina INDEL',
            'loader_func': loader.load_illumina_indel,
            'features': loader.ILLUMINA_FEATURES,
            'paper_rf': 0.9284
        },
        {
            'name': 'Ion SNV',
            'loader_func': loader.load_ion_snp,
            'features': loader.ION_FEATURES,
            'paper_rf': 0.9749
        },
        {
            'name': 'Ion INDEL',
            'loader_func': loader.load_ion_indel,
            'features': loader.ION_INDEL_FEATURES,
            'paper_rf': 0.9383
        }
    ]
    
    for config in models_config:
        print(f"\nProcessing: {config['name']}")
        print("-" * 30)
        
        # 1. Load and Preprocess Data
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
        
        # 3. Train Random Forest (Exact Parameters)
        print("  > Training H2O Random Forest...")
        rf_model = H2ORandomForestEstimator(
            ntrees=500,
            max_depth=50,
            stopping_rounds=5,
            stopping_tolerance=1e-2,  # Set to 0.001 based on interpretation
            stopping_metric='logloss',
            seed=42,
            score_each_iteration=False,
            score_tree_interval = 0
        )
        
        rf_model.train(
            x=config['features'], 
            y='target',
            training_frame=train_h2o,
            validation_frame=val_h2o
        )
        
        # 4. Evaluate
        perf = rf_model.model_performance(test_h2o)
        rf_auc = perf.auc()
        
        print(f"  > Done.")
        print(f"  > Our RF AUC:   {rf_auc:.4f}")
        print(f"  > Paper RF AUC: {config['paper_rf']:.4f}")
        
        print(f"  > Done.")
        print(f"  > Our RF AUC:   {rf_auc:.4f}")
        print(f"  > Paper RF AUC: {config['paper_rf']:.4f}")
        
        # SAVE MODEL (ADD THIS)
        model_dir = Path("../models")
        model_dir.mkdir(exist_ok=True)
        rf_path = h2o.save_model(rf_model, path=str(model_dir), 
                                 filename=f"rf_{safe_name.lower()}", force=True)
        print(f"  > Saved model: rf_{safe_name.lower()}")
        
        results.append({
            'Dataset': config['name'],
            'Our_RF_AUC': rf_auc,
            'Paper_RF_AUC': config['paper_rf'],
            'Difference': rf_auc - config['paper_rf']
        })

    # Summary Output
    print("\n" + "="*70)
    print("FINAL SUMMARY: RANDOM FOREST REPLICATION")
    print("="*70)
    
    df = pd.DataFrame(results)
    print(df.to_string(index=False, float_format=lambda x: "{:.4f}".format(x)))
    
    # Optional: Save just these results
    df.to_csv('results/random_forest_replication.csv', index=False)
    print(f"\nSaved results to /results/random_forest_replication.csv")

    h2o.cluster().shutdown()

if __name__ == "__main__":
    train_random_forests()