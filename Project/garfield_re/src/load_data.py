"""Load and preprocess GARFIELD-NGS datasets"""

import pandas as pd
import numpy as np
from pathlib import Path

class GarfieldDataLoader:
    """Load and preprocess GARFIELD-NGS training data"""
    
    ILLUMINA_FEATURES = [
        "BaseQRankSum", "ReadPosRankSum", "DP", "FS", "MQ", 
        "MQRankSum", "QD", "SOR", "QUAL", "GQ"
    ]
    
    ION_FEATURES = [
        "FDP", "FAO", "QD", "FSAF", "FSAR", "FXX", "LEN", "HRUN",
        "RBI", "VARB", "STB", "STBP", "PB", "PBP", "MLLD", "SSSB", "QUAL", "GQ"
    ]
    
    # Ion INDEL is missing PB and PBP, but has SSEN and SSEP
    ION_INDEL_FEATURES = [
        "FDP", "FAO", "QD", "FSAF", "FSAR", "FXX", "LEN", "HRUN",
        "RBI", "VARB", "STB", "STBP", "MLLD", "SSEN", "SSEP", "SSSB", "QUAL", "GQ"
    ]
    
    def __init__(self, data_dir="GARFIELD-NGS/datasets"):
        self.data_dir = Path(data_dir)
        if not self.data_dir.exists():
            raise FileNotFoundError(f"Data directory not found: {data_dir}")
    
    def load_dataset(self, filename):
        filepath = self.data_dir / filename
        print(f"Loading {filename}...")
        df = pd.read_csv(filepath, sep='\t')
        print(f"  Shape: {df.shape}")
        if 'Class' in df.columns:
            print(f"  Classes: {df['Class'].value_counts().to_dict()}")
        return df
    
    def load_illumina_snp(self):
        print("\n" + "="*60)
        print("Loading Illumina SNP Datasets")
        print("="*60)
        return {
            'train': self.load_dataset('ILM_SNP_Training.balance.txt'),
            'validation': self.load_dataset('ILM_SNP_Validation1.txt'),
            'test': self.load_dataset('ILM_SNP_Test.txt'),
            'pretrain': self.load_dataset('ILM_SNP_Pre_training.balance.txt')
        }
    
    def load_illumina_indel(self):
        print("\n" + "="*60)
        print("Loading Illumina INDEL Datasets")
        print("="*60)
        return {
            'train': self.load_dataset('ILM_INDEL_Training.txt'),
            'validation': self.load_dataset('ILM_INDEL_Validation1.txt'),
            'test': self.load_dataset('ILM_INDEL_Test.txt'),
            'pretrain': self.load_dataset('ILM_INDEL_Pre_training.txt')
        }
    
    def load_ion_snp(self):
        print("\n" + "="*60)
        print("Loading ION SNP Datasets")
        print("="*60)
        return {
            'train': self.load_dataset('ION_SNP_Training.balance.txt'),
            'validation': self.load_dataset('ION_SNP_Validation1.txt'),
            'test': self.load_dataset('ION_SNP_Test.txt'),
            'pretrain': self.load_dataset('ION_SNP_Pre_training.balance.txt')
        }
    
    def load_ion_indel(self):
        print("\n" + "="*60)
        print("Loading ION INDEL Datasets")
        print("="*60)
        return {
            'train': self.load_dataset('ION_INDEL_Training.txt'),
            'validation': self.load_dataset('ION_INDEL_Validation1.txt'),
            'test': self.load_dataset('ION_INDEL_Test.txt'),
            'pretrain': self.load_dataset('ION_INDEL_Pre_training.txt')
        }
    
    def preprocess_data(self, df, features):
        X = df[features].copy()
        y = (df['Class'] == 'T').astype(int)
        
        # Handle missing values
        missing = X.isnull().sum()
        if missing.sum() > 0:
            print(f"\nFilling {missing.sum()} missing values with median")
            X = X.fillna(X.median())
        
        # Handle infinite values
        X = X.replace([np.inf, -np.inf], np.nan)
        X = X.fillna(X.median())
        
        print(f"\nFinal: X={X.shape}, y={y.shape}")
        print(f"Class balance: {y.value_counts().to_dict()}")
        return X, y

def main():
    print("GARFIELD-NGS Data Loader")
    print("="*60)
    
    loader = GarfieldDataLoader()
    ilm_snp = loader.load_illumina_snp()
    
    X_train, y_train = loader.preprocess_data(
        ilm_snp['train'], 
        loader.ILLUMINA_FEATURES
    )
    
    print("\n✓ Data loaded successfully!")
    print(f"Training samples: {len(X_train)}")
    print(f"Features: {list(X_train.columns)}")

if __name__ == "__main__":
    main()