import pandas as pd
import numpy as np
from pathlib import Path
import os

class GarfieldDataLoader:
    def __init__(self):
        # ---------------------------------------------------------
        # ROBUST PATH FINDING
        # ---------------------------------------------------------
        # 1. Get the folder where THIS script (load_data.py) lives
        script_dir = Path(__file__).resolve().parent
        
        # 2. Look for 'data' folder relative to the script
        # If script is in 'src/', data should be in '../data'
        project_root = script_dir.parent
        data_path_1 = project_root / "data"
        
        # 3. Fallback: maybe we are running from root?
        data_path_2 = Path("data")
        
        if data_path_1.exists() and data_path_1.is_dir():
            self.data_dir = data_path_1
            print(f"✓ Found data at: {self.data_dir}")
        elif data_path_2.exists() and data_path_2.is_dir():
            self.data_dir = data_path_2
            print(f"✓ Found data at: {self.data_dir}")
        else:
            # Debug info to help you fix it if it fails
            raise FileNotFoundError(
                f"CRITICAL ERROR: Could not find 'data' folder.\n"
                f"  Checked: {data_path_1}\n"
                f"  Checked: {data_path_2}\n"
                f"  Current Working Dir: {os.getcwd()}"
            )

        # ---------------------------------------------------------
        # FEATURES
        # ---------------------------------------------------------
        self.ILLUMINA_FEATURES = [
            "BaseQRankSum", "ReadPosRankSum", "DP", "FS", "MQ",
            "MQRankSum", "QD", "SOR", "QUAL", "GQ"
        ]
        
        self.ION_FEATURES = [
            "FDP", "FAO", "QD", "FSAF", "FSAR", "FXX", "LEN", "HRUN",
            "RBI", "VARB", "STB", "STBP", "PB", "PBP", "MLLD", "SSSB", "QUAL", "GQ"
        ]
        
        self.ION_INDEL_FEATURES = [
            "FDP", "FAO", "QD", "FSAF", "FSAR", "FXX", "LEN", "HRUN",
            "RBI", "VARB", "STB", "STBP", "MLLD", "SSEN", "SSEP", "SSSB", "QUAL", "GQ"
        ]

    def load_dataset(self, filename):
        """Load a dataset from the data directory"""
        file_path = self.data_dir / filename
        
        if not file_path.exists():
            raise FileNotFoundError(f"File not found: {file_path}")
            
        print(f"Loading {filename}...")
        
        try:
            df = pd.read_csv(file_path, sep='\t')
        except:
            df = pd.read_csv(file_path, sep=r'\s+')
            
        df = df.replace('.', np.nan)
        
        all_features = set(self.ILLUMINA_FEATURES + self.ION_FEATURES + self.ION_INDEL_FEATURES)
        for col in df.columns:
            if col in all_features:
                df[col] = pd.to_numeric(df[col], errors='coerce')
        
        return df