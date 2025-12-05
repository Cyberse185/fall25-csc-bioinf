import pandas as pd

print("="*70)
print("PREPARING SNV AND INDEL DATASETS")
print("="*70)

# Load GIAB truth
print("\nLoading GIAB truth...")
with open('GIAB_positions.txt', 'r') as f:
    giab_positions = set(line.strip() for line in f)
print(f"Loaded {len(giab_positions):,} true variants")

def split_snv_indel(feature_file, output_prefix):
    """Split feature file into SNV and INDEL"""
    print(f"\nProcessing {feature_file}...")
    df = pd.read_csv(feature_file, sep='\t')
    
    # Create position for GIAB labeling
    df['position'] = df['CHROM'] + ':' + df['POS'].astype(str) + ':' + df['REF'] + ':' + df['ALT']
    df['target'] = df['position'].isin(giab_positions).astype(int)
    
    # Determine if SNV or INDEL
    # SNV: REF and ALT are both single nucleotides
    df['is_snv'] = (df['REF'].str.len() == 1) & (df['ALT'].str.len() == 1)
    
    # Split
    df_snv = df[df['is_snv']].copy()
    df_indel = df[~df['is_snv']].copy()
    
    # Drop helper columns
    df_snv = df_snv.drop(['position', 'is_snv'], axis=1)
    df_indel = df_indel.drop(['position', 'is_snv'], axis=1)
    
    # Save
    snv_file = f"{output_prefix}_SNV.tsv"
    indel_file = f"{output_prefix}_INDEL.tsv"
    
    df_snv.to_csv(snv_file, sep='\t', index=False)
    df_indel.to_csv(indel_file, sep='\t', index=False)
    
    print(f"  SNVs:   {len(df_snv):,} (TRUE: {df_snv['target'].sum():,}) → {snv_file}")
    print(f"  INDELs: {len(df_indel):,} (TRUE: {df_indel['target'].sum():,}) → {indel_file}")
    
    return len(df_snv), len(df_indel)

print("\n" + "="*70)
print("FULL COVERAGE SAMPLES")
print("="*70)

split_snv_indel('SRR098401_full_features.table', 'SRR098401_full')
split_snv_indel('SRR1611180_full_features.table', 'SRR1611180_full')
split_snv_indel('SRR3197786_full_features.table', 'SRR3197786_full')

print("\n" + "="*70)
print("30X COVERAGE SAMPLES")
print("="*70)

split_snv_indel('SRR098401_30x_features.table', 'SRR098401_30x')
split_snv_indel('SRR1611180_30x_features.table', 'SRR1611180_30x')
split_snv_indel('SRR3197786_30x_features.table', 'SRR3197786_30x')

print("\n✓ All datasets split into SNV and INDEL!")
print("\nNext: Run train_h2o_snv_indel.py to train separate models")