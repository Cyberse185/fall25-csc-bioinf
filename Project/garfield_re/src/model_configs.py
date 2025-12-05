"""
Model configurations from GARFIELD-NGS paper (Supplementary Table S2)
"""

MODEL_CONFIGS = {
    'illumina_snv': {
        'activation': 'Tanh',
        'hidden': [20, 60, 60, 50, 70],
        'rho': 0.956,
        'epsilon': 1.00e-10,
        'l1': 0.0022,
        'l2': 0.00293,
        'mini_batch_size': 1
    },
    'illumina_indel': {
        'activation': 'Rectifier',
        'hidden': [10, 70, 60, 10, 10],
        'rho': 0.961,
        'epsilon': 1.00e-9,
        'l1': 0.00015,
        'l2': 0.01311,
        'mini_batch_size': 1
    },
    'ion_snv': {
        'activation': 'Tanh',
        'hidden': [40, 90, 90, 30, 100],
        'rho': 0.978,
        'epsilon': 1.00e-10,
        'l1': 0.0128,
        'l2': 0.00224,
        'mini_batch_size': 1
    },
    'ion_indel': {
        'activation': 'Rectifier',
        'hidden': [100, 40, 50, 90, 90],
        'rho': 0.985,
        'epsilon': 1.00e-9,
        'l1': 0.00026,
        'l2': 0.05488,
        'mini_batch_size': 1
    }
}

FILTERING_THRESHOLDS = {
    'illumina_snv': 0.025,
    'illumina_indel': 0.630,
    'ion_snv': 0.139,
    'ion_indel': 0.320
}

TRAINING_PARAMS = {
    'epochs': 1000,
    'stopping_rounds': 5,
    'stopping_tolerance': 1e-3,
    'stopping_metric': 'logloss'
}

def get_model_config(model_type):
    if model_type not in MODEL_CONFIGS:
        raise ValueError(f"Unknown model type: {model_type}")
    config = MODEL_CONFIGS[model_type].copy()
    config.update(TRAINING_PARAMS)
    config['threshold'] = FILTERING_THRESHOLDS[model_type]
    return config