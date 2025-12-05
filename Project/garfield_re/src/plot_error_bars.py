import numpy as np
import matplotlib.pyplot as plt
import math

# ============================================================================
# 1. DATA INPUT (From your logs)
# ============================================================================

# Sample sizes from your log output
# Full: "TRUE: 762,734, FALSE: 188,169"
n_full_pos = 762734
n_full_neg = 188169

# 30X: "TRUE: 226,200, FALSE: 80,386"
n_30x_pos = 226200
n_30x_neg = 80386

# Your Results
models = ['LR', 'RF', 'GARFIELD']
auc_full = [0.7499, 0.7474, 0.7790]
auc_30x  = [0.7344, 0.7598, 0.7754]

# ============================================================================
# 2. CALCULATE ERROR BARS (Hanley & McNeil Method)
# ============================================================================
def calculate_se(auc, n_pos, n_neg):
    """
    Calculate Standard Error of AUC using Hanley & McNeil (1982)
    This accounts for the correlation between positive and negative cases.
    """
    q1 = auc / (2 - auc)
    q2 = 2 * auc**2 / (1 + auc)
    
    numerator = (auc * (1 - auc)) + \
                ((n_pos - 1) * (q1 - auc**2)) + \
                ((n_neg - 1) * (q2 - auc**2))
    
    denominator = n_pos * n_neg
    se = math.sqrt(numerator / denominator)
    return se

# Calculate error (95% CI = 1.96 * SE)
err_full = []
err_30x = []

print(f"{'Model':<10} {'Set':<5} {'AUC':<8} {'Error (95% CI)'}")
print("-" * 45)

for i, m in enumerate(models):
    # Full Error
    se_f = calculate_se(auc_full[i], n_full_pos, n_full_neg)
    ci_f = 1.96 * se_f
    err_full.append(ci_f)
    print(f"{m:<10} Full  {auc_full[i]:<8.4f} ±{ci_f:.5f}")
    
    # 30X Error
    se_3 = calculate_se(auc_30x[i], n_30x_pos, n_30x_neg)
    ci_3 = 1.96 * se_3
    err_30x.append(ci_3)
    print(f"{m:<10} 30X   {auc_30x[i]:<8.4f} ±{ci_3:.5f}")

# ============================================================================
# 3. PLOT WITH ERROR BARS
# ============================================================================
fig, ax = plt.subplots(figsize=(10, 6))

x = np.arange(len(models))
width = 0.35

# Plot Full
rects1 = ax.bar(x - width/2, auc_full, width, label='Full Coverage', 
                color='#2ca02c', alpha=0.9,
                yerr=err_full, capsize=5, error_kw={'ecolor': 'black', 'elinewidth': 1.5})

# Plot 30X
rects2 = ax.bar(x + width/2, auc_30x, width, label='30X Coverage', 
                color='#ff7f0e', alpha=0.9,
                yerr=err_30x, capsize=5, error_kw={'ecolor': 'black', 'elinewidth': 1.5})

ax.set_ylabel('AUROC', fontsize=12)
ax.set_title('SNV Model Robustness with 95% Confidence Intervals', fontsize=14, fontweight='bold')
ax.set_xticks(x)
ax.set_xticklabels(models, fontsize=11)
ax.legend(loc='lower right')
ax.grid(axis='y', alpha=0.3)

# Zoom in to make bars visible (Standard Error is tiny!)
ax.set_ylim([0.70, 0.80])

# Add text labels on top
def autolabel(rects):
    for rect in rects:
        height = rect.get_height()
        ax.annotate(f'{height:.3f}',
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 10),  # 10 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom', fontsize=9, fontweight='bold')

autolabel(rects1)
autolabel(rects2)

plt.tight_layout()
plt.savefig('snv_robustness_with_error_bars.png', dpi=300)
print("\n✓ Plot saved to snv_robustness_with_error_bars.png")