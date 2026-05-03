# ============================================================
#  HBV Multi-omics — Machine Learning Pipeline
#  Google Colab Script
#  Target: Q1 Journal Publication
# ============================================================
#
#  INPUT FILES (upload to Colab):
#    - cohort_D_ELISA_miRNA_ml_ready.csv   (n=77,  17 features)
#    - cohort_E_ALL_ml_ready.csv            (n=73,  28 features)
#
#  MODELS: Logistic Regression, Random Forest,
#          XGBoost, LightGBM
#
#  METHOD: Stratified 5-fold CV inside sklearn Pipeline
#          (scaler fitted on train fold only → no leakage)
#          Bootstrap 95% CI on AUC
#          SHAP values for best model
# ============================================================


# ============================================================
# CELL 1 — Install & Import
# ============================================================

# !pip install xgboost lightgbm shap -q

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
import shap
import warnings
warnings.filterwarnings('ignore')

from sklearn.pipeline          import Pipeline
from sklearn.preprocessing     import StandardScaler, label_binarize
from sklearn.linear_model      import LogisticRegression
from sklearn.ensemble          import RandomForestClassifier
from sklearn.model_selection   import (StratifiedKFold, cross_validate,
                                        cross_val_predict, GridSearchCV)
from sklearn.metrics           import (roc_auc_score, f1_score,
                                        accuracy_score, confusion_matrix,
                                        roc_curve, classification_report)
from sklearn.utils             import resample
from xgboost                   import XGBClassifier
from lightgbm                  import LGBMClassifier

np.random.seed(42)

# Publication plot style
plt.rcParams.update({
    'font.family':     'DejaVu Sans',
    'font.size':       11,
    'axes.titlesize':  12,
    'axes.labelsize':  11,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'figure.dpi':      150,
    'savefig.dpi':     300,
    'savefig.bbox':    'tight',
})

COLORS = {
    'LR':    '#4393C3',
    'RF':    '#2CA02C',
    'XGB':   '#D62728',
    'LGBM':  '#FF7F0E',
}

print("All libraries loaded successfully.")


# ============================================================
# CELL 2 — Load Data & Define Feature Sets
# ============================================================

# Load cohorts
cD = pd.read_csv('cohort_D_ELISA_miRNA_ml_ready.csv')
cE = pd.read_csv('cohort_E_ALL_ml_ready.csv')

print(f"Cohort D (ELISA+miRNA): {cD.shape}")
print(f"  Group 0 (Control): {(cD['Group']==0).sum()}")
print(f"  Group 1 (HBV):     {(cD['Group']==1).sum()}")
print(f"\nCohort E (ALL):         {cE.shape}")
print(f"  Group 0 (Control): {(cE['Group']==0).sum()}")
print(f"  Group 1 (HBV):     {(cE['Group']==1).sum()}")

# ── Feature sets ─────────────────────────────────────────────
# From R statistical analysis:
#   ALL 5 ELISA cytokines → significant
#   9/12 miRNAs → significant (drop miR-21, miR-125, miR-155)
#   2/11 SNPs   → significant (TGFb-509, TNFa-863)

ELISA_COLS = ['IL10_ELISA', 'TGFb_ELISA', 'IL6_level',
              'TNFA_ELISA', 'IFN_ELISA']

MIRNA_SIG  = ['mir_10', 'mir_17', 'mir_24', 'mir_26',
              'mir_122', 'mir_145', 'mir_146', 'mir_148', 'mir_221']

MIRNA_ALL  = ['mir_10', 'mir_17', 'mir_21', 'mir_24', 'mir_26',
              'mir_122', 'mir_125', 'mir_145', 'mir_146',
              'mir_148', 'mir_155', 'mir_221']

SNP_SIG    = ['TGFb_509_Geno', 'TNFA_863_Geno']

SNP_ALL    = ['IL10_1082_Geno', 'IL10_819_Geno', 'TGFb_800_Geno',
              'TGFb_509_Geno', 'TGFb_codon10_Geno', 'TGFb_codon25_Geno',
              'TNFA_863_Geno', 'TNFA_376_Geno', 'TNFA_308_Geno',
              'TNFA_857_Geno', 'TNFA_489_Geno']

# ── Define experiment sets ────────────────────────────────────
EXPERIMENTS = {
    'ELISA only':          (cD, ELISA_COLS),
    'miRNA only':          (cD, MIRNA_SIG),
    'ELISA + miRNA':       (cD, ELISA_COLS + MIRNA_SIG),
    'ELISA + miRNA (ALL)': (cD, ELISA_COLS + MIRNA_ALL),
    'ALL features':        (cE, SNP_SIG + ELISA_COLS + MIRNA_SIG),
}

# Main experiment for detailed analysis
X_main = cD[ELISA_COLS + MIRNA_SIG].values
y_main = cD['Group'].values

feature_names_main = (
    ['IL-10', 'TGF-β', 'IL-6', 'TNF-α', 'IFN-γ'] +
    ['miR-10', 'miR-17', 'miR-24', 'miR-26',
     'miR-122', 'miR-145', 'miR-146', 'miR-148', 'miR-221']
)

print(f"\nMain experiment features: {len(feature_names_main)}")
print(f"Feature names: {feature_names_main}")


# ============================================================
# CELL 3 — Define Models & CV
# ============================================================

CV = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

def make_pipeline(model):
    """Wrap model in Pipeline with StandardScaler.
    Scaler is fitted only on train fold → no data leakage."""
    return Pipeline([
        ('scaler', StandardScaler()),
        ('clf',    model)
    ])

MODELS = {
    'LR': make_pipeline(
        LogisticRegression(
            C=1.0, max_iter=1000,
            class_weight='balanced',
            random_state=42
        )
    ),
    'RF': make_pipeline(
        RandomForestClassifier(
            n_estimators=500,
            max_depth=None,
            min_samples_leaf=2,
            class_weight='balanced',
            random_state=42,
            n_jobs=-1
        )
    ),
    'XGB': make_pipeline(
        XGBClassifier(
            n_estimators=300,
            max_depth=4,
            learning_rate=0.05,
            subsample=0.8,
            colsample_bytree=0.8,
            eval_metric='logloss',
            use_label_encoder=False,
            random_state=42,
            verbosity=0
        )
    ),
    'LGBM': make_pipeline(
        LGBMClassifier(
            n_estimators=300,
            max_depth=4,
            learning_rate=0.05,
            subsample=0.8,
            colsample_bytree=0.8,
            class_weight='balanced',
            random_state=42,
            verbose=-1,
            n_jobs=-1
        )
    ),
}

print("Models defined:")
for name in MODELS:
    print(f"  {name}: {MODELS[name].named_steps['clf'].__class__.__name__}")


# ============================================================
# CELL 4 — Bootstrap AUC Function
# ============================================================

def bootstrap_auc(y_true, y_prob, n_boot=1000, ci=0.95, seed=42):
    """Compute bootstrap confidence interval for AUC."""
    rng   = np.random.RandomState(seed)
    aucs  = []
    n     = len(y_true)
    for _ in range(n_boot):
        idx = rng.choice(n, n, replace=True)
        if len(np.unique(y_true[idx])) < 2:
            continue
        aucs.append(roc_auc_score(y_true[idx], y_prob[idx]))
    aucs = np.array(aucs)
    alpha = (1 - ci) / 2
    return (np.percentile(aucs, alpha * 100),
            np.percentile(aucs, (1 - alpha) * 100))


# ============================================================
# CELL 5 — Main Cross-Validation (ELISA + miRNA significant)
# ============================================================

print("=" * 60)
print("MAIN EXPERIMENT: ELISA + significant miRNA (n=77)")
print("Stratified 5-fold CV | No data leakage")
print("=" * 60)

cv_results    = {}
oof_probs     = {}
oof_preds     = {}

for name, pipe in MODELS.items():
    # Cross-validated predicted probabilities (out-of-fold)
    oof_prob = cross_val_predict(
        pipe, X_main, y_main,
        cv=CV, method='predict_proba', n_jobs=-1
    )[:, 1]

    oof_pred = cross_val_predict(
        pipe, X_main, y_main,
        cv=CV, method='predict', n_jobs=-1
    )

    auc  = roc_auc_score(y_main, oof_prob)
    f1   = f1_score(y_main, oof_pred, average='macro')
    acc  = accuracy_score(y_main, oof_pred)
    ci_l, ci_u = bootstrap_auc(y_main, oof_prob)

    oof_probs[name] = oof_prob
    oof_preds[name] = oof_pred

    cv_results[name] = {
        'AUC':    round(auc,  4),
        'F1':     round(f1,   4),
        'Acc':    round(acc,  4),
        'CI_low': round(ci_l, 4),
        'CI_up':  round(ci_u, 4),
    }

    print(f"\n{name}:")
    print(f"  AUC  = {auc:.4f}  95% CI [{ci_l:.4f}, {ci_u:.4f}]")
    print(f"  F1   = {f1:.4f}")
    print(f"  Acc  = {acc:.4f}")

results_df = pd.DataFrame(cv_results).T
results_df.index.name = 'Model'
results_df.reset_index(inplace=True)
print("\n── Summary Table ──")
print(results_df.to_string(index=False))
results_df.to_csv('ML_results_main.csv', index=False)


# ============================================================
# CELL 6 — Feature Set Comparison
# ============================================================

print("\n" + "=" * 60)
print("FEATURE SET COMPARISON")
print("=" * 60)

# Use RF as reference model (good balance of performance/stability)
ref_model_name = 'RF'
ref_clf = RandomForestClassifier(
    n_estimators=500, min_samples_leaf=2,
    class_weight='balanced', random_state=42, n_jobs=-1
)

feat_comparison = {}

for exp_name, (cohort, feat_cols) in EXPERIMENTS.items():
    X = cohort[feat_cols].values
    y = cohort['Group'].values

    pipe = make_pipeline(ref_clf.__class__(
        n_estimators=500, min_samples_leaf=2,
        class_weight='balanced', random_state=42, n_jobs=-1
    ))

    oof_prob = cross_val_predict(
        pipe, X, y, cv=StratifiedKFold(5, shuffle=True, random_state=42),
        method='predict_proba', n_jobs=-1
    )[:, 1]

    auc = roc_auc_score(y, oof_prob)
    ci_l, ci_u = bootstrap_auc(y, oof_prob)

    feat_comparison[exp_name] = {
        'n': len(y),
        'n_features': len(feat_cols),
        'AUC': round(auc, 4),
        'CI': f"[{ci_l:.3f}, {ci_u:.3f}]"
    }
    print(f"  {exp_name:25s} n={len(y):3d}  "
          f"features={len(feat_cols):2d}  "
          f"AUC={auc:.4f}  CI {ci_l:.3f}–{ci_u:.3f}")

feat_df = pd.DataFrame(feat_comparison).T
feat_df.to_csv('ML_feature_comparison.csv')


# ============================================================
# CELL 7 — FIGURE A: ROC Curves (4 models)
# ============================================================

fig, ax = plt.subplots(figsize=(6, 5.5))

for name in MODELS:
    fpr, tpr, _ = roc_curve(y_main, oof_probs[name])
    auc  = cv_results[name]['AUC']
    ci_l = cv_results[name]['CI_low']
    ci_u = cv_results[name]['CI_up']
    ax.plot(fpr, tpr, color=COLORS[name], linewidth=2,
            label=f"{name}  AUC={auc:.3f} [{ci_l:.3f}–{ci_u:.3f}]")

ax.plot([0,1],[0,1], 'k--', linewidth=0.8, alpha=0.5)
ax.fill_between([0,1],[0,1],[0,1], alpha=0.04, color='gray')
ax.set_xlabel('False Positive Rate (1 – Specificity)')
ax.set_ylabel('True Positive Rate (Sensitivity)')
ax.set_title('ROC Curves — 5-fold CV\n(ELISA + miRNA panel, n=77)',
             fontsize=12, fontweight='bold')
ax.legend(loc='lower right', fontsize=9.5, framealpha=0.9)
ax.set_xlim([-0.02, 1.02])
ax.set_ylim([-0.02, 1.05])
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('FigureML_A_ROC_curves.png', dpi=300)
plt.savefig('FigureML_A_ROC_curves.pdf', dpi=300)
plt.show()
print("Figure A (ROC curves) saved.")


# ============================================================
# CELL 8 — FIGURE B: Confusion Matrices (2x2 grid)
# ============================================================

fig, axes = plt.subplots(2, 2, figsize=(8, 7))
axes = axes.flatten()

labels_cm = ['Control', 'HBV']

for idx, (name, ax) in enumerate(zip(MODELS, axes)):
    cm = confusion_matrix(y_main, oof_preds[name])
    cm_pct = cm.astype(float) / cm.sum(axis=1, keepdims=True) * 100

    sns.heatmap(cm, annot=False, fmt='d', ax=ax,
                cmap='Blues', linewidths=0.5,
                xticklabels=labels_cm,
                yticklabels=labels_cm,
                cbar=False)

    # Annotate with count + %
    for i in range(2):
        for j in range(2):
            ax.text(j + 0.5, i + 0.4,
                    f"{cm[i,j]}",
                    ha='center', va='center',
                    fontsize=14, fontweight='bold',
                    color='white' if cm_pct[i,j] > 50 else 'black')
            ax.text(j + 0.5, i + 0.65,
                    f"({cm_pct[i,j]:.1f}%)",
                    ha='center', va='center',
                    fontsize=9,
                    color='white' if cm_pct[i,j] > 50 else 'black')

    tn, fp, fn, tp = cm.ravel()
    sens = tp / (tp + fn) * 100
    spec = tn / (tn + fp) * 100

    ax.set_title(
        f"{name}  |  AUC={cv_results[name]['AUC']:.3f}\n"
        f"Sens={sens:.1f}%  Spec={spec:.1f}%",
        fontsize=10, fontweight='bold'
    )
    ax.set_xlabel('Predicted', fontsize=10)
    ax.set_ylabel('Actual',    fontsize=10)

plt.suptitle('Confusion Matrices — 5-fold CV (n=77)',
             fontsize=12, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig('FigureML_B_Confusion_matrices.png', dpi=300,
            bbox_inches='tight')
plt.savefig('FigureML_B_Confusion_matrices.pdf', dpi=300,
            bbox_inches='tight')
plt.show()
print("Figure B (Confusion matrices) saved.")


# ============================================================
# CELL 9 — FIGURE C: Model Comparison Bar Chart
# ============================================================

metrics = ['AUC', 'F1', 'Acc']
metric_labels = ['AUC', 'Macro F1', 'Accuracy']
x = np.arange(len(metrics))
width = 0.18

fig, ax = plt.subplots(figsize=(8, 5))

for i, (name, res) in enumerate(cv_results.items()):
    vals = [res['AUC'], res['F1'], res['Acc']]
    bars = ax.bar(x + i * width, vals, width,
                  label=name, color=COLORS[name],
                  alpha=0.85, edgecolor='white', linewidth=0.5)

    # Error bars for AUC only
    ax.errorbar(
        x[0] + i * width, res['AUC'],
        yerr=[[res['AUC'] - res['CI_low']],
              [res['CI_up']  - res['AUC']]],
        fmt='none', color='black', capsize=4, linewidth=1.2
    )

ax.set_xticks(x + width * 1.5)
ax.set_xticklabels(metric_labels, fontsize=11)
ax.set_ylabel('Score', fontsize=11)
ax.set_title('Model Performance Comparison — 5-fold CV\n'
             '(Error bars = 95% bootstrap CI on AUC)',
             fontsize=12, fontweight='bold')
ax.legend(title='Model', fontsize=10, title_fontsize=10,
          framealpha=0.9)
ax.set_ylim([0.5, 1.05])
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.grid(axis='y', alpha=0.3)
ax.axhline(y=1.0, color='gray', linestyle='--',
           linewidth=0.6, alpha=0.5)

plt.tight_layout()
plt.savefig('FigureML_C_Model_comparison.png', dpi=300)
plt.savefig('FigureML_C_Model_comparison.pdf', dpi=300)
plt.show()
print("Figure C (Model comparison) saved.")


# ============================================================
# CELL 10 — FIGURE D: Feature Set Comparison
# ============================================================

feat_names_plot = list(feat_comparison.keys())
feat_aucs = [feat_comparison[k]['AUC'] for k in feat_names_plot]
feat_ns   = [feat_comparison[k]['n']   for k in feat_names_plot]

# Parse CI
feat_ci_l, feat_ci_u = [], []
for k in feat_names_plot:
    ci_str = feat_comparison[k]['CI']
    lo, hi = ci_str.strip('[]').split(', ')
    feat_ci_l.append(float(lo))
    feat_ci_u.append(float(hi))

feat_ci_l = np.array(feat_ci_l)
feat_ci_u = np.array(feat_ci_u)
feat_aucs_arr = np.array(feat_aucs)

colors_feat = ['#4393C3', '#2CA02C', '#D62728', '#FF7F0E', '#9467BD']

fig, ax = plt.subplots(figsize=(8, 4.5))
bars = ax.barh(feat_names_plot, feat_aucs, color=colors_feat,
               alpha=0.85, edgecolor='white', linewidth=0.5)

ax.errorbar(feat_aucs_arr, range(len(feat_names_plot)),
            xerr=[feat_aucs_arr - feat_ci_l, feat_ci_u - feat_aucs_arr],
            fmt='none', color='black', capsize=4, linewidth=1.2)

for i, (auc, n) in enumerate(zip(feat_aucs, feat_ns)):
    ax.text(auc + 0.005, i, f"{auc:.3f}  (n={n})",
            va='center', fontsize=9.5)

ax.set_xlabel('AUC (5-fold CV)', fontsize=11)
ax.set_title('Feature Set Comparison — Random Forest\n'
             '(Error bars = 95% bootstrap CI)',
             fontsize=12, fontweight='bold')
ax.set_xlim([0.5, 1.12])
ax.axvline(x=0.7, color='gray', linestyle='--',
           linewidth=0.8, alpha=0.5, label='AUC=0.7 reference')
ax.axvline(x=0.9, color='green', linestyle='--',
           linewidth=0.8, alpha=0.5, label='AUC=0.9 reference')
ax.legend(fontsize=9, loc='lower right')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.grid(axis='x', alpha=0.3)

plt.tight_layout()
plt.savefig('FigureML_D_Feature_comparison.png', dpi=300)
plt.savefig('FigureML_D_Feature_comparison.pdf', dpi=300)
plt.show()
print("Figure D (Feature set comparison) saved.")


# ============================================================
# CELL 11 — SHAP Analysis (Best model = RF or XGB)
# ============================================================

print("\nRunning SHAP analysis on Random Forest...")

# Fit RF on full data (for SHAP — we use full data here
# because SHAP is for interpretation, not evaluation)
rf_full = Pipeline([
    ('scaler', StandardScaler()),
    ('clf', RandomForestClassifier(
        n_estimators=500, min_samples_leaf=2,
        class_weight='balanced', random_state=42, n_jobs=-1
    ))
])
rf_full.fit(X_main, y_main)

# Get scaled features (SHAP needs post-scaler data)
X_scaled = rf_full.named_steps['scaler'].transform(X_main)

# SHAP TreeExplainer
explainer  = shap.TreeExplainer(rf_full.named_steps['clf'])
shap_vals  = explainer.shap_values(X_scaled)

# For binary: shap_values is list[2]; take class 1 (HBV)
if isinstance(shap_vals, list):
    sv = shap_vals[1]
else:
    sv = shap_vals

shap_df = pd.DataFrame(sv, columns=feature_names_main)
mean_abs_shap = shap_df.abs().mean().sort_values(ascending=False)

print("\nTop features by mean |SHAP|:")
print(mean_abs_shap.to_string())

mean_abs_shap.to_csv('SHAP_feature_importance.csv', header=['mean_abs_SHAP'])


# ============================================================
# CELL 12 — FIGURE E: SHAP Summary Plot
# ============================================================

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Left: Bar plot (mean |SHAP|)
ax = axes[0]
top_features = mean_abs_shap.head(14)
colors_shap  = ['#D62728' if f in ['IL-10','IFN-γ','TNF-α','IL-6','TGF-β']
                else '#4393C3' for f in top_features.index]

bars = ax.barh(top_features.index[::-1],
               top_features.values[::-1],
               color=colors_shap[::-1],
               alpha=0.85, edgecolor='white')

ax.set_xlabel('Mean |SHAP value|', fontsize=11)
ax.set_title('Feature Importance\n(Mean |SHAP| — Random Forest)',
             fontsize=11, fontweight='bold')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.grid(axis='x', alpha=0.3)

# Legend
from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor='#D62728', alpha=0.85, label='Cytokine'),
    Patch(facecolor='#4393C3', alpha=0.85, label='miRNA'),
]
ax.legend(handles=legend_elements, fontsize=9,
          loc='lower right')

# Right: SHAP dot plot (beeswarm style — manual)
ax2 = axes[1]
top_n   = 12
top_idx = [list(feature_names_main).index(f)
           for f in mean_abs_shap.head(top_n).index
           if f in feature_names_main]
top_names = [feature_names_main[i] for i in top_idx]

for rank, (fname, fidx) in enumerate(zip(top_names, top_idx)):
    vals  = sv[:, fidx]
    feats = X_scaled[:, fidx]

    # Color by feature value (low=blue, high=red)
    norm  = plt.Normalize(feats.min(), feats.max())
    cmap  = plt.cm.RdBu_r
    colors_dot = cmap(norm(feats))

    # Jitter y slightly
    jitter = np.random.uniform(-0.2, 0.2, len(vals))
    sc = ax2.scatter(vals, rank + jitter,
                     c=feats, cmap='RdBu_r', s=15,
                     alpha=0.6, vmin=feats.min(), vmax=feats.max())

ax2.set_yticks(range(top_n))
ax2.set_yticklabels(top_names, fontsize=10)
ax2.axvline(x=0, color='black', linewidth=0.8)
ax2.set_xlabel('SHAP value (impact on HBV prediction)', fontsize=11)
ax2.set_title('SHAP Distribution per Feature\n(Red=high value, Blue=low value)',
              fontsize=11, fontweight='bold')
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.grid(axis='x', alpha=0.3)

cbar = plt.colorbar(sc, ax=ax2, shrink=0.6, pad=0.02)
cbar.set_label('Feature value\n(scaled)', fontsize=9)

plt.suptitle('SHAP Analysis — Random Forest Classifier\n'
             'Cytokine + miRNA Biomarker Panel (n=77)',
             fontsize=12, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig('FigureML_E_SHAP.png', dpi=300, bbox_inches='tight')
plt.savefig('FigureML_E_SHAP.pdf', dpi=300, bbox_inches='tight')
plt.show()
print("Figure E (SHAP) saved.")


# ============================================================
# CELL 13 — FIGURE F: XGBoost Feature Importance
#            (for comparison — paper can show both)
# ============================================================

xgb_full = Pipeline([
    ('scaler', StandardScaler()),
    ('clf', XGBClassifier(
        n_estimators=300, max_depth=4,
        learning_rate=0.05, subsample=0.8,
        colsample_bytree=0.8, eval_metric='logloss',
        use_label_encoder=False, random_state=42, verbosity=0
    ))
])
xgb_full.fit(X_main, y_main)

xgb_importance = pd.Series(
    xgb_full.named_steps['clf'].feature_importances_,
    index=feature_names_main
).sort_values(ascending=False)

explainer_xgb  = shap.TreeExplainer(xgb_full.named_steps['clf'])
shap_vals_xgb  = explainer_xgb.shap_values(X_scaled)
if isinstance(shap_vals_xgb, list):
    sv_xgb = shap_vals_xgb[1]
else:
    sv_xgb = shap_vals_xgb

mean_shap_xgb = pd.DataFrame(sv_xgb, columns=feature_names_main
                              ).abs().mean().sort_values(ascending=False)

fig, axes = plt.subplots(1, 2, figsize=(13, 5))

# XGB Gain importance
ax = axes[0]
top_xgb = xgb_importance.head(12)
colors_x = ['#D62728' if f in ['IL-10','IFN-γ','TNF-α','IL-6','TGF-β']
             else '#FF7F0E' for f in top_xgb.index]
ax.barh(top_xgb.index[::-1], top_xgb.values[::-1],
        color=colors_x[::-1], alpha=0.85, edgecolor='white')
ax.set_xlabel('Feature Importance (Gain)', fontsize=11)
ax.set_title('XGBoost Feature Importance', fontsize=11, fontweight='bold')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.grid(axis='x', alpha=0.3)

# XGB SHAP
ax2 = axes[1]
top_xgb_shap = mean_shap_xgb.head(12)
colors_xs = ['#D62728' if f in ['IL-10','IFN-γ','TNF-α','IL-6','TGF-β']
              else '#FF7F0E' for f in top_xgb_shap.index]
ax2.barh(top_xgb_shap.index[::-1], top_xgb_shap.values[::-1],
         color=colors_xs[::-1], alpha=0.85, edgecolor='white')
ax2.set_xlabel('Mean |SHAP value|', fontsize=11)
ax2.set_title('XGBoost — Mean |SHAP|', fontsize=11, fontweight='bold')
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.grid(axis='x', alpha=0.3)

plt.suptitle('XGBoost Feature Analysis — Cytokine + miRNA Panel',
             fontsize=12, fontweight='bold')
plt.tight_layout()
plt.savefig('FigureML_F_XGB_importance.png', dpi=300)
plt.savefig('FigureML_F_XGB_importance.pdf', dpi=300)
plt.show()
print("Figure F (XGBoost importance) saved.")


# ============================================================
# CELL 14 — Summary Table (paper-ready)
# ============================================================

print("\n" + "=" * 70)
print("FINAL RESULTS TABLE (for paper Methods/Results section)")
print("=" * 70)

summary = []
for name, res in cv_results.items():
    oof_pred = oof_preds[name]
    cm = confusion_matrix(y_main, oof_pred)
    tn, fp, fn, tp = cm.ravel()
    sens = tp / (tp + fn)
    spec = tn / (tn + fp)
    ppv  = tp / (tp + fp) if (tp+fp) > 0 else 0
    npv  = tn / (tn + fn) if (tn+fn) > 0 else 0

    summary.append({
        'Model':       name,
        'AUC':         f"{res['AUC']:.3f}",
        '95% CI':      f"[{res['CI_low']:.3f}–{res['CI_up']:.3f}]",
        'Macro F1':    f"{res['F1']:.3f}",
        'Accuracy':    f"{res['Acc']:.3f}",
        'Sensitivity': f"{sens:.3f}",
        'Specificity': f"{spec:.3f}",
        'PPV':         f"{ppv:.3f}",
        'NPV':         f"{npv:.3f}",
    })

summary_df = pd.DataFrame(summary)
print(summary_df.to_string(index=False))
summary_df.to_csv('ML_summary_table_final.csv', index=False)

print("\n✓ All ML analysis complete.")
print("Files saved:")
print("  ML_results_main.csv")
print("  ML_feature_comparison.csv")
print("  ML_summary_table_final.csv")
print("  SHAP_feature_importance.csv")
print("  FigureML_A_ROC_curves.png/pdf")
print("  FigureML_B_Confusion_matrices.png/pdf")
print("  FigureML_C_Model_comparison.png/pdf")
print("  FigureML_D_Feature_comparison.png/pdf")
print("  FigureML_E_SHAP.png/pdf")
print("  FigureML_F_XGB_importance.png/pdf")
