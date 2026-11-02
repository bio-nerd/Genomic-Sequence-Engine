import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA

sns.set_theme(style="ticks", palette="muted")
print("Generating synthetic TCGA expression data...")

np.random.seed(42)
n_samples_per_class, n_genes = 50, 500
classes = ['BRCA', 'LUAD', 'PRAD', 'KIRC']

data, labels = [], []
for i, label in enumerate(classes):
    gene_means = np.random.normal(loc=i * 2.0, scale=1.0, size=n_genes)
    samples = np.random.normal(loc=gene_means, scale=1.5, size=(n_samples_per_class, n_genes))
    data.append(samples)
    labels.extend([label] * n_samples_per_class)

X = np.vstack(data)
y = np.array(labels)

pca = PCA(n_components=2)
pcs = pca.fit_transform(X)
pca_df = pd.DataFrame(pcs, columns=['PC1', 'PC2'])
pca_df['Cancer_Type'] = y

plt.figure(figsize=(8, 6), dpi=300)
sns.scatterplot(x='PC1', y='PC2', hue='Cancer_Type', style='Cancer_Type', data=pca_df, s=70, alpha=0.85)
plt.title("TCGA Pan-Cancer Gene Expression Clustering (PCA)", fontsize=12, fontweight='bold')
plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}% Variance)")
plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}% Variance)")
plt.legend(title="Subtype", bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()

plot_path = "tcga_pca_plot.png"
plt.savefig(plot_path, dpi=300)
print(f"Success! High-res plot saved as '{plot_path}'.")
