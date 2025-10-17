# Inputs (set these paths):
NORMAL_PATH = "/Users/yeseniajuarez/Documents/genomic_analysis/Homework_3/normal_merged.csv"
TUMOR_PATH  = "/Users/yeseniajuarez/Documents/genomic_analysis/Homework_3/tumor_merged.csv"

import pandas as pd
import numpy as np

#read file inputs 
def read_inputs(normal_path, tumor_path):
    normal = pd.read_csv(normal_path)
    tumor  = pd.read_csv(tumor_path)
    for df in (normal, tumor):
        if "alt_seq" not in df.columns and "var_seq2" in df.columns:
            df["alt_seq"] = df["var_seq2"]
        if "left" in df.columns:
            df["left"] = pd.to_numeric(df["left"], errors="coerce")
    return normal, tumor

def group(df, tag):
    keep_cols = ["chrom","left","ref_seq","alt_seq","Patient_ID","VCF_ID"]
    sub = df[keep_cols].copy()
    grouped = (sub
               .groupby(["chrom","left","ref_seq","alt_seq"], dropna=False)
               .agg({"Patient_ID": lambda x: sorted(pd.unique(x)),
                     "VCF_ID":    lambda x: sorted(pd.unique(x))})
               .reset_index())
    grouped[f"{tag}#"] = grouped["Patient_ID"].apply(lambda lst: len(set(lst)))
    grouped = grouped.rename(columns={"Patient_ID": f"Patient_ID_{tag}",
                                      "VCF_ID":    f"VCF_ID_{tag}"})
    unique_patients = sub["Patient_ID"].nunique()
    return grouped, unique_patients

def merge(normal_grouped, tumor_grouped):
    keys = ["chrom","left","ref_seq","alt_seq"]
    aml = pd.merge(normal_grouped, tumor_grouped, on=keys, how="outer")
    unique_normal_variants = len(normal_grouped)
    unique_tumor_variants  = len(tumor_grouped)
    shared = pd.merge(normal_grouped[keys], tumor_grouped[keys], on=keys, how="inner").drop_duplicates()
    shared_variants = len(shared)
    return aml, unique_normal_variants, unique_tumor_variants, shared_variants

def parse_list(x):
    if isinstance(x, list):
        return x
    elif pd.isna(x):
        return []
    else:
        return [x]


def aml_expand(normal, tumor):

    all_cols = sorted(set(normal.columns).union(tumor.columns))
    normal_r = normal.reindex(columns=all_cols)
    tumor_r  = tumor.reindex(columns=all_cols)

    both = pd.concat([normal_r, tumor_r], axis=0, ignore_index=True)

    expand_cols = [c for c in ["SYMBOL","Gene","Feature","cDNA_position","BIOTYPE","Consequence"] if c in both.columns]
    parsed = {c: both[c].apply(parse_list) for c in expand_cols}
    lens   = pd.DataFrame({c: parsed[c].str.len() for c in expand_cols}) if expand_cols else pd.DataFrame(index=both.index)
    maxlen = (lens.max(axis=1) if not lens.empty else pd.Series(1, index=both.index)).fillna(1).astype(int)

    expanded = both.loc[both.index.repeat(maxlen)].copy()
    expanded.reset_index(drop=True, inplace=True)
    expanded["orig_idx"] = both.index.repeat(maxlen).to_numpy()
    expanded["pos"] = expanded.groupby("orig_idx").cumcount()

    for c in expand_cols:
        expanded[c] = expanded.apply(
            lambda r: parsed[c].iloc[r["orig_idx"]][r["pos"]]
            if r["pos"] < len(parsed[c].iloc[r["orig_idx"]]) else np.nan, axis=1
        )

    expanded = expanded.drop(columns=["orig_idx","pos"]).drop_duplicates().reset_index(drop=True)
    return expanded



normal, tumor = read_inputs(NORMAL_PATH, TUMOR_PATH)

#  subset, group to lists, add N#/T#, rename
normal_grouped, n_norm_patients = group(normal, "Normal")
tumor_grouped,  n_tum_patients  = group(tumor,  "Tumor")
normal_grouped.to_csv("AML_normal_grouped.csv", index=False)
tumor_grouped.to_csv("AML_tumor_grouped.csv", index=False)
print("1.1.1 Unique normal patients:", n_norm_patients)
print("1.1.2 Unique tumor patients:",  n_tum_patients)

# outer merge on variant keys
aml, u_norm_vars, u_tum_vars, shared_vars = merge(normal_grouped, tumor_grouped)
aml.to_csv("AML.csv", index=False)
print("1.2.1 Unique normal variants:", u_norm_vars)
print("1.2.2 Unique tumor variants:",  u_tum_vars)
print("1.2.3 Shared variants (common):", shared_vars)

# concat + expand CSQ lists 
expanded = aml_expand(normal, tumor)
expanded.to_csv("AML_Expand.csv", index=False)
print("1.3.1 Rows in AML_Expand.csv:", len(expanded))

gene_cols = [c for c in ["SYMBOL","Gene","Feature"] if c in expanded.columns]
tx_cols   = [c for c in ["chrom","left","right","ref_seq","alt_seq","Feature","cDNA_position","BIOTYPE"] if c in expanded.columns]
expanded[gene_cols].drop_duplicates().to_csv("AML_gene.csv", index=False)
expanded[tx_cols].drop_duplicates().to_csv("AML_tx.csv", index=False)

#Part 2: Random Forest
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sklearn

from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import accuracy_score, confusion_matrix
from sklearn.ensemble import RandomForestClassifier
from sklearn.datasets import load_iris

iris = load_iris()
df = pd.DataFrame(data=iris.data, columns=iris.feature_names)
df['target'] = iris.target

print(df)

X = df.iloc[:, :-1].values
y = df.iloc[:, -1].values

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

scaler = StandardScaler()
X_train = scaler.fit_transform(X_train)
X_test = scaler.transform(X_test)

classifier = RandomForestClassifier(n_estimators=10, max_depth=4, random_state=42)
classifier.fit(X_train, y_train)
y_pred = classifier.predict(X_test)

accuracy = accuracy_score(y_test, y_pred)
print(f'Accuracy: {accuracy * 100:.2f}%')

conf_matrix = confusion_matrix(y_test, y_pred)

for name, importance in zip(iris.feature_names, classifier.feature_importances_):
    print(name, importance)

plt.figure(figsize=(8, 6))
sns.heatmap(conf_matrix, annot=True, fmt='g', cmap='Blues', cbar=False, 
            xticklabels=iris.target_names, yticklabels=iris.target_names)

plt.title('Confusion Matrix Heatmap')
plt.xlabel('Predicted Labels')
plt.ylabel('True Labels')
plt.show()



#Part 3: K-Means Clustering
import numpy as np

from sklearn.datasets import load_digits

data, labels = load_digits(return_X_y=True)
(n_samples, n_features), n_digits = data.shape, np.unique(labels).size

print(f"# digits: {n_digits}; # samples: {n_samples}; # features {n_features}")

from time import time

from sklearn import metrics
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler


def bench_k_means(kmeans, name, data, labels):
    """Benchmark to evaluate the KMeans initialization methods.

    Parameters
    ----------
    kmeans : KMeans instance
        A :class:`~sklearn.cluster.KMeans` instance with the initialization
        already set.
    name : str
        Name given to the strategy. It will be used to show the results in a
        table.
    data : ndarray of shape (n_samples, n_features)
        The data to cluster.
    labels : ndarray of shape (n_samples,)
        The labels used to compute the clustering metrics which requires some
        supervision.
    """
    t0 = time()
    estimator = make_pipeline(StandardScaler(), kmeans).fit(data)
    fit_time = time() - t0
    results = [name, fit_time, estimator[-1].inertia_]

    # Define the metrics which require only the true labels and estimator
    # labels
    clustering_metrics = [
        metrics.homogeneity_score,
        metrics.completeness_score,
        metrics.v_measure_score,
        metrics.adjusted_rand_score,
        metrics.adjusted_mutual_info_score,
    ]
    results += [m(labels, estimator[-1].labels_) for m in clustering_metrics]

    # The silhouette score requires the full dataset
    results += [
        metrics.silhouette_score(
            data,
            estimator[-1].labels_,
            metric="euclidean",
            sample_size=300,
        )
    ]

    # Show the results
    formatter_result = (
        "{:9s}\t{:.3f}s\t{:.0f}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}"
    )
    print(formatter_result.format(*results))
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA

print(82 * "_")
print("init\t\ttime\tinertia\thomo\tcompl\tv-meas\tARI\tAMI\tsilhouette")

kmeans = KMeans(init="k-means++", n_clusters=n_digits, n_init=4, random_state=0)
bench_k_means(kmeans=kmeans, name="k-means++", data=data, labels=labels)

kmeans = KMeans(init="random", n_clusters=n_digits, n_init=4, random_state=0)
bench_k_means(kmeans=kmeans, name="random", data=data, labels=labels)

pca = PCA(n_components=n_digits).fit(data)
kmeans = KMeans(init=pca.components_, n_clusters=n_digits, n_init=1)
bench_k_means(kmeans=kmeans, name="PCA-based", data=data, labels=labels)

print(82 * "_")

import matplotlib.pyplot as plt

reduced_data = PCA(n_components=2).fit_transform(data)
kmeans = KMeans(init="k-means++", n_clusters=n_digits, n_init=4)
kmeans.fit(reduced_data)

# Step size of the mesh. Decrease to increase the quality of the VQ.
h = 0.02  # point in the mesh [x_min, x_max]x[y_min, y_max].

# Plot the decision boundary. For that, we will assign a color to each
x_min, x_max = reduced_data[:, 0].min() - 1, reduced_data[:, 0].max() + 1
y_min, y_max = reduced_data[:, 1].min() - 1, reduced_data[:, 1].max() + 1
xx, yy = np.meshgrid(np.arange(x_min, x_max, h), np.arange(y_min, y_max, h))

# Obtain labels for each point in mesh. Use last trained model.
Z = kmeans.predict(np.c_[xx.ravel(), yy.ravel()])

# Put the result into a color plot
Z = Z.reshape(xx.shape)
plt.figure(1)
plt.clf()
plt.imshow(
    Z,
    interpolation="nearest",
    extent=(xx.min(), xx.max(), yy.min(), yy.max()),
    cmap=plt.cm.Paired,
    aspect="auto",
    origin="lower",
)

plt.plot(reduced_data[:, 0], reduced_data[:, 1], "k.", markersize=2)
# Plot the centroids as a white X
centroids = kmeans.cluster_centers_
plt.scatter(
    centroids[:, 0],
    centroids[:, 1],
    marker="x",
    s=169,
    linewidths=3,
    color="w",
    zorder=10,
)
plt.title(
    "K-means clustering on the digits dataset (PCA-reduced data)\n"
    "Centroids are marked with white cross"
)
plt.xlim(x_min, x_max)
plt.ylim(y_min, y_max)
plt.xticks(())
plt.yticks(())
plt.show()