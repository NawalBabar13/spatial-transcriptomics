# Tutorial 1 — Basic Scanpy Spatial Analysis

**Library:** `scanpy`  
**Dataset:** 10x Visium — V1 Mouse Brain Sagittal Posterior  
**Reference:** [Scanpy Spatial Tutorial](https://scanpy-tutorials.readthedocs.io/en/latest/spatial/basic-analysis.html)

---

## 🎯 Objective

This tutorial introduces the foundational workflow for analyzing **10x Visium spatial transcriptomics data** using Scanpy. We load a mouse brain tissue section, assess data quality, cluster spots by transcriptional similarity, and project those clusters back onto the physical tissue image to reveal spatial gene expression patterns in the brain.

---

## 🧠 Biological Context

The dataset is a sagittal section of the **mouse brain posterior** captured using the 10x Visium platform. Each spot on the Visium slide is ~55 µm in diameter and captures RNA from approximately 1–10 cells directly beneath it. By clustering spots with similar gene expression profiles and then mapping them back to their spatial coordinates, we can identify distinct anatomical regions of the brain — such as the cortex, hippocampus, and cerebellum — without any prior anatomical annotation.

---

## 📂 Files in This Folder

| File | Description |
|------|-------------|
| `Tutorial1_Basic_Scanpy_Spatial.ipynb` | Full runnable notebook |
| `QC_Plots.png` | Histograms of total counts and genes per spot |
| `UMAP Plot (transcriptional similarity).png` | UMAP embedding colored by Leiden cluster |
| `spatial plot_1.png` | Leiden clusters overlaid on H&E tissue image |
| `spatial plot_2.png` | Specific gene expression overlaid on tissue |
| `rank_genes_groups_heatmap.png` | Top marker genes per cluster (heatmap) |

---

## 🔬 Methods & Pipeline

### Step 1 — Data Loading
```python
adata = sc.datasets.visium_sge(sample_id="V1_Mouse_Brain_Sagittal_Posterior")
adata.var_names_make_unique()
```
The `AnnData` object (`adata`) stores:
- `.X` — gene expression matrix (spots × genes)
- `.obs` — spot-level metadata (coordinates, QC metrics)
- `.var` — gene-level metadata
- `.uns['spatial']` — the high-resolution tissue image

### Step 2 — Quality Control
```python
sc.pp.calculate_qc_metrics(adata, inplace=True)
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.filter_cells(adata, min_genes=200)
```
We remove:
- **Genes** detected in fewer than 3 spots (likely noise or dropouts)
- **Spots** with fewer than 200 genes detected (likely empty or damaged spots)

**Output:** `QC_Plots.png` — histograms showing the distribution of total UMI counts and number of genes per spot before filtering.

### Step 3 — Normalization
```python
sc.pp.normalize_total(adata, inplace=True)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000)
```
- **Library-size normalization:** scales each spot so total counts = 10,000, correcting for sequencing depth differences
- **Log1p transformation:** log(count + 1), stabilizes variance and makes the data more normally distributed
- **Highly variable genes:** selects the 2,000 most informative genes for downstream analysis, reducing noise

### Step 4 — Dimensionality Reduction
```python
sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
```
- **PCA:** reduces the 2,000-gene matrix to 50 principal components capturing the main axes of variation
- **Neighbors graph:** builds a k-nearest-neighbor graph in PCA space (spots that are transcriptionally similar are connected)
- **UMAP:** further reduces to 2D for visualization while preserving local structure

### Step 5 — Leiden Clustering
```python
sc.tl.leiden(adata, key_added="clusters")
```
Leiden is a graph-based community detection algorithm. It partitions the neighbor graph into groups (clusters) of transcriptionally similar spots. These clusters often correspond to distinct cell types or tissue regions.

**Output:** `UMAP Plot (transcriptional similarity).png`

### Step 6 — Spatial Visualization
```python
sc.pl.spatial(adata, img_key="hires", color="clusters")
sc.pl.spatial(adata, img_key="hires", color=["Nrgn", "Camk2a"])
```
This projects the cluster labels back onto the physical tissue image, revealing the spatial organization of transcriptionally distinct populations.

**Outputs:** `spatial plot_1.png`, `spatial plot_2.png`

### Step 7 — Marker Gene Identification
```python
sc.tl.rank_genes_groups(adata, "clusters", method="t-test")
sc.pl.rank_genes_groups_heatmap(adata, n_genes=5, groupby="clusters")
```
For each cluster, we identify genes that are significantly more highly expressed compared to all other clusters (one-vs-rest t-test). The heatmap shows the top 5 marker genes per cluster.

**Output:** `rank_genes_groups_heatmap.png`

---

## 📊 Results

| Step | Output |
|------|--------|
| QC | Spots with <200 genes removed; genes in <3 spots removed |
| Normalization | Library-size normalized + log1p transformed |
| Clustering | Leiden algorithm identified distinct transcriptional clusters |
| Spatial mapping | Clusters map to anatomically coherent brain regions |
| Marker genes | Top differentially expressed genes identified per cluster |

---

## 💡 Key Takeaways

- Visium data combines gene expression with spatial coordinates, enabling anatomically-aware analysis
- Spots from the same brain region cluster together in UMAP space, confirming that transcriptional similarity reflects spatial proximity
- Marker genes for each cluster correspond to known brain region-specific genes (e.g. *Nrgn* for cortex, *Camk2a* for excitatory neurons)
- The spatial plot reveals clear regional structure — cortex, hippocampus, and subcortical regions are visually separable
