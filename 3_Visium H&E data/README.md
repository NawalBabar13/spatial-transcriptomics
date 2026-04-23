# Tutorial 3 — Visium H&E Data

**Library:** `squidpy`  
**Dataset:** Squidpy built-in Visium H&E mouse brain dataset  
**Reference:** [Squidpy Visium H&E Tutorial](https://squidpy.readthedocs.io/en/stable/notebooks/tutorials/tutorial_visium_hne.html)

---

## 🎯 Objective

This tutorial analyzes a **Visium dataset paired with an H&E (Hematoxylin and Eosin) stained tissue image** of the mouse brain. Beyond clustering, we focus on two key spatial analyses: **Moran's I spatial autocorrelation** (which genes are spatially patterned?) and **neighborhood enrichment** (which clusters physically neighbor each other?).

---

## 🩺 What Is H&E Staining?

H&E is the most widely used stain in histopathology:
- **Hematoxylin** — stains cell nuclei blue/purple (binds to DNA)
- **Eosin** — stains cytoplasm and extracellular matrix pink

The result is the classic pink-and-purple tissue image familiar from pathology slides. In the context of spatial transcriptomics, the H&E image provides an anatomical reference that can be visually compared to gene expression clusters — allowing you to verify that transcriptomically-defined clusters correspond to histologically-distinct tissue regions.

---

## 📂 Files in This Folder

| File | Description |
|------|-------------|
| `Tutorial3_Visium_HnE.ipynb` | Full runnable notebook |
| `View the H&E Tissue Image.png` | Raw H&E image with all spots plotted |
| `Visualize Clusters on H&E Image.png` | Leiden clusters overlaid on H&E tissue |
| `spatial_scatter.png` | Spatially variable gene expression on tissue |
| `Neighborhood Enrichment Analysis.png` | Heatmap of cluster co-localization scores |

---

## 🔬 Methods & Pipeline

### Step 1 — Data Loading
```python
adata = sq.datasets.visium_hne_adata()
```
Dataset contains:
- Mouse brain Visium gene expression (spots × genes)
- High-resolution H&E image
- Spot coordinates registered to the image

### Step 2 — View Tissue Image
```python
sq.pl.spatial_scatter(adata, shape="circle", size=0.5)
```
Renders the H&E tissue image with all capture spots overlaid. This confirms the spatial layout of the barcoded spots across the tissue section.

**Output:** `View the H&E Tissue Image.png`

### Step 3 — QC, Normalization, Clustering
Standard pipeline:
```python
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.normalize_total(adata, inplace=True)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000)
sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.leiden(adata, key_added="clusters")
```

### Step 4 — Spatial Cluster Visualization
```python
sq.pl.spatial_scatter(adata, color="clusters")
```
Overlays cluster identity onto each spot on the H&E image. Brain anatomical structures become clearly visible — the cortex, hippocampus, white matter tracts, and subcortical regions each contain distinct clusters.

**Output:** `Visualize Clusters on H&E Image.png`

### Step 5 — Spatial Neighborhood Graph
```python
sq.gr.spatial_neighbors(adata)
```
Constructs a graph where nodes are Visium spots and edges connect physically adjacent spots. This is distinct from the transcriptomic neighbor graph (used for clustering) — here adjacency is purely geometric.

### Step 6 — Moran's I Spatial Autocorrelation
```python
sq.gr.spatial_autocorr(adata, mode="moran")
adata.uns["moranI"].head(10)
```

**Moran's I** is a classical spatial statistics metric. For each gene, it computes a score between -1 and +1:

| Moran's I value | Interpretation |
|-----------------|----------------|
| Close to **+1** | Gene is highly spatially clustered (expressed in specific regions) |
| Close to **0** | Gene is randomly distributed across the tissue |
| Close to **-1** | Gene is spatially dispersed (checkerboard-like pattern) |

Genes with high Moran's I are called **spatially variable genes (SVGs)** — they carry the most spatial information and are biologically interesting because their expression is not random but tied to tissue anatomy.

```python
top_genes = adata.uns["moranI"].head(4).index.tolist()
sq.pl.spatial_scatter(adata, color=top_genes, ncols=2)
```

**Output:** `spatial_scatter.png`

### Step 7 — Neighborhood Enrichment Analysis
```python
sq.gr.nhood_enrichment(adata, cluster_key="clusters")
sq.pl.nhood_enrichment(adata, cluster_key="clusters")
```

Neighborhood enrichment tests whether pairs of clusters are **physically co-localized** more than expected by chance. For each pair of clusters (A, B), it counts how many spots of cluster B are neighbors of spots in cluster A, then compares to a permutation null distribution.

- **Positive enrichment score** → clusters A and B tend to be physically adjacent (e.g. cortex layer I sits next to layer II)
- **Negative enrichment score** → clusters A and B are spatially segregated

The result is a heatmap showing pairwise enrichment scores for all cluster combinations.

**Output:** `Neighborhood Enrichment Analysis.png`

---

## 📊 Results

| Analysis | Finding |
|----------|---------|
| H&E visualization | Brain anatomy visible — distinct regions correspond to Leiden clusters |
| Moran's I | Top SVGs show clear spatial zonation — strongly expressed in specific anatomical regions |
| Neighborhood enrichment | Adjacent cortical layers show mutual enrichment; cortex and cerebellum clusters show negative enrichment (spatially separated) |

---

## 💡 Key Takeaways

- H&E staining provides an anatomical ground truth that spatially-defined transcriptomic clusters can be directly compared against
- Moran's I is more rigorous than simply looking at spatial plots — it provides a ranked, statistical list of which genes carry spatial information
- Neighborhood enrichment reveals tissue architecture at the cluster level: which cell populations form boundaries, which intermingle, and which are segregated
- Together, spatially variable genes + neighborhood enrichment give a complete picture of tissue organization from both a molecular and structural perspective
