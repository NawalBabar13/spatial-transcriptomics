# Tutorial 2 — Visium Fluorescence Data

**Library:** `squidpy`  
**Dataset:** Squidpy built-in Visium fluorescence dataset  
**Reference:** [Squidpy Documentation](https://squidpy.readthedocs.io/en/stable/)

---

## 🎯 Objective

This tutorial extends the basic Scanpy workflow by introducing **Squidpy** — a library built on top of Scanpy specifically for spatial omics analysis. We analyze a Visium dataset paired with a **fluorescence microscopy image** (rather than an H&E stain), extract quantitative image features from the fluorescence channels, and apply spatial statistics to understand how cell populations are distributed across the tissue.

---

## 🔬 What Is a Fluorescence Visium Dataset?

Unlike H&E staining (which uses chemical dyes to color cell nuclei and cytoplasm), fluorescence Visium uses **immunofluorescence** — antibodies tagged with fluorescent dyes that light up specific proteins. This gives multiple image channels, each representing a different protein marker. The advantage is that you can correlate gene expression from the transcriptomic data with protein-level signals from the image.

---

## 📂 Files in This Folder

| File | Description |
|------|-------------|
| `Tutorial2_Visium_Fluorescence.ipynb` | Full runnable notebook |
| `Spacial_scatter_plot.png` | All spots visualized on the fluorescence tissue image |
| `Spatial Scatter Plot of Clusters.png` | Leiden clusters overlaid on tissue |
| `Ripley's Statistics (Spatial Distribution).png` | Ripley's L curves per cluster |

---

## 🔬 Methods & Pipeline

### Step 1 — Data Loading
```python
adata = sq.datasets.visium_fluo_adata()
```
The dataset includes:
- Gene expression matrix (spots × genes)
- Multi-channel fluorescence image stored in `adata.uns['spatial']`
- Spot coordinates mapped to the image

### Step 2 — Initial Spatial Visualization
```python
sq.pl.spatial_scatter(adata, shape="circle", size=0.5)
```
Squidpy's `spatial_scatter` renders spots as circles on the tissue image. This gives an immediate visual of tissue coverage and spot density.

**Output:** `Spacial_scatter_plot.png`

### Step 3 — QC, Normalization, Clustering
Same pipeline as Tutorial 1:
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
**Output:** `Spatial Scatter Plot of Clusters.png`

### Step 5 — Image Feature Extraction
```python
sq.im.calculate_image_features(
    adata,
    sq.im.ImageContainer.from_adata(adata),
    features="summary",
    key_added="img_features"
)
```
Squidpy extracts **quantitative features from the fluorescence image** for each spot's local image patch. Summary features include:
- **Mean intensity** per channel — how bright the fluorescence signal is in that spot
- **Standard deviation** — variability of signal
- **Percentile values** — robust descriptors of signal distribution

These features are stored in `adata.obsm["img_features"]` and can be correlated with gene expression clusters to link transcriptomic and proteomic signals.

### Step 6 — Spatial Neighborhood Graph
```python
sq.gr.spatial_neighbors(adata)
```
Builds a graph where each spot is a node, and edges connect physically adjacent spots. This graph is the foundation for all spatial statistics in Squidpy. By default it uses a hexagonal grid structure (matching the Visium spot layout).

### Step 7 — Ripley's L Statistic
```python
sq.gr.ripley(adata, cluster_key="clusters", mode="L")
sq.pl.ripley(adata, cluster_key="clusters", mode="L")
```
**Ripley's L** is a spatial statistics test that asks: *"Are spots of this cluster more spatially clustered than expected if they were randomly distributed?"*

- The L curve is plotted against a radius `r`
- **L(r) > 0:** spots are more clustered than random at that radius (spatial aggregation)
- **L(r) = 0:** spots are randomly distributed (complete spatial randomness)
- **L(r) < 0:** spots are more dispersed than random (spatial regularity)

**Output:** `Ripley's Statistics (Spatial Distribution).png`

---

## 📊 Results

| Analysis | Finding |
|----------|---------|
| Spatial scatter | Spots uniformly cover the tissue section |
| Clustering | Leiden identifies transcriptionally distinct populations |
| Image features | Fluorescence intensity varies by cluster, linking protein and RNA signals |
| Ripley's L | Certain clusters show significant spatial aggregation (L > 0), indicating they occupy defined tissue regions |

---

## 💡 Key Takeaways

- Squidpy's `spatial_scatter` gives finer control over spatial visualization than Scanpy's `pl.spatial`
- Image feature extraction bridges transcriptomic clusters with the physical appearance of tissue under fluorescence — a multi-modal approach
- Ripley's L provides a rigorous statistical test for spatial clustering, going beyond visual inspection
- Fluorescence images carry additional biological information (protein expression) that can complement and validate transcriptome-based clusters
