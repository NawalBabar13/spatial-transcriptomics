# Tutorial 4 — Xenium Single-Cell Spatial Data

**Library:** `squidpy`  
**Dataset:** Squidpy built-in Xenium mouse brain dataset  
**Reference:** [Squidpy Xenium Tutorial](https://squidpy.readthedocs.io/en/stable/notebooks/tutorials/tutorial_xenium.html)

---

## 🎯 Objective

This tutorial analyzes data from the **10x Genomics Xenium platform**, which achieves **single-cell spatial resolution**. Unlike Visium (which captures RNA from ~55 µm spots containing multiple cells), Xenium directly images individual RNA molecules inside individual cells, giving you the precise location of every cell in the tissue along with its gene expression profile.

---

## 🔭 Xenium vs Visium — Key Differences

| Feature | Visium | Xenium |
|---------|--------|--------|
| Resolution | ~55 µm spots | Single cell (~10 µm) |
| Cells per measurement | 1–10 cells/spot | 1 cell per barcode |
| Gene coverage | Whole transcriptome | Targeted panel (~300–500 genes) |
| Image type | H&E or fluorescence | DAPI + fluorescence |
| Spatial coordinates | Grid layout (regular) | Arbitrary cell centroids |
| Data size | Thousands of spots | Tens of thousands of cells |

Xenium trades whole-transcriptome coverage for dramatically higher spatial resolution, making it ideal for studying cell-type organization, cell-cell communication, and tissue microenvironments.

---

## 📂 Files in This Folder

| File | Description |
|------|-------------|
| `Tutorial4_Xenium.ipynb` | Full runnable notebook |
| `spatial_scatter_ Raw Spatial Layout .png` | All single cells plotted at their tissue coordinates |
| `UMAP of Single Cells.png` | UMAP embedding of single cells colored by Leiden cluster |
| `Spatially Variable Genes_1.png` | Top spatially variable genes (Moran's I) — set 1 |
| `Spatially Variable Genes_2.png` | Top spatially variable genes (Moran's I) — set 2 |
| `Neighborhood Enrichment_heatmap.png` | Pairwise cell-type co-localization heatmap |

---

## 🔬 Methods & Pipeline

### Step 1 — Data Loading
```python
adata = sq.datasets.xenium()
print(f"Cells: {adata.n_obs}, Genes: {adata.n_vars}")
```
The AnnData object stores:
- `.X` — expression matrix (cells × targeted genes)
- `.obsm['spatial']` — x/y centroid coordinates for each cell
- `.uns['spatial']` — DAPI/fluorescence image of the tissue

### Step 2 — Raw Spatial Layout
```python
sq.pl.spatial_scatter(adata, shape=None, size=1)
```
Note: we use `shape=None` and `size=1` because cells are individual points (not large 55 µm Visium spots). This plot shows thousands of individual cells scattered across the tissue section.

**Output:** `spatial_scatter_ Raw Spatial Layout .png`

### Step 3 — QC and Normalization
```python
sc.pp.filter_genes(adata, min_cells=5)
sc.pp.filter_cells(adata, min_genes=5)
sc.pp.normalize_total(adata, inplace=True)
sc.pp.log1p(adata)
```
Stricter filtering thresholds than Visium because Xenium targets a small gene panel — a cell with fewer than 5 detected genes from a targeted panel likely has low RNA capture quality.

### Step 4 — PCA, UMAP, Leiden Clustering
```python
sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.leiden(adata, key_added="clusters")
```
Same dimensionality reduction workflow — clusters now represent **cell types** (individual cells, not tissue spots).

**Output:** `UMAP of Single Cells.png`

### Step 5 — Mapping Clusters Back to Space
```python
sq.pl.spatial_scatter(adata, color="clusters", shape=None, size=1)
```
Each dot on the tissue image is a single cell, colored by its cluster identity. At this resolution you can see:
- Laminar organization of cortical layers
- Sharp boundaries between brain regions
- Small cellular niches invisible at Visium resolution

### Step 6 — Spatial Neighborhood Graph (Delaunay)
```python
sq.gr.spatial_neighbors(adata, coord_type="generic", delaunay=True)
```
For Xenium we use:
- `coord_type="generic"` — because cell coordinates are arbitrary (not a regular hexagonal grid like Visium)
- `delaunay=True` — uses **Delaunay triangulation** to connect each cell to its geometrically natural neighbors. This is more appropriate than a fixed radius for irregularly distributed single cells.

### Step 7 — Moran's I Spatial Autocorrelation
```python
sq.gr.spatial_autocorr(adata, mode="moran")
top_genes = adata.uns["moranI"].head(4).index.tolist()
sq.pl.spatial_scatter(adata, color=top_genes, shape=None, size=1, ncols=2)
```
At single-cell resolution, Moran's I identifies genes whose expression is spatially structured. The results are far more refined than at Visium spot resolution — you can see expression patterns at the scale of individual cell layers.

**Outputs:** `Spatially Variable Genes_1.png`, `Spatially Variable Genes_2.png`

### Step 8 — Neighborhood Enrichment
```python
sq.gr.nhood_enrichment(adata, cluster_key="clusters")
sq.pl.nhood_enrichment(adata, cluster_key="clusters")
```
At single-cell resolution, neighborhood enrichment reveals **cell-type co-localization** patterns:
- Which cell types form physical boundaries with each other?
- Which cell types are intermingled?
- Are there unexpected spatial relationships between cell populations?

**Output:** `Neighborhood Enrichment_heatmap.png`

---

## 📊 Results

| Analysis | Finding |
|----------|---------|
| Raw spatial layout | Tens of thousands of individual cells visible, denser in grey matter regions |
| UMAP clustering | Distinct cell type populations visible as separate UMAP lobes |
| Spatial clusters | Cortical layers, hippocampal subfields, and white matter tracts resolved at single-cell level |
| Moran's I | SVGs show sharp, layer-specific expression gradients at cellular resolution |
| Neighborhood enrichment | Cortical neurons of adjacent layers show high mutual enrichment; inhibitory and excitatory neurons show distinct co-localization patterns |

---

## 💡 Key Takeaways

- Xenium's single-cell resolution resolves tissue architecture that is invisible at Visium resolution — individual cortical layers, microenvironmental niches, and rare cell populations
- Delaunay triangulation is the appropriate neighborhood construction method for cells with arbitrary spatial positions
- Moran's I at single-cell resolution reveals gene expression gradients that mirror histological boundaries
- Neighborhood enrichment at cell-type level provides a quantitative map of tissue microarchitecture — essentially a molecular description of cell-cell proximity
- The targeted gene panel limits global transcriptomic insight but enables much higher spatial resolution and throughput

---

## 🔁 Comparison: Tutorials 3 vs 4

| Aspect | Tutorial 3 (Visium H&E) | Tutorial 4 (Xenium) |
|--------|------------------------|---------------------|
| Resolution | Spots (multi-cell) | Single cells |
| Gene coverage | Whole transcriptome | Targeted ~300–500 genes |
| Spatial graph | Hexagonal grid | Delaunay triangulation |
| Cluster interpretation | Tissue region types | Cell types |
| Moran's I resolution | Regional patterns | Cell-layer precision |
