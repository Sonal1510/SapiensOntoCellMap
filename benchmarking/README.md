# Benchmarking Suite

## PBMC3k Benchmark (`benchmark_pbmc3k.py`)

Standard benchmark: Zheng et al. 2017, 2,638 PBMCs, 9 clusters, 8 cell types.

### Quick start (full pipeline)

```bash
# Install optional deps
pip install scanpy celltypist

# Run
source SapiensOntoCellMap_env/bin/activate
python benchmarking/benchmark_pbmc3k.py
```

### If you already have Seurat FindAllMarkers CSV

```bash
python benchmarking/benchmark_pbmc3k.py \
  --deg_csv path/to/seurat_markers.csv \
  --deg_format seurat
```

### Skip CellTypist

```bash
python benchmarking/benchmark_pbmc3k.py --no_celltypist
```

### Outputs (in `benchmarking/results/`)

| File | Description |
|------|-------------|
| `pbmc3k_degs.csv` | Scanpy-format DEGs |
| `sapiensonto_out/` | Full SapiensOntoCellMap output (CSV + HTML) |
| `celltypist_results.csv` | CellTypist cluster-level labels |
| `benchmark_results.csv` | Per-cluster metrics (all tools) |
| `benchmark_summary.txt` | Aggregate accuracy scores |
| `benchmark_per_cluster_heatmap.png` | Figure 1: accuracy heatmap |
| `benchmark_tool_comparison.png` | Figure 2: bar chart comparison |
| `benchmark_confidence_vs_accuracy.png` | Figure 3: confidence vs. accuracy |

### Metrics

| Metric | Description |
|--------|-------------|
| Top-1 accuracy | Top predicted cell type string-matches ground truth |
| CL exact accuracy | Predicted CL ID == ground truth CL ID |
| Hierarchical accuracy | Predicted cell type is in the correct broad lineage |
| Broad-type accuracy | `Broad_Type` column matches ground truth lineage |
| Concordance | % clusters where SapiensOntoCellMap and CellTypist agree |

### Adding scType (R) results

1. Run scType in R on the PBMC3k dataset
2. Export cluster → cell type mapping to CSV
3. Fill in `MANUAL_COMPARATORS` in `benchmark_pbmc3k.py`:

```python
MANUAL_COMPARATORS = {
    "scType": {
        "0": "CD4+ T cells",
        "1": "CD14+ Monocytes",
        # ...
    }
}
```

### Extending to other datasets

- **Human Lung Cell Atlas (HLCA):** subset to ~5K cells, run same pipeline
- **In-house skin scRNA-seq:** use `run_scrnaseq_samples.sh` + define custom `GROUND_TRUTH` dict
