"""
PBMC3k Dataset Downloader
==========================
Downloads and preprocesses the PBMC3k dataset via the Scanpy built-in loader,
which fetches the Seurat-tutorial-processed version with pre-annotated Louvain
cluster labels (CD4 T cells, CD14+ Monocytes, etc.).

Requires: pip install scanpy
"""
from __future__ import annotations

import logging
from pathlib import Path

from benchmarking.download.datasets.base_dataset_downloader import BaseDatasetDownloader

logger = logging.getLogger(__name__)


class PBMC3kDownloader(BaseDatasetDownloader):
    """
    Downloads PBMC3k (Zheng et al. 2017) via scanpy.datasets.pbmc3k_processed().

    The processed dataset includes:
    - 2,638 cells × 1,838 genes (highly variable)
    - adata.obs['louvain']: Seurat tutorial cell type labels (ground truth)
    - Log-normalized raw counts rebuilt from pbmc3k() for DEG + CellTypist

    Output file: <output_dir>/pbmc3k_processed.h5ad
    """

    @property
    def dataset_name(self) -> str:
        return "PBMC3k"

    def _target_path(self) -> Path:
        return self.output_dir / "pbmc3k_processed.h5ad"

    def _download(self) -> Path:
        try:
            import scanpy as sc
        except ImportError as exc:
            raise RuntimeError(
                "scanpy is required to download PBMC3k. "
                "Install with: pip install scanpy"
            ) from exc

        logger.info("Fetching PBMC3k processed (Seurat tutorial labels) via scanpy...")
        adata = sc.datasets.pbmc3k_processed()

        logger.info("Rebuilding log-normalized layer from raw counts...")
        adata_raw = sc.datasets.pbmc3k()
        sc.pp.filter_cells(adata_raw, min_genes=200)
        sc.pp.filter_genes(adata_raw, min_cells=3)
        adata_raw.var["mt"] = adata_raw.var_names.str.startswith("MT-")
        sc.pp.calculate_qc_metrics(
            adata_raw, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
        )
        adata_raw = adata_raw[adata_raw.obs.pct_counts_mt < 5].copy()
        sc.pp.normalize_total(adata_raw, target_sum=1e4)
        sc.pp.log1p(adata_raw)

        # Transfer Louvain labels from processed to log-normalized
        common_cells = adata.obs_names.intersection(adata_raw.obs_names)
        adata_out = adata_raw[common_cells].copy()
        adata_out.obs["louvain"] = adata.obs.loc[common_cells, "louvain"]
        adata_out.raw = adata_out.copy()

        target = self._target_path()
        adata_out.write_h5ad(target)
        logger.info(f"Saved PBMC3k to {target} ({len(adata_out)} cells)")
        return target
