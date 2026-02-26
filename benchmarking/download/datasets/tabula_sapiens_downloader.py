"""
Tabula Sapiens Dataset Downloader
===================================
Downloads skin and blood subsets from the Tabula Sapiens atlas via CELLxGENE.

Reference: The Tabula Sapiens Consortium, Science 2022 (doi:10.1126/science.abl4896)
Data access: https://cellxgene.cziscience.com/collections/e5f58829-1a66-40b5-a624-9046778e74f5

Requires: pip install cellxgene-census
"""
from __future__ import annotations

import logging
from pathlib import Path

from benchmarking.download.datasets.base_dataset_downloader import BaseDatasetDownloader

logger = logging.getLogger(__name__)

# CellxGene collection ID for Tabula Sapiens (stable)
TABULA_SAPIENS_COLLECTION_ID = "e5f58829-1a66-40b5-a624-9046778e74f5"

# Tissues to include in the benchmark subset
BENCHMARK_TISSUES = ("blood", "skin of body")

# Maximum cells per cell type to include (avoids class imbalance)
MAX_CELLS_PER_TYPE = 300


class TabulaSapiensDownloader(BaseDatasetDownloader):
    """
    Downloads a benchmark-sized subset of Tabula Sapiens (skin + blood).

    Fetches up to MAX_CELLS_PER_TYPE cells per annotated cell type via the
    CELLxGENE Census API. The full atlas (~500K cells) is NOT downloaded;
    only the targeted subset is fetched to disk.

    Output file: <output_dir>/tabula_sapiens_benchmark.h5ad
    """

    @property
    def dataset_name(self) -> str:
        return "Tabula Sapiens (skin + blood)"

    def _target_path(self) -> Path:
        return self.output_dir / "tabula_sapiens_benchmark.h5ad"

    def _download(self) -> Path:
        try:
            import cellxgene_census
        except ImportError as exc:
            raise RuntimeError(
                "cellxgene-census is required. "
                "Install with: pip install cellxgene-census"
            ) from exc
        try:
            import scanpy as sc
        except ImportError as exc:
            raise RuntimeError(
                "scanpy is required. Install with: pip install scanpy"
            ) from exc

        logger.info(f"Opening CELLxGENE Census...")
        tissue_filter = " or ".join(
            f'tissue_general=="{t}"' for t in BENCHMARK_TISSUES
        )
        obs_filter = (
            f'({tissue_filter}) '
            f'and collection_id=="{TABULA_SAPIENS_COLLECTION_ID}" '
            f'and is_primary_data==True'
        )

        with cellxgene_census.open_soma() as census:
            logger.info(f"Fetching subset (filter: {obs_filter[:80]}...)")
            adata = cellxgene_census.get_anndata(
                census,
                organism="Homo sapiens",
                obs_value_filter=obs_filter,
                obs_column_names=["cell_type", "tissue_general", "dataset_id"],
            )

        # Subsample to MAX_CELLS_PER_TYPE per annotated cell type
        import pandas as pd
        groups = adata.obs.groupby("cell_type", group_keys=False)
        sampled_idx = groups.apply(
            lambda g: g.sample(n=min(len(g), MAX_CELLS_PER_TYPE), random_state=42)
        ).index
        adata = adata[sampled_idx].copy()

        # Log-normalize
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

        target = self._target_path()
        adata.write_h5ad(target)
        logger.info(
            f"Saved Tabula Sapiens subset to {target} "
            f"({len(adata)} cells, {adata.obs['cell_type'].nunique()} cell types)"
        )
        return target
