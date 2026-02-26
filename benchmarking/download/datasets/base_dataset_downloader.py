"""
Base Dataset Downloader
=======================
Abstract base class for all benchmark dataset downloaders.

Design principle: strictly separate from BioDataDownloader (src/download/),
which is reserved for marker database and reference file acquisition only.
"""
from __future__ import annotations

import logging
from abc import ABC, abstractmethod
from pathlib import Path

logger = logging.getLogger(__name__)


class BaseDatasetDownloader(ABC):
    """
    Abstract base for benchmark dataset downloaders.

    Subclasses implement _download() to fetch one specific dataset.
    The public download() method handles idempotency: if the expected
    output file already exists it is returned immediately without re-downloading.

    Parameters
    ----------
    output_dir : str | Path
        Directory where downloaded files are stored.
    """

    def __init__(self, output_dir: str | Path) -> None:
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def download(self, force: bool = False) -> Path:
        """
        Download dataset if not already present.

        Parameters
        ----------
        force : bool
            If True, re-download even when the target file already exists.

        Returns
        -------
        Path
            Absolute path to the downloaded/cached file.
        """
        target = self._target_path()
        if target.exists() and not force:
            logger.info(f"[{self.__class__.__name__}] Using cached file: {target}")
            return target
        logger.info(f"[{self.__class__.__name__}] Downloading to: {target}")
        return self._download()

    @property
    @abstractmethod
    def dataset_name(self) -> str:
        """Human-readable dataset name for logging."""

    @abstractmethod
    def _target_path(self) -> Path:
        """Expected local path of the primary downloaded artifact."""

    @abstractmethod
    def _download(self) -> Path:
        """Perform the download. Must return local path on success."""
