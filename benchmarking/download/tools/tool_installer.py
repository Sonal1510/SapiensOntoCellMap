"""
Benchmarking Tool Installer / Availability Checker
====================================================
Verifies that required annotation tools are installed and accessible.
Provides clear install instructions when tools are missing.

This module does NOT install tools automatically — it checks availability
and raises informative errors so the user can install the correct version
for their environment.
"""
from __future__ import annotations

import logging
import shutil
import subprocess
import sys
from dataclasses import dataclass

logger = logging.getLogger(__name__)


@dataclass
class ToolSpec:
    """Specification for one benchmarking tool."""
    name: str
    interface: str                    # "python" or "r"
    python_package: str | None        # pip package name (python tools only)
    python_import: str | None         # import name to verify (python tools only)
    r_package: str | None             # R package name (R tools only)
    install_instructions: str         # shown when tool is missing


TOOL_REGISTRY: dict[str, ToolSpec] = {
    "celltypist": ToolSpec(
        name="CellTypist",
        interface="python",
        python_package="celltypist",
        python_import="celltypist",
        r_package=None,
        install_instructions="pip install celltypist",
    ),
    "scsa": ToolSpec(
        name="SCSA",
        interface="python",
        python_package="scsa",
        python_import="scsa",
        r_package=None,
        install_instructions=(
            "pip install scsa\n"
            "  OR: git clone https://github.com/bioinfo-ibms-pumc/SCSA && pip install -e SCSA/"
        ),
    ),
    "sctype": ToolSpec(
        name="scType",
        interface="r",
        python_package=None,
        python_import=None,
        r_package="HGNChelper",      # indirect check — scType has no CRAN package
        install_instructions=(
            "In R: install.packages(c('HGNChelper', 'openxlsx'))\n"
            "  scType is loaded from GitHub at runtime by run_sctype.R"
        ),
    ),
    "singler": ToolSpec(
        name="SingleR",
        interface="r",
        python_package=None,
        python_import=None,
        r_package="SingleR",
        install_instructions=(
            "In R: BiocManager::install(c('SingleR', 'celldex'))"
        ),
    ),
}


class ToolInstaller:
    """
    Checks availability of benchmarking tools and reports missing ones.

    Parameters
    ----------
    tools : list[str]
        Tool keys from TOOL_REGISTRY to check. E.g. ['celltypist', 'sctype'].
    """

    def __init__(self, tools: list[str]) -> None:
        unknown = [t for t in tools if t not in TOOL_REGISTRY]
        if unknown:
            raise ValueError(f"Unknown tools: {unknown}. Known: {list(TOOL_REGISTRY)}")
        self.tools = tools

    def check_all(self, raise_on_missing: bool = False) -> dict[str, bool]:
        """
        Check availability of all specified tools.

        Returns
        -------
        dict[str, bool]
            {tool_key: is_available}
        """
        results = {}
        for key in self.tools:
            spec = TOOL_REGISTRY[key]
            if spec.interface == "python":
                available = self._check_python(spec)
            else:
                available = self._check_r(spec)
            results[key] = available
            if not available:
                logger.warning(
                    f"[ToolInstaller] {spec.name} not found. "
                    f"Install with:\n  {spec.install_instructions}"
                )
                if raise_on_missing:
                    raise RuntimeError(
                        f"{spec.name} is required but not installed. "
                        f"Install with:\n  {spec.install_instructions}"
                    )
        return results

    @staticmethod
    def _check_python(spec: ToolSpec) -> bool:
        try:
            __import__(spec.python_import)
            return True
        except ImportError:
            return False

    @staticmethod
    def _check_r(spec: ToolSpec) -> bool:
        if shutil.which("Rscript") is None:
            logger.warning("[ToolInstaller] Rscript not found in PATH — R tools unavailable.")
            return False
        result = subprocess.run(
            ["Rscript", "-e", f"library({spec.r_package})"],
            capture_output=True,
            text=True,
        )
        return result.returncode == 0
