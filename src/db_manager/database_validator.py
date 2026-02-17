#!/usr/bin/python3
"""
Author          : Sonal Rashmi
Date            : 16/02/2026
Description     : Enforces biological schema rules on the master marker database.
                  Runs after all parsers and concatenation but before saving.
                  Quarantines rows that fail validation rather than silently dropping them.
"""
import pandas as pd
from typing import Tuple, Dict


class DatabaseValidator:
    """Enforces biological schema rules on the master marker database."""

    # === COLUMN RULES ===
    RULES = {
        'db_tissue_id': {
            'type': 'ontology_id',
            'allowed_prefixes': ('UBERON:', 'UBERON_'),
            'nullable': True,
            'description': 'Must be a UBERON ontology ID or NaN'
        },
        'db_cell_id': {
            'type': 'ontology_id',
            'allowed_prefixes': ('CL:', 'CL_'),
            'nullable': True,
            'description': 'Must be a Cell Ontology (CL) ID or NaN'
        },
        'tissue_id': {
            'type': 'ontology_id',
            'allowed_prefixes': ('UBERON:', 'UBERON_'),
            'nullable': True,
            'description': 'Must be a UBERON ontology ID or NaN'
        },
        'cell_id': {
            'type': 'ontology_id',
            'allowed_prefixes': ('CL:', 'CL_'),
            'nullable': True,
            'description': 'Must be a Cell Ontology (CL) ID or NaN'
        },
        'gene': {
            'type': 'gene_symbol',
            'nullable': False,
            'description': 'Must be a valid gene symbol (non-null)'
        },
        'cell_name': {
            'type': 'cell_name',
            'nullable': True,
            'description': 'Must not contain obsolete prefix'
        },
        'source_type': {
            'type': 'controlled_vocab',
            'allowed_values': ['Computational', 'Experiment', 'Literature',
                               'Single-Cell Sequencing', 'Review', 'Company'],
            'nullable': False,
            'description': 'Must be one of the allowed source types'
        },
        'database': {
            'type': 'non_empty_string',
            'nullable': False,
            'description': 'Must be a non-empty database name'
        }
    }

    # Mapping for normalizing source_type to controlled vocabulary
    SOURCE_TYPE_MAP = {
        'computational': 'Computational',
        'experiment': 'Experiment',
        'literature': 'Literature',
        'single-cell sequencing': 'Single-Cell Sequencing',
        'review': 'Review',
        'company': 'Company',
    }

    def validate(self, df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame, Dict]:
        """
        Validate the combined DataFrame.

        Returns:
            clean_df: rows that pass all rules
            quarantine_df: rows that fail (with 'violation_reason' column)
            report: summary dict of all checks
        """
        print("\n--- Running Database Validation ---")
        df = df.copy()
        report: Dict[str, Dict] = {}
        violation_reasons: pd.Series = pd.Series([''] * len(df), index=df.index)

        # --- 1. Normalize source_type (fix before checking) ---
        fixes = self._normalize_source_type(df['source_type'])
        report['source_type_normalized'] = {
            'fixed': int(fixes),
            'description': 'source_type values normalized to controlled vocabulary'
        }

        # --- 2. Check ontology prefix columns ---
        for col in ['db_tissue_id', 'db_cell_id', 'tissue_id', 'cell_id']:
            if col not in df.columns:
                continue
            rule = self.RULES[col]
            bad_mask = self._check_ontology_prefix(df[col], rule['allowed_prefixes'])
            count = int(bad_mask.sum())
            if count > 0:
                # Nullify invalid values (don't quarantine the row)
                df.loc[bad_mask, col] = None
            report[f'{col}_prefix'] = {
                'violations': count,
                'action': 'nullified',
                'description': rule['description']
            }

        # --- 3. Check gene symbol ---
        gene_null_mask, gene_space_mask, gene_non_ascii_mask = self._check_gene_symbol(df['gene'])

        # Quarantine rows with null gene
        null_gene_count = int(gene_null_mask.sum())
        if null_gene_count > 0:
            violation_reasons.loc[gene_null_mask] += 'null_gene;'
        report['gene_null'] = {
            'violations': null_gene_count,
            'action': 'quarantined',
            'description': 'gene is NaN'
        }

        # Flag rows with space-containing genes (likely full names, keep but flag)
        space_gene_count = int(gene_space_mask.sum())
        if space_gene_count > 0:
            if 'gene_quality' not in df.columns:
                df['gene_quality'] = None
            df.loc[gene_space_mask, 'gene_quality'] = 'possible_full_name'
        report['gene_has_spaces'] = {
            'violations': space_gene_count,
            'action': 'flagged (gene_quality column)',
            'description': 'gene contains spaces (may be full name, not symbol)'
        }

        # Flag non-ASCII genes
        non_ascii_count = int(gene_non_ascii_mask.sum())
        if non_ascii_count > 0:
            if 'gene_quality' not in df.columns:
                df['gene_quality'] = None
            # Don't overwrite existing flags
            ascii_only = gene_non_ascii_mask & df['gene_quality'].isna()
            df.loc[ascii_only, 'gene_quality'] = 'non_ascii'
            both = gene_non_ascii_mask & gene_space_mask
            df.loc[both, 'gene_quality'] = 'possible_full_name;non_ascii'
        report['gene_non_ascii'] = {
            'violations': non_ascii_count,
            'action': 'flagged (gene_quality column)',
            'description': 'gene contains non-ASCII characters'
        }

        # --- 4. Check obsolete cell_name ---
        obsolete_mask = self._check_obsolete_cell_name(df['cell_name'])
        obsolete_count = int(obsolete_mask.sum())
        if obsolete_count > 0:
            violation_reasons.loc[obsolete_mask] += 'obsolete_cell_name;'
        report['cell_name_obsolete'] = {
            'violations': obsolete_count,
            'action': 'flagged',
            'description': 'cell_name starts with "obsolete" (should have been resolved by parser)'
        }

        # --- 5. Check controlled vocabulary for source_type ---
        if 'source_type' in df.columns:
            allowed = set(self.RULES['source_type']['allowed_values'])
            bad_source = df['source_type'].notna() & ~df['source_type'].isin(allowed)
            bad_source_count = int(bad_source.sum())
            if bad_source_count > 0:
                violation_reasons.loc[bad_source] += 'invalid_source_type;'
            report['source_type_invalid'] = {
                'violations': bad_source_count,
                'action': 'quarantined',
                'description': 'source_type not in controlled vocabulary after normalization'
            }

            # Check null source_type
            null_source = df['source_type'].isna()
            null_source_count = int(null_source.sum())
            if null_source_count > 0:
                violation_reasons.loc[null_source] += 'null_source_type;'
            report['source_type_null'] = {
                'violations': null_source_count,
                'action': 'quarantined',
                'description': 'source_type is NaN'
            }

        # --- 6. Check non-empty database ---
        if 'database' in df.columns:
            empty_db = df['database'].isna() | (df['database'].astype(str).str.strip() == '')
            empty_db_count = int(empty_db.sum())
            if empty_db_count > 0:
                violation_reasons.loc[empty_db] += 'empty_database;'
            report['database_empty'] = {
                'violations': empty_db_count,
                'action': 'quarantined',
                'description': 'database column is empty or NaN'
            }

        # --- Split into clean and quarantine ---
        quarantine_mask = violation_reasons.str.len() > 0
        quarantine_df = df.loc[quarantine_mask].copy()
        quarantine_df['violation_reason'] = violation_reasons.loc[quarantine_mask]
        clean_df = df.loc[~quarantine_mask].copy()

        # Drop gene_quality column from quarantine (not needed there)
        # but keep it in clean_df for downstream analysis

        report['summary'] = {
            'total_rows': len(df),
            'clean_rows': len(clean_df),
            'quarantined_rows': len(quarantine_df)
        }

        return clean_df, quarantine_df, report

    def _check_ontology_prefix(self, series: pd.Series, allowed_prefixes: tuple) -> pd.Series:
        """
        Check non-null values start with allowed prefix.
        Returns a boolean mask of violating rows.
        """
        non_null = series.notna()
        if not non_null.any():
            return pd.Series(False, index=series.index)
        starts_valid = series.str.startswith(allowed_prefixes, na=False)
        return non_null & ~starts_valid

    def _check_gene_symbol(self, series: pd.Series) -> Tuple[pd.Series, pd.Series, pd.Series]:
        """
        Check gene column for issues.
        Returns:
            null_mask: rows where gene is NaN
            space_mask: rows where gene contains spaces (possible full names)
            non_ascii_mask: rows where gene contains non-ASCII characters
        """
        null_mask = series.isna()
        space_mask = series.astype(str).str.contains(r'\s', na=False) & ~null_mask
        non_ascii_mask = series.astype(str).str.contains(r'[^\x00-\x7F]', na=False) & ~null_mask
        return null_mask, space_mask, non_ascii_mask

    def _check_obsolete_cell_name(self, series: pd.Series) -> pd.Series:
        """
        Detect cell_name starting with "obsolete".
        These should have been resolved by base_parser Step 6.1.
        """
        return series.str.contains(r'(?i)^obsolete\b', na=False)

    def _normalize_source_type(self, series: pd.Series) -> int:
        """
        Normalize source_type to controlled vocabulary using case-insensitive mapping.
        Modifies the series in place (since it's a column reference).
        Returns count of values that were fixed.
        """
        non_null = series.notna()
        if not non_null.any():
            return 0

        lowered = series.str.lower().str.strip()
        mapped = lowered.map(self.SOURCE_TYPE_MAP)
        needs_fix = non_null & mapped.notna() & (series != mapped)
        fixed_count = int(needs_fix.sum())
        if fixed_count > 0:
            series.loc[needs_fix] = mapped.loc[needs_fix]
        return fixed_count

    def _print_report(self, report: Dict) -> None:
        """Print summary table of all validation checks."""
        print("\n" + "=" * 70)
        print("DATABASE VALIDATION REPORT")
        print("=" * 70)

        summary = report.pop('summary', {})

        print(f"\n{'Check':<35} {'Violations':>12} {'Action':<25}")
        print("-" * 70)

        for check_name, details in report.items():
            violations = details.get('violations', details.get('fixed', 0))
            action = details.get('action', 'n/a')
            print(f"{check_name:<35} {violations:>12} {action:<25}")

        print("-" * 70)
        print(f"\nTotal rows:        {summary.get('total_rows', 'N/A')}")
        print(f"Clean rows:        {summary.get('clean_rows', 'N/A')}")
        print(f"Quarantined rows:  {summary.get('quarantined_rows', 'N/A')}")
        print("=" * 70)
