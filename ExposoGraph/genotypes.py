"""Genotype catalog loaded from genotypes.json.

activity_fraction multiplies an enzyme rate. risk_fold multiplies an
exposure scenario score. Callers must not treat one as the other.
"""

from __future__ import annotations

import json
import warnings
from functools import lru_cache
from pathlib import Path
from typing import Any

_PATH = Path(__file__).resolve().parent / "data" / "genotypes.json"


@lru_cache(maxsize=1)
def load_genotypes() -> dict[str, Any]:
    """Return the genotype catalog."""
    return json.loads(_PATH.read_text(encoding="utf-8"))


def activity_fraction(gene: str, diplotype: str, *, allow_generic: bool = True) -> float:
    """Return the catalog activity fraction for one gene and diplotype.

    Gene phenotype aliases win over the generic PM/IM/NM scale. An
    unrecognized label warns and returns 1.0. A numeric diplotype between
    0 and 2 is a caller override. ``allow_generic=False`` skips the
    generic scale, which is how the interaction engine treated a phenotype
    that was not in that gene's table.
    """
    catalog = load_genotypes()
    key = diplotype.lower().strip()
    gene_block = catalog["genes"].get(gene.upper().strip(), {})
    for phenotype in gene_block.get("phenotypes", {}).values():
        aliases = {str(alias).lower() for alias in phenotype.get("aliases", [])}
        if key in aliases:
            return float(phenotype["activity_fraction"])
    for rule in gene_block.get("contains", []):
        if str(rule["token"]).lower() in key:
            return float(rule["activity_fraction"])
    if not allow_generic:
        return 1.0
    scale_key = catalog["label_aliases"].get(key)
    if scale_key is not None:
        return float(catalog["standard_scale"][scale_key])
    try:
        value = float(diplotype)
    except (TypeError, ValueError):
        value = None
    if value is not None and 0.0 <= value <= 2.0:
        return value
    warnings.warn(
        f"Unrecognized diplotype '{diplotype}' for gene '{gene}'; defaulting to NM (1.0)",
        stacklevel=2,
    )
    return 1.0


def risk_fold_rows() -> list[tuple[str, str, dict[str, float]]]:
    """Return exposure matches that carry a per-member risk fold."""
    rows: list[tuple[str, str, dict[str, float]]] = []
    for gene, block in load_genotypes()["genes"].items():
        for phenotype in block.get("phenotypes", {}).values():
            folds = phenotype.get("risk_folds")
            match = phenotype.get("exposure_match")
            if isinstance(folds, dict) and match:
                rows.append((gene, str(match), {str(k): float(v) for k, v in folds.items()}))
    return rows
