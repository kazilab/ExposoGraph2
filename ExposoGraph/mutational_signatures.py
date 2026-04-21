"""COSMIC SBS mutational-signature endpoints for ExposoGraph 2.0.

Provides a validation endpoint that links predicted per-class risk scores to
COSMIC mutational signatures (Manuscript v6, Discussion):

    "mutational signatures provide a potential validation endpoint...
    linking ExposoGraph 2.0 predictions to tumor mutational signatures
    (e.g., SBS4 for tobacco, SBS22 for aristolochic acid, experimental
    signatures from COSMIC) represents an important future direction".

The module reads a curated JSON catalogue under :mod:`ExposoGraph.data` and
offers three main entry points:

* :func:`get_mutational_signature_catalog` - the raw JSON catalogue.
* :func:`predict_expected_signatures` - given a risk map keyed by carcinogen
  class, return the SBS signatures predicted to be elevated with weights
  proportional to the input risk magnitudes.
* :func:`score_signature_attribution` - compare predicted signature weights
  against observed tumour signature activities (e.g. SigProfiler output) and
  return per-signature residuals and an overall concordance metric.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

_DATA_FILE = Path(__file__).with_name("data") / "mutational_signatures.json"
_CATALOG_CACHE: dict[str, Any] | None = None

# Weight applied to "secondary" signature attributions relative to "primary"
# signature entries. Primary signatures are experimentally attributed to the
# aetiology; secondary signatures are reported associations.
_SECONDARY_WEIGHT: float = 0.35


def _load_catalog() -> dict[str, Any]:
    global _CATALOG_CACHE
    if _CATALOG_CACHE is None:
        with _DATA_FILE.open("r", encoding="utf-8") as handle:
            _CATALOG_CACHE = json.load(handle)
    return _CATALOG_CACHE


@dataclass
class SignaturePrediction:
    """Predicted relative weight of a single COSMIC SBS signature."""

    signature: str
    weight: float
    aetiology: str
    confidence: str
    contributing_classes: list[str]


@dataclass
class SignatureAttributionScore:
    """Comparison between predicted and observed signature activities."""

    signature: str
    predicted_weight: float
    observed_weight: float
    residual: float


@dataclass
class SignatureAttributionSummary:
    """Aggregate concordance summary across a signature set."""

    per_signature: list[SignatureAttributionScore]
    total_abs_residual: float
    mean_abs_residual: float
    spearman_like_agreement: float


def get_mutational_signature_catalog() -> dict[str, Any]:
    """Return a deep-copy-safe view of the signature JSON catalogue."""
    catalog: dict[str, Any] = json.loads(json.dumps(_load_catalog()))
    return catalog


def carcinogen_class_to_signatures(
    carcinogen_class: str,
) -> dict[str, list[str]]:
    """Return the primary / secondary SBS signatures for a carcinogen class.

    Unknown classes return empty lists rather than raising, so callers can
    safely iterate across a risk map that may include classes outside the
    curated attribution table.
    """
    mapping = _load_catalog().get("carcinogen_class_map", {})
    entry = mapping.get(carcinogen_class, {})
    return {
        "primary": list(entry.get("primary", [])),
        "secondary": list(entry.get("secondary", [])),
    }


def predict_expected_signatures(
    risk_map: dict[str, float],
    *,
    secondary_weight: float = _SECONDARY_WEIGHT,
    normalize: bool = True,
) -> list[SignaturePrediction]:
    """Translate a {carcinogen_class: risk_score} map into expected SBS weights.

    The weight for each signature is::

        weight(sig) = Σ_{class where sig is primary}   risk[class]
                    + secondary_weight * Σ_{class where sig is secondary} risk[class]

    When ``normalize`` is True, the weights are scaled to sum to 1.0 (i.e.
    proportions of total predicted signature burden) to match the relative
    activity reporting used by tools such as SigProfilerAssignment.
    """
    catalog = _load_catalog()
    class_map = catalog.get("carcinogen_class_map", {})
    signature_meta = catalog.get("signatures", {})

    accumulated: dict[str, float] = {}
    contributors: dict[str, set[str]] = {}

    for carcinogen_class, score in risk_map.items():
        if not score or score <= 0:
            continue
        entry = class_map.get(carcinogen_class)
        if not entry:
            continue
        for sig in entry.get("primary", []):
            accumulated[sig] = accumulated.get(sig, 0.0) + float(score)
            contributors.setdefault(sig, set()).add(carcinogen_class)
        for sig in entry.get("secondary", []):
            accumulated[sig] = accumulated.get(sig, 0.0) + secondary_weight * float(score)
            contributors.setdefault(sig, set()).add(carcinogen_class)

    if normalize and accumulated:
        total = sum(accumulated.values())
        if total > 0:
            accumulated = {sig: value / total for sig, value in accumulated.items()}

    predictions: list[SignaturePrediction] = []
    for sig, weight in sorted(accumulated.items(), key=lambda item: -item[1]):
        meta = signature_meta.get(sig, {})
        predictions.append(
            SignaturePrediction(
                signature=sig,
                weight=round(weight, 4),
                aetiology=str(meta.get("aetiology", "unknown")),
                confidence=str(meta.get("confidence", "unknown")),
                contributing_classes=sorted(contributors.get(sig, set())),
            )
        )
    return predictions


def score_signature_attribution(
    predicted: dict[str, float] | list[SignaturePrediction],
    observed: dict[str, float],
) -> SignatureAttributionSummary:
    """Compare predicted signature weights against observed tumour activities.

    ``observed`` is expected to be a {signature: relative_activity} map; both
    sides are normalized internally so scale mismatches (e.g. proportions vs
    absolute mutation counts) do not bias the residuals.
    """
    if isinstance(predicted, list):
        predicted = {entry.signature: entry.weight for entry in predicted}

    def _normalize(values: dict[str, float]) -> dict[str, float]:
        clean = {k: float(v) for k, v in values.items() if v and float(v) > 0}
        total = sum(clean.values())
        if total <= 0:
            return {}
        return {k: v / total for k, v in clean.items()}

    pred_norm = _normalize(predicted)
    obs_norm = _normalize(observed)

    signatures = sorted(set(pred_norm) | set(obs_norm))
    per_signature: list[SignatureAttributionScore] = []
    for sig in signatures:
        p = pred_norm.get(sig, 0.0)
        o = obs_norm.get(sig, 0.0)
        per_signature.append(
            SignatureAttributionScore(
                signature=sig,
                predicted_weight=round(p, 4),
                observed_weight=round(o, 4),
                residual=round(p - o, 4),
            )
        )

    total_abs = sum(abs(entry.residual) for entry in per_signature)
    mean_abs = total_abs / len(per_signature) if per_signature else 0.0

    # Rank-order agreement: fraction of signature pairs that share the same
    # ordering between predicted and observed.
    ranked = sorted(per_signature, key=lambda e: -e.predicted_weight)
    correct = 0
    pairs = 0
    for i in range(len(ranked)):
        for j in range(i + 1, len(ranked)):
            pairs += 1
            if ranked[i].observed_weight >= ranked[j].observed_weight:
                correct += 1
    rank_agreement = correct / pairs if pairs > 0 else 1.0

    return SignatureAttributionSummary(
        per_signature=per_signature,
        total_abs_residual=round(total_abs, 4),
        mean_abs_residual=round(mean_abs, 4),
        spearman_like_agreement=round(rank_agreement, 4),
    )


def summary_to_dict(summary: SignatureAttributionSummary) -> dict[str, Any]:
    """Convert :class:`SignatureAttributionSummary` to a JSON-friendly dict."""
    return {
        "per_signature": [asdict(entry) for entry in summary.per_signature],
        "total_abs_residual": summary.total_abs_residual,
        "mean_abs_residual": summary.mean_abs_residual,
        "spearman_like_agreement": summary.spearman_like_agreement,
    }


__all__ = [
    "SignatureAttributionScore",
    "SignatureAttributionSummary",
    "SignaturePrediction",
    "carcinogen_class_to_signatures",
    "get_mutational_signature_catalog",
    "predict_expected_signatures",
    "score_signature_attribution",
    "summary_to_dict",
]
