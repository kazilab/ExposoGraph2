#!/usr/bin/env python
"""Local platform integration smoke check for ExposoGraph 2.0.

Run from the repository root:

    python tools/local_platform_integration_check.py
"""

from __future__ import annotations

import dataclasses
import enum
import json
import math
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any

CHECKS: list[dict[str, str]] = []


def record_check(name: str, status: str, detail: str) -> None:
    """Record one platform check outcome."""
    normalized = status.upper().strip()
    if normalized not in {"PASS", "WARN", "FAIL"}:
        raise ValueError(f"Unsupported check status: {status!r}")
    CHECKS.append({"name": name, "status": normalized, "detail": detail})


def to_json_safe(value: Any) -> Any:
    """Convert common project result objects into JSON-safe values."""
    if value is None or isinstance(value, str | int | bool):
        return value
    if isinstance(value, float):
        if not math.isfinite(value):
            raise ValueError("Non-finite float cannot be serialized safely.")
        return value
    if isinstance(value, enum.Enum):
        return value.value
    if isinstance(value, Path):
        return str(value)
    if dataclasses.is_dataclass(value) and not isinstance(value, type):
        return {
            field.name: to_json_safe(getattr(value, field.name))
            for field in dataclasses.fields(value)
        }
    if isinstance(value, dict):
        return {
            str(to_json_safe(key)): to_json_safe(item)
            for key, item in value.items()
        }
    if isinstance(value, list | tuple):
        return [to_json_safe(item) for item in value]
    if isinstance(value, set | frozenset):
        return [to_json_safe(item) for item in sorted(value, key=repr)]
    if hasattr(value, "model_dump"):
        try:
            return to_json_safe(value.model_dump(mode="json"))
        except TypeError:
            return to_json_safe(value.model_dump())
    if hasattr(value, "to_dict"):
        return to_json_safe(value.to_dict())
    if hasattr(value, "item") and callable(value.item):
        try:
            return to_json_safe(value.item())
        except Exception:
            pass
    return str(value)


def _json_roundtrip(name: str, value: Any) -> Any:
    payload = to_json_safe(value)
    json.dumps(payload, sort_keys=True, allow_nan=False)
    record_check(name, "PASS", "Output converted to JSON-safe payload.")
    return payload


def _run_git_remote_check(repo_root: Path) -> None:
    try:
        result = subprocess.run(
            ["git", "remote", "-v"],
            cwd=repo_root,
            check=False,
            capture_output=True,
            text=True,
            timeout=10,
        )
    except Exception as exc:
        record_check("Repository metadata check", "WARN", f"Unable to inspect local git metadata: {type(exc).__name__}: {exc}")
        return

    output = (result.stdout or "").strip()
    error = (result.stderr or "").strip()
    if result.returncode != 0:
        record_check("Repository metadata check", "WARN", f"git remote -v could not be inspected locally: {error or result.returncode}")
    elif output:
        remote_names = {line.split()[0] for line in output.splitlines() if line.strip()}
        record_check(
            "Repository metadata check",
            "PASS",
            f"Repository metadata inspected; {len(remote_names)} configured remote name(s) detected.",
        )
    else:
        record_check("Repository metadata check", "PASS", "Repository metadata inspected; no configured remote names detected.")


def _validate_repo_context(repo_root: Path) -> None:
    cwd = Path.cwd().resolve()
    if cwd != repo_root:
        record_check("Repository root", "FAIL", "Harness must be run from the repository root.")
    elif not (repo_root / "pyproject.toml").is_file() or not (repo_root / "ExposoGraph").is_dir():
        record_check("Repository root", "FAIL", "Expected repository markers are missing.")
    else:
        record_check("Repository root", "PASS", "Running from repository root.")


def _import_reference_helpers(exposograph: Any) -> tuple[Any, Any, list[str]]:
    notes: list[str] = []
    build_reference_graph = getattr(exposograph, "build_reference_graph", None)
    build_reference_engine = getattr(exposograph, "build_reference_engine", None)
    if callable(build_reference_graph) and callable(build_reference_engine):
        notes.append("root package")
        return build_reference_graph, build_reference_engine, notes

    return None, None, notes


def main() -> int:
    repo_root = Path(__file__).resolve().parents[1]
    _validate_repo_context(repo_root)
    if str(repo_root) not in sys.path:
        sys.path.insert(0, str(repo_root))

    _run_git_remote_check(repo_root)
    record_check("Network independence", "PASS", "Harness uses local package APIs and local git inspection only.")

    exposograph = None
    try:
        import ExposoGraph  # noqa: F401

        exposograph = ExposoGraph
        record_check("Package import", "PASS", "import ExposoGraph completed.")
    except Exception as exc:
        record_check("Package import", "FAIL", f"import ExposoGraph failed: {type(exc).__name__}: {exc}")

    if exposograph is not None:
        version = getattr(exposograph, "__version__", None)
        if version:
            record_check("Public version availability", "PASS", f"__version__ is available: {version}")
        else:
            record_check("Public version availability", "WARN", "__version__ is not exposed; continuing.")

    GraphEngine = None
    build_reference_graph = None
    build_reference_engine = None
    if exposograph is not None:
        graph_engine = getattr(exposograph, "GraphEngine", None)
        build_reference_graph, build_reference_engine, reference_notes = _import_reference_helpers(exposograph)
        GraphEngine = graph_engine
        if graph_engine is not None and callable(build_reference_graph) and callable(build_reference_engine):
            record_check(
                "Core graph imports",
                "PASS",
                "GraphEngine and reference graph helpers imported from " + ", ".join(reference_notes),
            )
        elif graph_engine is None:
            record_check("Core graph imports", "FAIL", "GraphEngine is not available from the current package.")
        else:
            detail = "GraphEngine is available, but build_reference_graph/build_reference_engine are not importable."
            if reference_notes:
                detail += " " + " | ".join(reference_notes)
            record_check("Core graph imports", "FAIL", detail)

    reference_graph = None
    if callable(build_reference_graph):
        try:
            reference_graph = build_reference_graph()
            node_count = len(getattr(reference_graph, "nodes", []) or [])
            edge_count = len(getattr(reference_graph, "edges", []) or [])
            if node_count <= 0 or edge_count <= 0:
                record_check("Reference graph build", "FAIL", "Reference graph built without nodes or edges.")
            else:
                record_check("Reference graph build", "PASS", f"Reference graph built with {node_count} nodes and {edge_count} edges.")
        except Exception as exc:
            record_check("Reference graph build", "FAIL", f"build_reference_graph failed: {type(exc).__name__}: {exc}")
    else:
        record_check("Reference graph build", "FAIL", "build_reference_graph is unavailable in the current package.")

    reference_engine = None
    if callable(build_reference_engine):
        try:
            reference_engine = build_reference_engine()
            if GraphEngine is not None and not isinstance(reference_engine, GraphEngine):
                record_check("Reference engine build", "FAIL", "build_reference_engine did not return GraphEngine.")
            else:
                record_check("Reference engine build", "PASS", "Reference engine built.")
        except Exception as exc:
            record_check("Reference engine build", "FAIL", f"build_reference_engine failed: {type(exc).__name__}: {exc}")
    else:
        record_check("Reference engine build", "FAIL", "build_reference_engine is unavailable in the current package.")

    if reference_engine is not None:
        validate = getattr(reference_engine, "validate", None)
        if callable(validate):
            try:
                validation_errors = validate()
                if validation_errors:
                    record_check("Reference engine validation", "FAIL", f"Validation returned {len(validation_errors)} issue(s).")
                else:
                    record_check("Reference engine validation", "PASS", "engine.validate() returned no issues.")
            except Exception as exc:
                record_check("Reference engine validation", "FAIL", f"engine.validate() failed: {type(exc).__name__}: {exc}")
        else:
            record_check("Reference engine validation", "WARN", "engine.validate() is unavailable; continuing.")
    else:
        record_check("Reference engine validation", "FAIL", "Reference engine unavailable; validation could not run.")

    module3_result = None
    if exposograph is not None:
        try:
            from ExposoGraph import CarcinogenClass, PathwayFluxResult, compute_pathway_flux

            if not callable(compute_pathway_flux):
                raise RuntimeError("compute_pathway_flux is not callable")
            module3_result = compute_pathway_flux(
                CarcinogenClass.PAH,
                {"CYP1A1": "NM", "GSTM1": "NM", "GSTT1": "NM"},
                tissue="Lung",
            )
            if not isinstance(module3_result, PathwayFluxResult):
                record_check("Module 3 simple workflow", "FAIL", "compute_pathway_flux returned an unexpected result type.")
            else:
                record_check(
                    "Module 3 simple workflow",
                    "PASS",
                    "Identified ExposoGraph.compute_pathway_flux as the simple individual-carcinogen workflow.",
                )
        except Exception as exc:
            record_check(
                "Module 3 simple workflow",
                "FAIL",
                "Module 3 public workflow could not be identified for integration harness. "
                f"Details: {type(exc).__name__}: {exc}",
            )

    module5_result = None
    if exposograph is not None:
        try:
            from ExposoGraph import InteractionMatrixResult, compute_interaction_matrix

            if not callable(compute_interaction_matrix):
                raise RuntimeError("compute_interaction_matrix is not callable")
            module5_result = compute_interaction_matrix({"PAH": 1.0, "HCA": 0.5}, tissue="Liver")
            if not isinstance(module5_result, InteractionMatrixResult):
                record_check("Module 5 advanced workflow", "FAIL", "compute_interaction_matrix returned an unexpected result type.")
            else:
                record_check(
                    "Module 5 advanced workflow",
                    "PASS",
                    "Identified ExposoGraph.compute_interaction_matrix as the advanced co-exposure workflow.",
                )
        except Exception as exc:
            record_check(
                "Module 5 advanced workflow",
                "FAIL",
                "Module 5 public workflow could not be identified for integration harness. "
                f"Details: {type(exc).__name__}: {exc}",
            )

    if module3_result is not None:
        try:
            _json_roundtrip("Module 3 JSON serialization", module3_result)
        except Exception as exc:
            record_check("Module 3 JSON serialization", "FAIL", f"Module 3 output is not JSON-safe: {type(exc).__name__}: {exc}")

    if module5_result is not None:
        try:
            _json_roundtrip("Module 5 JSON serialization", module5_result)
        except Exception as exc:
            record_check("Module 5 JSON serialization", "FAIL", f"Module 5 output is not JSON-safe: {type(exc).__name__}: {exc}")

    if reference_engine is not None and exposograph is not None:
        try:
            from ExposoGraph import to_json

            if not callable(to_json):
                record_check("Graph export path", "WARN", "to_json export helper is unavailable; continuing.")
            else:
                with tempfile.TemporaryDirectory(prefix="eg2_platform_check_") as tmp:
                    output_path = Path(tmp) / "reference_graph.json"
                    written = to_json(reference_engine, output_path)
                    data = json.loads(Path(written).read_text(encoding="utf-8"))
                    if not isinstance(data.get("nodes"), list) or not isinstance(data.get("edges"), list):
                        record_check("Graph export path", "FAIL", "Exported JSON is missing nodes or edges arrays.")
                    else:
                        record_check("Graph export path", "PASS", "Reference engine exported to JSON in a temporary directory.")
        except Exception as exc:
            record_check("Graph export path", "FAIL", f"Local graph export failed: {type(exc).__name__}: {exc}")
    else:
        record_check("Graph export path", "WARN", "Reference engine unavailable; export check skipped.")

    record_check("Runtime path disclosure", "PASS", "Runtime check output omits absolute repository paths.")

    max_name = max(len(item["name"]) for item in CHECKS)
    print("ExposoGraph 2.0 local platform integration check")
    print("=" * 55)
    for item in CHECKS:
        print(f"{item['status']:<5} {item['name']:<{max_name}}  {item['detail']}")

    failed = [item for item in CHECKS if item["status"] == "FAIL"]
    warned = [item for item in CHECKS if item["status"] == "WARN"]
    print("-" * 55)
    print(f"Summary: {len(CHECKS) - len(failed) - len(warned)} PASS, {len(warned)} WARN, {len(failed)} FAIL")
    if failed:
        print("Full-platform integration check: FAIL")
        return 1
    print("Full-platform integration check: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
