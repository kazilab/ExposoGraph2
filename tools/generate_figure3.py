"""Regenerate Figure 3 from the local Figure 3 notebook source."""

from __future__ import annotations

from pathlib import Path
import warnings

import matplotlib
import nbformat as nbf

from create_figure3_notebook import NOTEBOOK_PATH, build_notebook


def _write_notebook(path: Path, notebook: nbf.NotebookNode) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        nbf.write(notebook, handle)


def main() -> None:
    matplotlib.use("Agg")
    warnings.filterwarnings(
        "ignore",
        message="FigureCanvasAgg is non-interactive, and thus cannot be shown",
        category=UserWarning,
    )
    notebook = build_notebook()
    _write_notebook(NOTEBOOK_PATH, notebook)
    print(f"Regenerated notebook: {NOTEBOOK_PATH}")

    namespace: dict[str, object] = {"__name__": "__figure3_generation__"}
    for index, cell in enumerate(notebook.cells, start=1):
        if cell.cell_type != "code":
            continue
        print(f"Executing Figure 3 notebook code cell {index}")
        source = "".join(cell.source)
        exec(compile(source, f"{NOTEBOOK_PATH}#cell-{index}", "exec"), namespace)


if __name__ == "__main__":
    main()
