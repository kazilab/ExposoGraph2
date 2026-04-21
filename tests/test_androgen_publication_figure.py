"""Tests for the androgen publication figure renderer."""

from pathlib import Path
from types import SimpleNamespace

import matplotlib

from ExposoGraph import build_androgen_module_graph, render_androgen_publication_figure

matplotlib.use("Agg")


def test_render_androgen_publication_figure_writes_svg(tmp_path: Path):
    graph = build_androgen_module_graph()
    figure, saved_paths = render_androgen_publication_figure(
        graph,
        output_dir=tmp_path,
        showcase_summary=SimpleNamespace(node_count=107, edge_count=132),
        formats=("svg",),
    )
    try:
        output_path = saved_paths["svg"]
        assert output_path.exists()
        assert output_path.stat().st_size > 0

        svg = output_path.read_text(encoding="utf-8")
        assert "Androgen receptor" in svg
        assert "AR proliferative" in svg
        assert "UGT2B17 CN deletion" in svg
    finally:
        import matplotlib.pyplot as plt

        plt.close(figure)
