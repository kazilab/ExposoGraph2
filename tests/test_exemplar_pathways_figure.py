"""Tests for the Fig. 3 exemplar pathways renderer."""

from pathlib import Path

import matplotlib

from ExposoGraph import (
    build_androgen_module_graph,
    build_full_legends_graph,
    render_exemplar_pathways_figure,
)

matplotlib.use("Agg")


def test_render_exemplar_pathways_figure_writes_svg(tmp_path: Path):
    figure, saved_paths = render_exemplar_pathways_figure(
        showcase_graph=build_full_legends_graph(),
        androgen_graph=build_androgen_module_graph(),
        output_dir=tmp_path,
        formats=("svg",),
    )
    try:
        output_path = saved_paths["svg"]
        assert output_path.exists()
        assert output_path.stat().st_size > 0

        svg = output_path.read_text(encoding="utf-8")
        assert "Benzo[a]pyrene" in svg
        assert "4-Aminobiphenyl" in svg
        assert "Testosterone" in svg
        assert "Aflatoxin B1" in svg
        assert "SRD5A2 V89L" in svg
    finally:
        import matplotlib.pyplot as plt

        plt.close(figure)


