"""Smoke tests for tissue-specific subgraph extraction."""

from ExposoGraph import (
    DEFAULT_THRESHOLD_SWEEP,
    TissueType,
    get_available_gtex_genes,
    get_available_gtex_tissues,
    get_tissue_expression,
    get_tissue_weights,
    get_top_carcinogen_classes_for_tissue,
)


def test_liver_top_carcinogen_classes_list_is_populated():
    top = get_top_carcinogen_classes_for_tissue("Liver")

    assert len(top) >= 3
    assert all(isinstance(item, tuple) and len(item) == 2 for item in top)
    # Liver handles PAHs prominently via CYP1A1/3A4.
    assert any("PAH" in name for name, _ in top)


def test_tissue_expression_lookup_returns_mapping_for_liver():
    expr = get_tissue_expression(TissueType.LIVER)

    assert isinstance(expr, dict)
    assert all(isinstance(v, (int, float)) for v in expr.values())


def test_hpa_tissue_expression_uses_expected_ntpm_sentinel_values():
    assert get_tissue_expression(TissueType.LIVER)["CYP2C9"] == 675.4
    assert get_tissue_expression(TissueType.KIDNEY)["CCBL1"] == 7.6
    assert get_tissue_expression(TissueType.LUNG)["HLA_DPB1"] == 261.7
    # Colon is the mean of HPA GTEx detail rows: Colon - Sigmoid and Colon - Transverse.
    assert get_tissue_expression(TissueType.COLON)["CYP1A1"] == 1.9


def test_tissue_weights_are_max_normalized_from_ntpm_values():
    tissues = get_available_gtex_tissues()
    genes = get_available_gtex_genes()

    expression_by_tissue = {tissue: get_tissue_expression(tissue) for tissue in tissues}
    weights_by_tissue = {tissue: get_tissue_weights(tissue) for tissue in tissues}

    for gene in genes:
        gene_values = [expression_by_tissue[tissue][gene] for tissue in tissues]
        max_value = max(gene_values)
        for tissue in tissues:
            value = expression_by_tissue[tissue][gene]
            expected = (
                0.0
                if max_value == 0 or value / max_value < 0.01
                else round(value / max_value, 4)
            )
            assert weights_by_tissue[tissue][gene] == expected


def test_default_threshold_sweep_is_set():
    assert DEFAULT_THRESHOLD_SWEEP == (0.10, 0.25, 0.50)
