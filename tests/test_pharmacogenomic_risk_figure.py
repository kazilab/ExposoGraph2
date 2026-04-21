from ExposoGraph import build_full_legends_engine
from ExposoGraph.pharmacogenomic_risk_figure import (
    build_pharmacogenomic_risk_class_profiles,
    build_pharmacogenomic_risk_gene_profiles,
    render_pharmacogenomic_risk_figure,
)


def test_risk_gene_profiles_cover_scored_reference_genes():
    engine = build_full_legends_engine()
    profiles = build_pharmacogenomic_risk_gene_profiles(engine)
    profile_by_id = {profile.gene_id: profile for profile in profiles}

    assert len(profiles) >= 10
    assert "CYP2E1" in profile_by_id
    assert profile_by_id["CYP2E1"].activity_score == 1.3
    assert profile_by_id["CYP2E1"].class_count >= 2
    assert "Solvent" in profile_by_id["CYP2E1"].carcinogen_groups


def test_risk_class_profiles_capture_activation_and_detox_support():
    engine = build_full_legends_engine()
    profiles = build_pharmacogenomic_risk_class_profiles(engine)
    profile_by_group = {
        profile.carcinogen_group: profile for profile in profiles
    }

    assert profile_by_group["PAH"].activation_score > 0
    assert profile_by_group["PAH"].repair_score > 0
    assert profile_by_group["Solvent"].detoxification_score > 0
    assert profile_by_group["Alkylating"].scored_gene_count == 0


def test_risk_figure_export_writes_png(tmp_path):
    engine = build_full_legends_engine()
    figure, outputs, gene_profiles, class_profiles = render_pharmacogenomic_risk_figure(
        engine,
        output_dir=tmp_path,
        formats=("png",),
        dpi=120,
    )

    assert figure is not None
    assert outputs["png"].exists()
    assert gene_profiles
    assert class_profiles
