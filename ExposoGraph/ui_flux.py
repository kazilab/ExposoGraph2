"""Tab: quantitative flux engine controls."""

from __future__ import annotations

import json

import streamlit as st

from .flux_engine import (
    CarcinogenClass,
    FluxTissueWeightSource,
    PathwayFluxResult,
    compute_pathway_flux,
)

_DEFAULT_GENOTYPES = {
    "CYP1A1": "NM",
    "CYP1B1": "NM",
    "CYP1A2": "NM",
    "CYP2A6": "NM",
    "CYP2A13": "NM",
    "CYP2E1": "NM",
    "CYP3A4": "NM",
    "EPHX1": "NM",
    "GSTM1": "NM",
    "GSTT1": "NM",
    "GSTP1": "NM",
    "NAT1": "NM",
    "NAT2": "NM",
    "ALDH2": "*1/*1",
    "NQO1": "NM",
}


def _enzyme_rows(result: PathwayFluxResult, *, activation: bool) -> list[dict[str, object]]:
    enzymes = result.activation_enzymes if activation else result.detox_enzymes
    return [
        {
            "enzyme": enzyme.enzyme,
            "flux": enzyme.flux,
            "fraction": enzyme.fraction,
            "genotype": enzyme.genotype_modifier,
            "tissue": enzyme.tissue_weight,
            "induction": enzyme.induction_modifier,
            "confidence": enzyme.confidence,
            "model": enzyme.model_kind,
        }
        for enzyme in enzymes
    ]


def render() -> None:
    """Render the quantitative flux engine tab."""
    st.markdown("#### Quantitative Flux Engine")

    col_a, col_b, col_c = st.columns(3)
    with col_a:
        carcinogen_class = st.selectbox(
            "Carcinogen class",
            [item.value for item in CarcinogenClass],
            index=0,
            key="flux_carcinogen_class",
        )
        tissue = st.selectbox(
            "Tissue",
            ["Liver", "Lung", "Bladder", "Breast", "Colon", "Kidney", "Prostate", "Esophagus"],
            index=1,
            key="flux_tissue",
        )
    with col_b:
        tissue_source = st.radio(
            "Tissue weights",
            [FluxTissueWeightSource.CURATED.value, FluxTissueWeightSource.GTEX.value],
            horizontal=True,
            key="flux_tissue_source",
        )
        concentration = st.number_input(
            "Substrate concentration (uM)",
            min_value=0.0,
            value=0.0,
            step=0.01,
            key="flux_substrate_concentration",
        )
        qivive = st.checkbox("Apply QIVIVE scaling", value=False, key="flux_qivive")
    with col_c:
        smoking = st.checkbox("Smoking induction", value=False, key="flux_lifestyle_smoking")
        alcohol_heavy = st.checkbox(
            "Heavy alcohol induction",
            value=False,
            key="flux_lifestyle_alcohol",
        )
        dioxin = st.checkbox("Dioxin/AhR induction", value=False, key="flux_lifestyle_dioxin")

    genotype_text = st.text_area(
        "Genotypes / phenotypes JSON",
        value=json.dumps(_DEFAULT_GENOTYPES, indent=2),
        height=260,
        key="flux_genotype_json",
    )

    try:
        genotypes = json.loads(genotype_text)
        if not isinstance(genotypes, dict):
            raise ValueError("genotypes must be a JSON object")
    except (json.JSONDecodeError, ValueError) as exc:
        st.error(f"Could not parse genotype JSON: {exc}")
        return

    lifestyle = {
        "smoking": smoking,
        "alcohol_heavy": alcohol_heavy,
        "dioxin_exposed": dioxin,
    }
    substrate_conc = concentration if concentration > 0 else None

    result = compute_pathway_flux(
        carcinogen_class,
        {str(gene): str(value) for gene, value in genotypes.items()},
        tissue,
        substrate_conc,
        tissue_weight_source=tissue_source,
        lifestyle=lifestyle,
        qivive=qivive,
    )

    metric_a, metric_b, metric_c, metric_d = st.columns(4)
    metric_a.metric("Activation flux", f"{result.total_activation:.4g}")
    metric_b.metric("Detox flux", f"{result.total_detox:.4g}")
    metric_c.metric("Activation/detox", f"{result.net_ratio:.4g}")
    metric_d.metric("log2 score", f"{result.susceptibility_score_log2:.4g}")

    st.caption(
        f"{result.risk_classification.value} risk class · {result.model_kind} · "
        f"{result.parameter_source} · {result.unit_note}"
    )
    if result.warnings:
        st.warning(", ".join(result.warnings))

    col_left, col_right = st.columns(2)
    with col_left:
        st.markdown("**Activation terms**")
        st.dataframe(
            _enzyme_rows(result, activation=True),
            use_container_width=True,
            hide_index=True,
        )
    with col_right:
        st.markdown("**Detox / repair terms**")
        st.dataframe(
            _enzyme_rows(result, activation=False),
            use_container_width=True,
            hide_index=True,
        )

    st.markdown("**PBPK steady-state concentrations**")
    ss_a, ss_b, ss_c, ss_d = st.columns(4)
    ss_a.metric(
        "Central substrate",
        f"{result.steady_state_concentrations_uM.get('central_substrate_uM', 0):.4g} uM",
    )
    ss_b.metric(
        "Tissue substrate",
        f"{result.steady_state_concentrations_uM.get('tissue_substrate_uM', 0):.4g} uM",
    )
    ss_c.metric(
        "Reactive intermediate",
        f"{result.steady_state_concentrations_uM.get('reactive_intermediate_uM', 0):.4g} uM",
    )
    ss_d.metric(
        "Detoxified metabolite",
        f"{result.steady_state_concentrations_uM.get('detoxified_metabolite_uM', 0):.4g} uM",
    )
    with st.expander("PBPK steady-state model context"):
        st.json(result.steady_state_model)
    if result.induction_factors_used:
        st.markdown("**Induction factors used**")
        st.json(result.induction_factors_used)
    if result.qivive_applied:
        st.markdown("**QIVIVE context**")
        st.json(result.qivive_context)
