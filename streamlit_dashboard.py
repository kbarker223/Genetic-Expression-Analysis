import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px

st.set_page_config(page_title="Gene Therapy Off-Target Safety", layout="wide")

# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

@st.cache_data
def load_pca_summary():
    df = pd.read_csv("data/pca_summary.csv")
    return df

@st.cache_data
def load_expression():
    return pd.read_parquet("data/gene_expression_matrix.parquet")

@st.cache_data
def load_disease_associations():
    return pd.read_parquet("data/disease_associations.parquet")

@st.cache_data
def load_sample_by_gene():
    return pd.read_parquet("data/sample_by_gene.parquet")

@st.cache_data
def compute_risk_table(pca_df, sbg_df, da_df):
    """Build the per-treatment risk summary table."""
    control_row = sbg_df[sbg_df["sample_type"] == "Control"]
    control_mean = control_row.drop(columns="sample_type").mean()

    disease_genes = set(da_df["gene_id"].unique())

    # Map gene ENSG IDs to association scores (max per gene)
    gene_score = da_df.groupby("gene_id")["association_score"].max()

    rows = []
    treatments = pca_df[~pca_df["sample_type"].isin(["Control", "Base"])]["sample_type"].tolist()

    for tx in treatments:
        tx_rows = sbg_df[sbg_df["sample_type"] == tx]
        if tx_rows.empty:
            continue
        tx_mean = tx_rows.drop(columns="sample_type").mean()

        # log2 fold-change relative to control (add pseudocount)
        fc = np.log2((tx_mean + 0.01) / (control_mean + 0.01))

        disturbed = fc[fc.abs() > 1].index.tolist()
        breadth = len(disturbed)

        # Evidence-weighted burden: sum of association scores for disturbed disease genes
        disturbed_disease = [g for g in disturbed if g in disease_genes]
        burden = gene_score.reindex(disturbed_disease).sum()

        high_conf = len([g for g in disturbed_disease if gene_score.get(g, 0) >= 0.5])

        # Therapeutic area coverage
        ta_set = set()
        for g in disturbed_disease:
            gene_tas = da_df[da_df["gene_id"] == g]["therapeutic_areas"].explode()
            ta_set.update(gene_tas.dropna().tolist())
        ta_coverage = len(ta_set)

        pca_row = pca_df[pca_df["sample_type"] == tx].iloc[0]
        rows.append({
            "treatment": tx,
            "breadth": breadth,
            "evidence_burden": round(burden, 3),
            "high_conf_perturbations": high_conf,
            "therapeutic_area_coverage": ta_coverage,
            "mean_dist": pca_row["mean_dist"],
            "log2FC_SNHG14": pca_row["log2FC"],
            "pareto_optimal": pca_row["pareto_optimal"],
        })

    risk_df = pd.DataFrame(rows)
    # Composite rank: normalize each risk dimension and sum
    for col in ["breadth", "evidence_burden", "high_conf_perturbations", "therapeutic_area_coverage", "mean_dist"]:
        mx = risk_df[col].max()
        risk_df[f"{col}_norm"] = risk_df[col] / mx if mx > 0 else 0
    norm_cols = [c for c in risk_df.columns if c.endswith("_norm")]
    risk_df["composite_risk"] = risk_df[norm_cols].mean(axis=1)
    risk_df = risk_df.drop(columns=norm_cols)
    risk_df["composite_rank"] = risk_df["composite_risk"].rank(ascending=True).astype(int)
    risk_df = risk_df.sort_values("composite_risk", ascending=False).reset_index(drop=True)
    return risk_df

@st.cache_data
def get_top_disturbed_genes(treatment, sbg_df, ge_df, da_df, n=10):
    """Return top N disease-linked disturbed genes for a treatment."""
    control_mean = sbg_df[sbg_df["sample_type"] == "Control"].drop(columns="sample_type").mean()
    tx_mean = sbg_df[sbg_df["sample_type"] == treatment].drop(columns="sample_type").mean()
    fc = np.log2((tx_mean + 0.01) / (control_mean + 0.01))

    disturbed = fc[fc.abs() > 1]
    disease_genes = set(da_df["gene_id"].unique())
    disturbed_disease = disturbed[disturbed.index.isin(disease_genes)]

    gene_score = da_df.groupby("gene_id")["association_score"].max()
    result = pd.DataFrame({
        "gene_id": disturbed_disease.index,
        "log2FC": disturbed_disease.values,
        "association_score": gene_score.reindex(disturbed_disease.index).values,
    })
    # merge gene names
    gene_names = ge_df[["gene_id", "gene_name"]].drop_duplicates()
    result = result.merge(gene_names, on="gene_id", how="left")
    result = result.sort_values("association_score", ascending=False).head(n)
    return result.reset_index(drop=True)


# ---------------------------------------------------------------------------
# Load data once
# ---------------------------------------------------------------------------

pca_df = load_pca_summary()
sbg_df = load_sample_by_gene()
ge_df = load_expression()
da_df = load_disease_associations()

with st.spinner("Computing risk table…"):
    risk_df = compute_risk_table(pca_df, sbg_df, da_df)

treatments = sorted(risk_df["treatment"].tolist())

# Initialize session state
if "selected_treatment" not in st.session_state:
    st.session_state.selected_treatment = treatments[0]

# ---------------------------------------------------------------------------
# Sidebar navigation
# ---------------------------------------------------------------------------

page = st.sidebar.radio("View", ["Pareto Plot", "Treatment Drilldown", "Risk Ranking"])

# Global treatment selector in sidebar — syncs with drilldown page
st.sidebar.divider()
st.session_state.selected_treatment = st.sidebar.selectbox(
    "Selected treatment",
    treatments,
    index=treatments.index(st.session_state.selected_treatment),
    key="sidebar_treatment",
)

# ---------------------------------------------------------------------------
# Page: Pareto Plot
# ---------------------------------------------------------------------------

if page == "Pareto Plot":
    st.title("Pareto Plot — Efficacy vs Off-Target Burden")

    plot_df = pca_df[~pca_df["sample_type"].isin(["Control", "Base"])].copy()
    plot_df = plot_df.merge(
        risk_df[["treatment", "composite_risk"]],
        left_on="sample_type", right_on="treatment", how="left",
    )

    fig = px.scatter(
        plot_df,
        x="log2FC",
        y="mean_dist",
        color="pareto_optimal",
        color_discrete_map={True: "#e63946", False: "#457b9d"},
        hover_name="sample_type",
        hover_data={
            "log2FC": ":.3f",
            "mean_dist": ":.1f",
            "std_dist": ":.1f",
            "composite_risk": ":.3f",
            "pareto_optimal": True,
        },
        labels={
            "log2FC": "log₂FC (SNHG14 suppression)",
            "mean_dist": "Mean PCA distance from Control",
        },
        title="",
    )
    fig.update_traces(marker_size=10)
    fig.update_layout(legend_title_text="Pareto optimal")

    # Annotate Pareto-optimal points
    for _, row in plot_df[plot_df["pareto_optimal"]].iterrows():
        fig.add_annotation(
            x=row["log2FC"], y=row["mean_dist"],
            text=row["sample_type"],
            showarrow=True, arrowhead=2, ax=20, ay=-20,
            font=dict(size=10),
        )

    st.plotly_chart(fig, use_container_width=True)

    st.caption(
        "X-axis: more negative = stronger SNHG14 suppression (efficacy). "
        "Y-axis: larger = more off-target transcriptional perturbation. "
        "Red = Pareto-optimal (best efficacy for given off-target burden)."
    )

# ---------------------------------------------------------------------------
# Page: Treatment Drilldown
# ---------------------------------------------------------------------------

elif page == "Treatment Drilldown":
    tx = st.session_state.selected_treatment
    st.title(f"Treatment Drilldown — {tx}")

    pca_row = pca_df[pca_df["sample_type"] == tx].iloc[0]
    risk_row = risk_df[risk_df["treatment"] == tx].iloc[0]

    col1, col2, col3, col4 = st.columns(4)
    col1.metric("SNHG14 log₂FC", f"{pca_row['log2FC']:.3f}")
    col2.metric("PCA distance (mean ± std)", f"{pca_row['mean_dist']:.1f} ± {pca_row['std_dist']:.1f}")
    col3.metric("Disturbed genes (|log₂FC|>1)", f"{risk_row['breadth']:,}")
    col4.metric("Composite risk score", f"{risk_row['composite_risk']:.3f}")

    col5, col6 = st.columns(2)
    col5.metric("Evidence-weighted burden", f"{risk_row['evidence_burden']:.3f}")
    col6.metric("High-confidence perturbations", f"{risk_row['high_conf_perturbations']:,}")

    st.subheader("Top disturbed disease-linked genes")
    with st.spinner("Loading gene data…"):
        top_genes = get_top_disturbed_genes(tx, sbg_df, ge_df, da_df, n=10)

    if top_genes.empty:
        st.info("No disease-linked disturbed genes found for this treatment.")
    else:
        fig = px.bar(
            top_genes,
            x="log2FC",
            y="gene_name",
            orientation="h",
            color="log2FC",
            color_continuous_scale="RdBu_r",
            hover_data={"gene_id": True, "association_score": ":.3f"},
            labels={"log2FC": "log₂FC vs Control", "gene_name": "Gene"},
        )
        fig.update_layout(yaxis={"categoryorder": "total ascending"}, coloraxis_showscale=False)
        st.plotly_chart(fig, use_container_width=True)

        st.dataframe(
            top_genes[["gene_name", "gene_id", "log2FC", "association_score"]].rename(columns={
                "gene_name": "Gene", "gene_id": "Ensembl ID",
                "log2FC": "log₂FC", "association_score": "Max disease association score",
            }),
            hide_index=True,
            use_container_width=True,
        )

# ---------------------------------------------------------------------------
# Page: Risk Ranking
# ---------------------------------------------------------------------------

elif page == "Risk Ranking":
    st.title("Risk Ranking — All Treatments")

    display_df = risk_df.copy()
    display_df["pareto_label"] = display_df["pareto_optimal"].map({True: "Pareto optimal", False: "Non-Pareto"})

    fig = px.bar(
        display_df,
        x="composite_risk",
        y="treatment",
        orientation="h",
        color="pareto_label",
        color_discrete_map={"Pareto optimal": "#e63946", "Non-Pareto": "#457b9d"},
        hover_data={
            "breadth": True,
            "evidence_burden": ":.3f",
            "high_conf_perturbations": True,
            "therapeutic_area_coverage": True,
            "log2FC_SNHG14": ":.3f",
            "composite_rank": True,
            "pareto_label": False,
        },
        labels={"composite_risk": "Composite risk score", "treatment": "Treatment"},
    )
    fig.update_layout(
        yaxis={"categoryorder": "total ascending"},
        legend_title_text="",
    )
    st.plotly_chart(fig, use_container_width=True)

    st.subheader("Full risk table")
    st.dataframe(
        display_df[[
            "composite_rank", "treatment", "composite_risk",
            "breadth", "evidence_burden", "high_conf_perturbations",
            "therapeutic_area_coverage", "log2FC_SNHG14", "mean_dist", "pareto_optimal",
        ]].rename(columns={
            "composite_rank": "Rank",
            "treatment": "Treatment",
            "composite_risk": "Composite risk",
            "breadth": "Disturbed genes",
            "evidence_burden": "Evidence burden",
            "high_conf_perturbations": "High-conf perturbations",
            "therapeutic_area_coverage": "TA coverage",
            "log2FC_SNHG14": "log₂FC SNHG14",
            "mean_dist": "PCA distance",
            "pareto_optimal": "Pareto optimal",
        }),
        hide_index=True,
        use_container_width=True,
    )
