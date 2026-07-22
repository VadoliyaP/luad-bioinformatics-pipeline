import streamlit as st
import pandas as pd
import plotly.express as px

# Page configuration
st.set_page_config(page_title="LUAD BioTarget Explorer", layout="wide")

st.title("🧬 Lung Adenocarcinoma (LUAD) Target Identification Pipeline")
st.markdown("### Co-expression Module Detection & Druggability Profiling")

# Load annotated druggability data
@st.cache_data
def load_data():
    return pd.read_csv("target_druggability_annotated.csv")

try:
    df = load_data()
    
    # ---------------------------------------------------------
    # Sidebar Filters
    # ---------------------------------------------------------
    st.sidebar.header("Filter & Search Targets")
    
    # Search box for gene symbol
    search_gene = st.sidebar.text_input("Search Gene Symbol:", "").strip().upper()
    
    # Multiselect for Modules
    selected_module = st.sidebar.multiselect(
        "Module Color:", 
        options=df["ModuleColor"].unique(), 
        default=df["ModuleColor"].unique()
    )
    
    # Multiselect for Druggability Tiers
    selected_tier = st.sidebar.multiselect(
        "Druggability Tier:", 
        options=df["Druggability_Tier"].unique(), 
        default=df["Druggability_Tier"].unique()
    )

    # Apply filters
    filtered_df = df[
        (df["ModuleColor"].isin(selected_module)) & 
        (df["Druggability_Tier"].isin(selected_tier))
    ]
    
    if search_gene:
        filtered_df = filtered_df[filtered_df["GeneSymbol"].str.upper().str.contains(search_gene, na=False)]

    # ---------------------------------------------------------
    # Key Metrics Display
    # ---------------------------------------------------------
    col1, col2, col3 = st.columns(3)
    col1.metric("Total Targets Displayed", len(filtered_df))
    
    tier_1_2_count = len(filtered_df[filtered_df["Druggability_Tier"].str.startswith(("Tier 1", "Tier 2"), na=False)])
    col2.metric("Tier 1 & 2 (Validated / Druggable)", tier_1_2_count)
    
    tier_3_count = len(filtered_df[filtered_df["Druggability_Tier"].str.startswith("Tier 3", na=False)])
    col3.metric("Tier 3 (Novel / Undrugged)", tier_3_count)

    st.markdown("---")

    # ---------------------------------------------------------
    # Visualization Section
    # ---------------------------------------------------------
    st.subheader("📊 Target Distribution by Module & Druggability Tier")
    
    if not filtered_df.empty:
        fig = px.histogram(
            filtered_df, 
            x="ModuleColor", 
            color="Druggability_Tier", 
            barmode="group",
            title="Druggability Breakdown per WGCNA Co-expression Module",
            labels={"ModuleColor": "WGCNA Module", "count": "Gene Count"},
            color_discrete_sequence=px.colors.qualitative.Set2
        )
        # Replaced deprecated use_container_width=True with width="stretch"
        st.plotly_chart(fig, width="stretch")
    else:
        st.info("No targets match the current filter selection.")

    # ---------------------------------------------------------
    # Interactive Data Table & Download
    # ---------------------------------------------------------
    st.subheader("📋 Annotated Target Explorer Table")
    
    # Column selector for display customization
    display_cols = ["GeneSymbol", "ModuleColor", "Druggability_Tier", "Categories", "Sample_Drugs"]
    available_cols = [c for c in display_cols if c in filtered_df.columns]
    
    st.dataframe(filtered_df[available_cols], width="stretch")

    # CSV Export Button
    csv = filtered_df.to_csv(index=False).encode('utf-8')
    st.download_button(
        label="📥 Download Filtered Target List (CSV)",
        data=csv,
        file_name="luad_filtered_targets.csv",
        mime="text/csv"
    )

except FileNotFoundError:
    st.error("⚠️ Could not find 'target_druggability_annotated.csv'. Please make sure Phase 4 script has completed successfully.")