import pandas as pd
import numpy as np

# 1. Define Known LUAD FDA-Approved Benchmark Targets (Offline Controls)
KNOWN_LUAD_TARGETS = {
    "EGFR": {"Tier": "Tier 1 (Approved Target)", "Sample_Drugs": "Osimertinib, Erlotinib, Gefitinib"},
    "KRAS": {"Tier": "Tier 1 (Approved Target)", "Sample_Drugs": "Sotorasib, Adagrasib"},
    "ALK":  {"Tier": "Tier 1 (Approved Target)", "Sample_Drugs": "Alectinib, Crizotinib, Lorlatinib"},
    "MET":  {"Tier": "Tier 1 (Approved Target)", "Sample_Drugs": "Capmatinib, Tepotinib"},
    "ERBB2": {"Tier": "Tier 1 (Approved Target)", "Sample_Drugs": "Trastuzumab Deruxtecan"},
    "BRAF": {"Tier": "Tier 1 (Approved Target)", "Sample_Drugs": "Dabrafenib + Trametinib"},
    "RET":  {"Tier": "Tier 1 (Approved Target)", "Sample_Drugs": "Selpercatinib, Pralsetinib"},
    "ROS1": {"Tier": "Tier 1 (Approved Target)", "Sample_Drugs": "Entrectinib, Crizotinib"},
    "CDK4": {"Tier": "Tier 1 (Approved Target)", "Sample_Drugs": "Palbociclib, Abemaciclib"},
    "CDK6": {"Tier": "Tier 1 (Approved Target)", "Sample_Drugs": "Palbociclib, Ribociclib"}
}

def assign_enriched_tier(row):
    gene = str(row['GeneSymbol']).upper()
    drug_count = row.get('Drug_Count', 0)
    gene_family = row.get('Gene_Family', 'Unknown')
    kme = row.get('kME', 0)
    
    # Check 1: Known LUAD Benchmark Drivers (Tier 1)
    if gene in KNOWN_LUAD_TARGETS:
        return KNOWN_LUAD_TARGETS[gene]['Tier'], KNOWN_LUAD_TARGETS[gene]['Sample_Drugs']
    
    # Check 2: API Hits with Active Drug Counts (Tier 1)
    if drug_count > 0:
        return "Tier 1 (Approved/Investigational Target)", row.get('Sample_Drugs', 'API Listed')
        
    # Check 3: Druggable Gene Families (Tier 2)
    druggable_families = ['Kinase', 'GPCR', 'Ion Channel', 'Transporter', 'ABC Transporter', 'Solute Carrier']
    if any(fam in str(gene_family) for fam in druggable_families) or gene.startswith(('MAPK', 'SLC', 'ABCC', 'GPR', 'CDK')):
        return "Tier 2 (Druggable Gene Family)", "Novel Family Hit (No Direct Approved Small Molecule)"
        
    # Check 4: Uncharacterized / Novel targets (Tier 3)
    return "Tier 3 (Novel / Unknown Target)", "None"

# Load current annotated dataset
df = pd.read_csv("wgcna_gene_modules.csv")

# Apply refined tiering logic
tier_results = df.apply(assign_enriched_tier, axis=1)
df['Druggability_Tier'] = [t[0] for t in tier_results]
df['Sample_Drugs'] = [t[1] for t in tier_results]

# Save updated results
df.to_csv("target_druggability_annotated.csv", index=False)
print("Updated Druggability Summary:\n", df['Druggability_Tier'].value_counts())