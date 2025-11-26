#! /usr/bin/env python

# make confidence interval around the mean of the data
import argparse
import warnings
import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.colors as pc
from functools import reduce
from multiprocessing import Pool, cpu_count
warnings.simplefilter(action='ignore')

# Define consistent color mapping for transcript types
TRANSCRIPT_TYPE_COLORS = {
    "protein_coding": "#1f77b4",          # blue
    "retained_intron": "#d62728",         # red (updated)
    "nonsense_mediated_decay": "#2ca02c", # green
    "non_stop_decay": "#ff7f0e",          # orange (swapped)
    "TEC": "#9467bd",                     # purple
    "processed_transcript": "#8c564b",    # brown
    "protein_coding_CDS_not_defined": "#e377c2",  # pink
    "protein_coding_LoF": "#7f7f7f",      # gray
    "transcribed_processed_pseudogene": "#bcbd22", # olive
    "transcribed_unprocessed_pseudogene": "#17becf", # teal
}

# Get color with fallback to default color palette
def get_transcript_color(transcript_type, default_colors=pc.qualitative.Plotly):
    """Get color for transcript type with fallback to default color scheme."""
    if transcript_type in TRANSCRIPT_TYPE_COLORS:
        return TRANSCRIPT_TYPE_COLORS[transcript_type]
    
    # Use hash of transcript type for consistent color assignment
    import hashlib
    hash_val = int(hashlib.md5(str(transcript_type).encode()).hexdigest(), 16)
    return default_colors[hash_val % len(default_colors)]

# Parse command-line arguments
argparser = argparse.ArgumentParser()

argparser.add_argument("--goi", help="Gene of interest")
argparser.add_argument("--gtf", help="GTF file")
argparser.add_argument("--kallisto", help="Kallisto file")
argparser.add_argument("--outdir", help="Output directory")

args = argparser.parse_args()

highlight_goi = args.goi
gtf_file = args.gtf
kallisto_file = args.kallisto
outdir = args.outdir


# Step 1: Parse the GTF file
def parse_gtf_attributes(attr_str):
    attrs = {}
    for attr in attr_str.strip().split(';'):
        if attr.strip():
            try:
                parts = attr.strip().split('=', 1)
                if len(parts) == 1:
                    parts = attr.strip().split(None, 1)
                if len(parts) == 2:
                    key, value = parts
                    attrs[key] = value.strip('"\'')
            except ValueError:
                continue
    return attrs

def process_gtf_line(line):
    if line.startswith('#'):
        return None
    fields = line.strip().split('\t')
    if len(fields) < 9 or fields[2] != 'transcript':
        return None

    attrs = parse_gtf_attributes(fields[8])
    transcript_info = {
        'transcript_id': attrs.get('transcript_id', ''),
        'transcript_name': attrs.get('transcript_name', 'unknown'),
        'transcript_type': attrs.get('transcript_type', ''),
        'gene_name': attrs.get('gene_name', 'unknown'),
        'gene_id': attrs.get('gene_id', ''),
        'gene_type': attrs.get('gene_type', ''),
        'gene_biotype': attrs.get('gene_biotype', ''),
        'transcript_source': attrs.get('transcript_source', ''),
        'transcript_id_version': attrs.get('transcript_id', ''),
    }

    if transcript_info['gene_biotype']:
        transcript_info['gene_type'] = transcript_info['gene_biotype']

    if 'transcript_version' in attrs and '.' not in transcript_info['transcript_id_version']:
        transcript_info['transcript_id_version'] += '.' + attrs.get('transcript_version', '1')

    if not transcript_info['transcript_type']:
        transcript_info['transcript_type'] = attrs.get('transcript_biotype', 'unknown')

    if transcript_info['transcript_source'] == 'serafina_manual':
        transcript_info['transcript_id_version'] = attrs.get('transcript_id', '')

    return transcript_info

def load_transcripts_from_gtf(gtf_file, num_workers=None):
    num_workers = num_workers or max(1, cpu_count() - 1)
    with open(gtf_file) as f:
        lines = f.readlines()

    with Pool(num_workers) as pool:
        results = pool.map(process_gtf_line, lines)

    transcripts = [r for r in results if r is not None]
    transcripts_df = pd.DataFrame(transcripts)
    transcripts_df = transcripts_df[transcripts_df['gene_type'] == 'protein_coding']
    print("Transcripts loaded from GTF.")
    return transcripts_df

# Step 2: Parse Kallisto data
def parse_kallisto_line(line):
    acc, tissue, db, work, f, n_reads, n_aligned, p_aligned = line.strip().split("\t")
    df = pd.read_csv(f, sep="\t", header=0, na_filter=False,
                     names=f"transcript_id,length,eff_length,est_counts,{acc}".split(","))[["transcript_id", f"{acc}"]]
    df[f"{acc}"] += 1
    return df

def parse_kallisto(kallisto_file, num_workers=None):
    num_workers = num_workers or max(1, cpu_count() - 1)
    with open(kallisto_file) as f:
        lines = f.readlines()

    with Pool(num_workers) as pool:
        dfs = pool.map(parse_kallisto_line, lines)

    tpm = reduce(lambda left, right: pd.merge(left, right, on='transcript_id', how='outer'), dfs)
    tpm.fillna(0, inplace=True)
    tpm["agg_exp"] = tpm.iloc[:, 1:].median(axis=1)
    
    print(f"Parsed kallisto file: {kallisto_file}")
    return tpm

def merge_expression(kal_data, gtf_data):
    """Merge expression data with GTF metadata."""
    kal_data["expressed"] = kal_data.iloc[:, 1:-1].apply(lambda row: (row >= 1).sum() >= 0, axis=1)

    # Merge Kallisto data with GTF annotations using a right merge to keep only protein_coding genes
    transcript_exp = kal_data.merge(gtf_data, on="transcript_id", how="right")

    # Compute total isoforms per gene
    transcript_exp["total_isoforms"] = transcript_exp.groupby("gene_id")["transcript_id"].transform("count")

    # Compute expression proportions
    transcript_exp["biotype_expression"] = transcript_exp.groupby(["gene_id", "transcript_type"])["agg_exp"].transform("sum")
    transcript_exp["total_gene_expression"] = transcript_exp.groupby("gene_id")["agg_exp"].transform("sum")
    transcript_exp["biotype_expression_proportion"] = transcript_exp["biotype_expression"] / transcript_exp["total_gene_expression"]

    transcript_exp = transcript_exp[['gene_id', 'gene_name', 'transcript_id', 'transcript_type', 
                           'total_isoforms', 'agg_exp', 'biotype_expression',
                           'biotype_expression_proportion', 'total_gene_expression', 'expressed']]
    
    transcript_exp.drop_duplicates(inplace=True)
    
    print("Expression data merged with GTF data.")

    return transcript_exp

def add_transcript_type_stats(expression_df):
    """Compute gtf_transcript_type_count and gtf_transcript_type_proportion for expression_df while keeping all rows."""
    
    # Step 1: Group and count transcript types per gene
    transcript_counts = expression_df.groupby(["gene_id", "gene_name", "transcript_type"]).size().reset_index(name="gtf_transcript_type_count")

    # Step 2: Pivot the table to get transcript counts per type
    transcript_counts_pivot = transcript_counts.pivot_table(
        index=["gene_id", "gene_name"], columns="transcript_type", values="gtf_transcript_type_count", fill_value=0
    ).reset_index()

    # Step 4: Melt to long format for easier merging
    transcript_counts_long = transcript_counts_pivot.melt(
        id_vars=["gene_id", "gene_name"],
        var_name="transcript_type",
        value_name="gtf_transcript_type_count"
    )

    expression_df = expression_df.merge(transcript_counts_long[["gene_id", "transcript_type", "gtf_transcript_type_count"]], on=["gene_id", "transcript_type"], how="left")

    # Step 5: Calculate transcript type proportions
    expression_df["gtf_transcript_type_proportion"] = (
        expression_df["gtf_transcript_type_count"] / expression_df["total_isoforms"]
    )
    
    print("Transcript types annotated.")

    return expression_df

def create_combined_plot(transcript_df, bin_labels, outdir, highlight_goi=None):

    # Extract unique transcript types
    unique_transcript_types = transcript_df["transcript_type"].unique()
    
    # Create consistent color mapping using our helper function
    transcript_type_colors = {t_type: get_transcript_color(t_type) for t_type in unique_transcript_types}

    # Create subplot layout
    fig_combined = make_subplots(
        rows=2, cols=2,
        subplot_titles=(
            "Average Proportion of Isoforms per Transcript Type",
            "Histogram of Total Gene Expression",
            "Average Proportion of Expression per Transcript Type",
            "Proportion of Expressed Transcripts per Transcript Type"
        ),
        shared_xaxes=False
    )

    # Plot 1: Average Proportion of Isoforms per Transcript Type
    for t_type in unique_transcript_types:
        print(t_type)

        # PLOT 1: Average Proportion of Isoforms per Transcript Type According to GTF Annotations
        subset = transcript_df[transcript_df["transcript_type"] == t_type][["expression_bin", "gene_id", "num_genes_in_expression_bin", "gtf_transcript_type_proportion"]].drop_duplicates()

        def normalize_gtf_proportion(group):
            if group.empty:
                return 0

            group = group.drop_duplicates(subset=["gene_id"])
            avg_proportion = group["gtf_transcript_type_proportion"].sum()  # Compute mean proportion
            num_genes = group["num_genes_in_expression_bin"].iloc[0] # Take the first num_genes value
            return avg_proportion / num_genes  # Compute normalized value
        
        normalized = subset.groupby("expression_bin").apply(normalize_gtf_proportion).reset_index(name="normalized_gtf_proportion")

        #print(normalized)

        fig_combined.add_trace(
            go.Scatter(
                x=normalized.index,
                y=normalized["normalized_gtf_proportion"],
                mode="lines+markers",
                name=t_type,
                line=dict(color=transcript_type_colors.get(t_type, "gray")),
                marker=dict(color=transcript_type_colors.get(t_type, "gray")),
                legendgroup=t_type,
            ),
            row=1, col=1
        )

        print(normalized)

        # PLOT 3: Average Proportion of Expression per Transcript Type

        subset = transcript_df[transcript_df["transcript_type"] == t_type][["expression_bin", "gene_id", "num_genes_in_expression_bin", "biotype_expression_proportion", "gene_name"]].drop_duplicates()

        def normalize_expression_proportion(group):
            if group.empty:
                return 0

            group = group.drop_duplicates(subset=["gene_id"])
            avg_proportion = group["biotype_expression_proportion"].sum()
            num_genes = group["num_genes_in_expression_bin"].iloc[0]

            return avg_proportion / num_genes
        
        normalized = subset.groupby("expression_bin").apply(normalize_expression_proportion).reset_index(name="normalized_expression_proportion")


        fig_combined.add_trace(
            go.Scatter(
                x=normalized.index,
                y=normalized["normalized_expression_proportion"],
                mode="lines+markers",
                name=f"{t_type} (Expression)",
                line=dict(color=transcript_type_colors.get(t_type, "gray")),
                marker=dict(color=transcript_type_colors.get(t_type, "gray"), opacity=0.75),
                legendgroup=t_type,
                showlegend=False,
            ),
            row=2, col=1
        )

        # Add this code right after:
        if t_type == "protein_coding" or t_type == "retained_intron":
            # For the confidence band, we need to calculate percentiles per bin
            # Add safety check for empty groups
            lower = subset.groupby("expression_bin").apply(
                lambda x: np.percentile(x["biotype_expression_proportion"], 5) if len(x) > 0 else np.nan
            )
            upper = subset.groupby("expression_bin").apply(
                lambda x: np.percentile(x["biotype_expression_proportion"], 95) if len(x) > 0 else np.nan
            )

            print(lower)
            print(upper)
            
            # Filter out NaN values before plotting
            valid_indices = ~(np.isnan(lower) | np.isnan(upper))
            if valid_indices.any():  # Only proceed if there are valid points
                plot_indices = normalized.index[valid_indices]
                plot_lower = lower[valid_indices]
                plot_upper = upper[valid_indices]
                
                # Get the biotype color and add transparency
                biotype_color = transcript_type_colors.get(t_type, "gray")
                # Convert to rgba with 0.2 opacity
                rgba_color = f"rgba({int(biotype_color[1:3], 16)},{int(biotype_color[3:5], 16)},{int(biotype_color[5:7], 16)},0.2)"
                
                fig_combined.add_trace(
                    go.Scatter(
                        x=list(plot_indices) + list(plot_indices)[::-1],
                        y=list(plot_upper) + list(plot_lower)[::-1],
                        fill='toself',
                        fillcolor=rgba_color,  # Use biotype color with opacity
                        line=dict(color='rgba(255,255,255,0)'),
                        name=f'90% Range ({t_type})',
                        showlegend=True,
                        legendgroup=t_type,
                    ),
                    row=2, col=1
                )

        if highlight_goi:
            highlight_subset = subset[subset["gene_name"] == highlight_goi]
            
            # Group by expression_bin and calculate the proportion directly
            if highlight_subset.shape[0] > 0:
                highlight = highlight_subset.groupby("expression_bin")["biotype_expression_proportion"].first().reset_index()
                print(highlight)

                fig_combined.add_trace(
                go.Scatter(
                    x=highlight.index,
                    y=highlight["biotype_expression_proportion"],
                    mode="markers",
                    marker_symbol="star",
                    name=f"Highlight: {highlight_goi} ({t_type})",
                    marker=dict(size=10, color=transcript_type_colors.get(t_type, "gray"), opacity=1),
                    legendgroup="highlight",
                ),
                row=2, col=1
            )

        # PLOT 4: Proportion of Expressed Transcripts Stratified by Total Gene Expression

        subset = transcript_df[transcript_df["transcript_type"] == t_type][["expression_bin", "gene_id", "num_genes_in_expression_bin", "biotype_expression_proportion"]].drop_duplicates()

        def normalize_expression_proportion(group):
            if group.empty:
                return 0

            group = group.drop_duplicates(subset=["gene_id"])
            avg_proportion = group["biotype_expression_proportion"].sum()
            num_genes = group["num_genes_in_expression_bin"].iloc[0]

            return avg_proportion / num_genes
        
        normalized = subset.groupby("expression_bin").apply(normalize_expression_proportion).reset_index(name="normalized_expression_proportion")

        fig_combined.add_trace(
            go.Scatter(
                x=normalized.index,
                y=normalized["normalized_expression_proportion"],
                mode="lines+markers",
                name=f"{t_type} (Expression)",
                line=dict(color=transcript_type_colors.get(t_type, "gray")),
                marker=dict(color=transcript_type_colors.get(t_type, "gray"), opacity=0.75),
                legendgroup=t_type,
                showlegend=False,
            ),
            row=2, col=2
        )

        print(normalized)

    # Plot 2: Histogram of Total Gene Expression (instead of isoforms)
    fig_combined.add_trace(
        go.Histogram(
            x=transcript_df.groupby("gene_id")["total_gene_expression"].first(),
            name="Total Gene Expression",
            opacity=0.7,
            marker=dict(color="blue"),
        ),
        row=1, col=2
    )

    # Update layout
    # Set x-axis category order explicitly
    _, expression_bin_labels = bin_labels

    # Ensure the correct bin labels are used for the x-axes
    fig_combined.update_xaxes(
        categoryorder="array",
        categoryarray=[str(lbl) for lbl in expression_bin_labels],
        tickvals=list(range(len(expression_bin_labels))),
        ticktext=expression_bin_labels,
        title_text="Binned Gene Expression",
        row=1, col=1,
    )
    fig_combined.update_xaxes(
        categoryorder="array",
        categoryarray=[str(lbl) for lbl in expression_bin_labels],
        tickvals=list(range(len(expression_bin_labels))),
        ticktext=expression_bin_labels,
        title_text="Binned Gene Expression",
        row=2, col=1,
    )
    fig_combined.update_xaxes(
        categoryorder="array",
        categoryarray=[str(lbl) for lbl in expression_bin_labels],  # Fixed closing bracket here
        tickvals=list(range(len(expression_bin_labels))),
        ticktext=expression_bin_labels,
        title_text="Binned Gene Expression",
        row=2, col=2,
    )

    fig_combined.update_yaxes(type="log", row=1, col=2)

    fig_combined.update_layout(
        title="Combined Visualization: Isoforms, Expression Histogram, and Expressed Transcripts",
        xaxis_title="Binned Gene Expression",
        xaxis2_title="Total Gene Expression (Histogram)",
        xaxis3_title="Binned Gene Expression",
        xaxis4_title="Binned Gene Expression",
        yaxis_title="Proportion of Isoforms",
        yaxis2_title="Log Count of Genes",
        yaxis3_title="Average Proportion of Expression by Transcript Type",
        yaxis4_title="Average Proportion of Annotated Transcripts<br>Expressed by Transcript Type",
        barmode="overlay",
        template="plotly_white",
        width=1500,
        height=1000
    )

    # Save the plot
    fig_combined.write_html(f"{outdir}/{highlight_goi}_combined_visualization.html")
    fig_combined.write_image(f"{outdir}/{highlight_goi}_combined_visualization.png", scale=2)
    fig_combined.write_image(f"{outdir}/{highlight_goi}_combined_visualization.svg", scale=2)  # Save as SVG

def calculate_protein_coding_to_retained_intron_ratio(expression_df):
    """Calculate the ratio of protein_coding to retained_intron biotype expression proportion for each gene."""
    
    # Get protein_coding and retained_intron proportions per gene
    gene_ratios = []
    
    for gene_id in expression_df['gene_id'].unique():
        gene_data = expression_df[expression_df['gene_id'] == gene_id]
        gene_name = gene_data['gene_name'].iloc[0]
        total_gene_expression = gene_data['total_gene_expression'].iloc[0]
        
        # Get protein_coding proportion
        pc_data = gene_data[gene_data['transcript_type'] == 'protein_coding']
        pc_proportion = pc_data['biotype_expression_proportion'].iloc[0] if len(pc_data) > 0 else 0
        
        # Get retained_intron proportion
        ri_data = gene_data[gene_data['transcript_type'] == 'retained_intron']
        ri_proportion = ri_data['biotype_expression_proportion'].iloc[0] if len(ri_data) > 0 else 0
        
        # Calculate ratio (handle division by zero)
        if ri_proportion > 0:
            ratio = pc_proportion / ri_proportion
        elif pc_proportion > 0:
            ratio = float('inf')  # Only protein_coding, no retained_intron
        else:
            ratio = 0  # Neither type present
        
        gene_ratios.append({
            'gene_id': gene_id,
            'gene_name': gene_name,
            'total_gene_expression': total_gene_expression,
            'protein_coding_proportion': pc_proportion,
            'retained_intron_proportion': ri_proportion,
            'pc_to_ri_ratio': ratio
        })
    
    # Convert to DataFrame and sort by ratio (descending)
    ratio_df = pd.DataFrame(gene_ratios)
    ratio_df = ratio_df.sort_values('pc_to_ri_ratio', ascending=False)
    
    # Add rank column
    ratio_df['rank'] = range(1, len(ratio_df) + 1)
    
    return ratio_df

def create_violin_plot(expression_df, outdir, highlight_goi='SRA1'):
    transcript_types = ['protein_coding', 'retained_intron']
    colors = {
        'protein_coding': '#1f77b4',
        'retained_intron': '#ff7f0e'
    }

    # Filter and simplify data
    data = expression_df[
        expression_df['transcript_type'].isin(transcript_types)
    ].drop_duplicates(subset=['gene_id', 'transcript_type'])

    data = data[['gene_id', 'gene_name', 'transcript_type', 'biotype_expression_proportion', 'total_gene_expression']]

    # Bin for jitter layer
    # Get rows for GOI
    goi_rows = expression_df.loc[expression_df['gene_name'] == highlight_goi]

    if goi_rows.empty:
        raise ValueError(f"Gene of interest '{highlight_goi}' not found in expression_df.")

    # Get the unique expression bin for the GOI
    goi_bins = goi_rows['expression_bin'].unique()
    if len(goi_bins) != 1:
        raise ValueError(f"GOI '{highlight_goi}' is assigned to multiple or no expression bins: {goi_bins}")

    goi_bin = goi_bins[0]  # e.g., "40-50"

    # Now filter all genes in the same bin
    bin_data = expression_df[
        (expression_df['transcript_type'].isin(transcript_types)) &
        (expression_df['expression_bin'] == goi_bin)]

    # Create the plot
    fig, ax = plt.subplots(figsize=(8, 6))

    for i, ttype in enumerate(transcript_types):
        subset = data[data['transcript_type'] == ttype].copy()
        subset['x'] = i

        sns.violinplot(
            y='biotype_expression_proportion',
            x='x',
            data=subset,
            cut=0,
            bw=0.15,
            width=0.8,
            inner=None,
            color=colors[ttype],
            linewidth=1.0,
            ax=ax
        )

        # Mask the right half
        ax.add_patch(plt.Rectangle(
            (i, -0.05), 0.5, 1.15,
            facecolor='white',
            zorder=3
        ))


        # IQR box + median line
        q1 = subset['biotype_expression_proportion'].quantile(0.25)
        q3 = subset['biotype_expression_proportion'].quantile(0.75)
        median = subset['biotype_expression_proportion'].median()

        box_width = 0.05
        ax.add_patch(plt.Rectangle(
            (i - box_width, q1),
            box_width,
            q3 - q1,
            facecolor='gray',
            edgecolor='black',
            linewidth=1.5,
            zorder=4
        ))

        ax.plot(
            [i - box_width, i],
            [median, median],
            color='black',
            linewidth=2,
            zorder=5
        )

        # Strip plot for 100–200 bin (jittered dots)
        strip_subset = bin_data[bin_data['transcript_type'] == ttype]
        jitter = np.random.uniform(0.02, 0.04, size=len(strip_subset))

        for j, (_, row) in enumerate(strip_subset.iterrows()):
            is_goi = row['gene_name'] == highlight_goi
            ax.scatter(
                i + jitter[j],
                row['biotype_expression_proportion'],
                marker='*' if is_goi else 'o',
                s=120 if is_goi else 15,
                edgecolor='black',
                linewidth=1.5 if is_goi else 0.5,
                color='gold' if is_goi else colors[ttype],
                zorder=6 if is_goi else 5,
                label=f"{highlight_goi} ({ttype})" if is_goi else None
            )

    # Customize axes
    ax.set_xticks(range(len(transcript_types)))
    ax.set_xticklabels(transcript_types)
    ax.set_xlim(-0.5, len(transcript_types) - 0.5)
    ax.set_ylim(-0.05, 1.05)  # Give visual breathing room
    ax.set_ylabel("Biotype Expression Proportion")
    ax.set_title("Half-Violin Plot of Expression Proportions")

    # Remove duplicate GOI legend entries
    handles, labels = ax.get_legend_handles_labels()
    seen = set()
    deduped = [(h, l) for h, l in zip(handles, labels) if l and not (l in seen or seen.add(l))]
    if deduped:
        ax.legend(*zip(*deduped), loc='best')

    plt.tight_layout()
    plt.savefig(f"{outdir}/{highlight_goi}_seaborn_half_violin.png", dpi=300)
    plt.savefig(f"{outdir}/{highlight_goi}_seaborn_half_violin.svg")
    plt.savefig(f"{outdir}/{highlight_goi}_seaborn_half_violin.pdf")
    plt.close()

    print(f"✓ Seaborn half-violin plot saved to: {outdir}/{highlight_goi}_seaborn_half_violin.*")

def main():
    # Main script execution
    gtf_data = load_transcripts_from_gtf(gtf_file, num_workers=10)
    kal_data = parse_kallisto(kallisto_file, num_workers=20)

    # Merge expression data with correct isoform binning
    expression_proportions = merge_expression(kal_data, gtf_data)
    weird_expression = expression_proportions.loc[(expression_proportions["transcript_type"] == "retained_intron") & (expression_proportions["biotype_expression_proportion"] > 0.2)]
    weird_expression.to_csv("weird_expression.tsv", sep="\t", header=True)

    num_isoforms_per_gene = expression_proportions.groupby("gene_id")["total_isoforms"].first()
    expression_per_gene = expression_proportions.groupby("gene_id")["total_gene_expression"].first()

    # Comment out the qcut code
    # expression_bins, expression_bin_labels = pd.qcut(expression_per_gene, q=10, retbins=True, duplicates="drop", precision=0)
    
    # Define custom bin edges
    bin_edges_0_100 = list(range(0, 101, 10))  # [0, 10, 20, ..., 100]
    bin_edges_100_500 = list(range(100, 501, 100))  # [100, 150, 200, ..., 500]
    custom_bin_edges = sorted(list(set(bin_edges_0_100 + bin_edges_100_500 + [float('inf')])))  # Combine and add 500+

    # Create human-readable bin labels
    custom_bin_labels = []
    for i in range(len(custom_bin_edges) - 1):
        if custom_bin_edges[i+1] == float('inf'):
            custom_bin_labels.append(f"{custom_bin_edges[i]}+")
        else:
            custom_bin_labels.append(f"{custom_bin_edges[i]}-{custom_bin_edges[i+1]}")

    # Use pd.cut to create bins with custom edges
    expression_bins = pd.cut(expression_per_gene, bins=custom_bin_edges, labels=custom_bin_labels, right=False)

    # Create a DataFrame mapping `gene_id` to its expression bin
    bin_mapping = pd.DataFrame({
        "gene_id": expression_per_gene.index, 
        "expression_bin": expression_bins
    }).reset_index(drop=True)

    # Create a DataFrame for expression bin counts
    expression_bin_counts = pd.value_counts(bin_mapping["expression_bin"]).reset_index(name="num_genes_in_expression_bin")

    # Merge bin mapping into expression_proportions
    expression_proportions = expression_proportions.merge(bin_mapping, on="gene_id", how="left")

    # Merge bin counts into expression_proportions
    expression_proportions = expression_proportions.merge(expression_bin_counts, on="expression_bin", how="left")
    
    # Add transcript type stats
    expression_proportions = add_transcript_type_stats(expression_proportions).dropna()
    
    # Add percentile rank of each transcript's expression proportion within its bin
    expression_proportions["biotype_expr_percentile_in_bin"] = (
        expression_proportions
        .groupby("expression_bin")["biotype_expression_proportion"]
        .transform(lambda x: x.rank(pct=True))
    )


    # Debugging: Print first few rows to verify
    print(expression_proportions[["gene_id", "expression_bin", "num_genes_in_expression_bin"]].head().to_string())
    print(f"Number of bins created: {len(custom_bin_labels)}")
    print(f"Bin labels: {custom_bin_labels}")
    expression_proportions.to_csv(f"{outdir}/expression_proportions.csv", index=False)

    # Calculate protein_coding to retained_intron ratio and save results
    gene_ratios = calculate_protein_coding_to_retained_intron_ratio(expression_proportions)
    
    # Save the ranked genes to a file
    output_file = f"{outdir}/genes_ranked_by_pc_to_ri_ratio.tsv"
    gene_ratios.to_csv(output_file, sep="\t", index=False)
    print(f"Saved gene rankings to: {output_file}")
    
    # Print top 10 genes with highest ratio
    print("\nTop 10 genes with highest protein_coding:retained_intron ratio:")
    print(gene_ratios.head(10)[['rank', 'gene_name', 'pc_to_ri_ratio', 'protein_coding_proportion', 'retained_intron_proportion']].to_string(index=False))
    
    # Print top 10 genes with lowest ratio (excluding infinite values)
    finite_ratios = gene_ratios[gene_ratios['pc_to_ri_ratio'] != float('inf')]
    if len(finite_ratios) > 0:
        print("\nTop 10 genes with lowest protein_coding:retained_intron ratio:")
        print(finite_ratios.tail(10)[['rank', 'gene_name', 'pc_to_ri_ratio', 'protein_coding_proportion', 'retained_intron_proportion']].to_string(index=False))

    # Pass the bin labels to the plot function
    create_combined_plot(expression_proportions, [custom_bin_labels, custom_bin_labels], outdir, highlight_goi)
    create_violin_plot(expression_proportions, outdir, highlight_goi)

if __name__ == "__main__":
    main()
