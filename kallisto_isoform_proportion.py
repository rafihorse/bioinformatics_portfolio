#!/usr/bin/env python
import pandas as pd
import argparse
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import matplotlib.colors as mcolors
import functools
from functools import reduce
from concurrent.futures import ProcessPoolExecutor

# Define the theme for Plotly figures
theme = 'plotly_white'

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--input_list', help="List of TPM files")
parser.add_argument('-t', '--gtf', help="GTF file with transcript annotations")
parser.add_argument('-g', '--gene_list', help="File containing list of genes to analyze")
parser.add_argument('--topN', type=int, help="If set, will only take top N transcripts per gene")
parser.add_argument('--min_exp', type=float, help="If set, will only include transcripts with mean expression above this threshold")
parser.add_argument('--outfiles', type=str, help="Custom naming convention for output files")
args = parser.parse_args()

def read_tpm_file(file_line):
    acc, tissue, db, work, f, n_reads, n_aligned, p_aligned = file_line.strip().split("\t")
    dataset = f"{acc}+{tissue.replace(' ', '.')}+{db}"
    df = pd.read_csv(f, sep="\t", usecols=["target_id", "tpm"], na_filter=False)
    df["tpm"] += 1
    df.rename(columns={'tpm': dataset, 'target_id': 'transcript_id'}, inplace=True)
    return df

def get_transcript_ids_for_genes(gtf_file, gene_list):
    """Extract transcript IDs for the genes in the gene list."""
    transcripts = load_transcripts_from_gtf(gtf_file)
    gene_transcripts = transcripts[transcripts['gene_name'].isin(gene_list)]
    return set(gene_transcripts['transcript_id'])

def read_filtered_tpm_file(file_line, transcript_ids):
    """Read and filter a single TPM file for relevant transcript IDs."""
    acc, tissue, db, work, f, n_reads, n_aligned, p_aligned = file_line.strip().split("\t")
    dataset = f"{acc}+{tissue.replace(' ', '.')}+{db}"
    df = pd.read_csv(f, sep="\t", usecols=["target_id", "tpm"], na_filter=False)
    df = df[df["target_id"].isin(transcript_ids)]  # Subset to relevant transcript IDs
    df["tpm"] += 1
    df.rename(columns={'tpm': dataset, 'target_id': 'transcript_id'}, inplace=True)
    return df

def load_tpm_data(input_list, transcript_ids, num_workers=60, chunk_size=50):
    """Load TPM data efficiently by processing files in chunks and merging on transcript_id."""
    with open(input_list, 'r') as files:
        file_lines = files.readlines()
    
    # Process files in chunks to reduce memory overhead
    all_dfs = []
    for i in range(0, len(file_lines), chunk_size):
        chunk = file_lines[i:i + chunk_size]
        with ProcessPoolExecutor(max_workers=num_workers) as executor:
            dfs = list(executor.map(functools.partial(read_filtered_tpm_file, transcript_ids=transcript_ids), chunk))
        # Merge all DataFrames in the chunk on 'transcript_id'
        merged_chunk = functools.reduce(lambda left, right: pd.merge(left, right, on="transcript_id", how="outer"), dfs)
        all_dfs.append(merged_chunk)

    # Merge all chunks on 'transcript_id'
    tpm = functools.reduce(lambda left, right: pd.merge(left, right, on="transcript_id", how="outer"), all_dfs)
    # tpm.to_csv("tpm_combined.tsv", index=False, header=True, sep="\t")
    tpm["mean_exp"] = tpm.iloc[:, 1:].mean(axis=1)  # Calculate mean expression across datasets
    return tpm

def parse_gtf_attributes(attr_str):
    attrs = {}
    for attr in attr_str.strip().split(';'):
        if attr.strip():
            try:
                # Find the first space or equals sign
                parts = attr.strip().split('=', 1)
                if len(parts) == 1:  # No equals sign found
                    parts = attr.strip().split(None, 1)
                
                if len(parts) == 2:
                    key, value = parts
                    # Remove quotes if present
                    attrs[key] = value.strip('"\'')
            except ValueError:
                continue
    return attrs

def load_transcripts_from_gtf(gtf_file):
    transcripts = []
    with open(gtf_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9 or fields[2] != 'transcript':
                continue
                
            attrs = parse_gtf_attributes(fields[8])
            transcript_info = {
                'transcript_id': attrs.get('transcript_id', ''),
                'transcript_name': attrs.get('transcript_name', 'unknown'),
                'transcript_type': attrs.get('transcript_type', ''),
                'gene_name': attrs.get('gene_name', 'unknown'),
                'gene_id': attrs.get('gene_id', ''),
                'transcript_source': attrs.get('transcript_source', ''),
                'transcript_id_version': attrs.get('transcript_id', ''),
            }
            
            # Add version if available, but avoid double extension
            if 'transcript_version' in attrs and '.' not in transcript_info['transcript_id_version']:
                transcript_info['transcript_id_version'] += '.' + attrs.get('transcript_version', '1')
            
            if not transcript_info['transcript_type']:
                transcript_info['transcript_type'] = attrs.get('transcript_biotype', 'unknown')
                
            if transcript_info['transcript_source'] == 'serafina_manual':
                transcript_info['transcript_id_version'] = attrs.get('transcript_id', '')
                
            transcripts.append(transcript_info)
    
    return pd.DataFrame(transcripts)

def process_gene(args):
    gene_name, tpm_data, transcript_info, topN, min_exp, file_name = args
    
    # Filter for gene's transcripts
    gene_transcripts = transcript_info[transcript_info['gene_name'] == gene_name]
    print("HELLO")
    print(gene_transcripts)
    # gene_transcripts.to_csv(f"{gene_name}_transcript_info.tsv", index=False, header=True, sep="\t")

    matching_transcripts = gene_transcripts['transcript_id'].unique()
    print(matching_transcripts)
    filter_tpm = tpm_data[tpm_data["transcript_id"].isin(matching_transcripts)]

    # filter_tpm.to_csv(f"{gene_name}_expression_data_raw.tsv", index=False, header=True, sep="\t")
    
    if filter_tpm.empty:
        print(f"No data found for gene {gene_name}")
        return None
        
    filter_tpm = filter_tpm.sort_values(by="mean_exp", ascending=False)
    
    if min_exp is not None:
        filter_tpm = filter_tpm[filter_tpm["mean_exp"] >= min_exp]
        if filter_tpm.empty:
            print(f"No transcripts with mean expression >= {min_exp} for gene {gene_name}")
            return None
    
    if topN:
        filter_tpm = filter_tpm.head(topN)
    
    # Sort by major transcript
    major_transcript_index = filter_tpm["mean_exp"].idxmax()
    columns_to_include = filter_tpm.columns[1:-1]
    major_sort_order = filter_tpm.loc[major_transcript_index, columns_to_include].sort_values(ascending=True).index
    sorted_data = pd.concat([filter_tpm.iloc[:, 0], filter_tpm[major_sort_order], filter_tpm.iloc[:, -1:]], axis=1)

    # sorted_data.to_csv(f"{gene_name}_expression_data_unmelted.tsv", index=False, header=True, sep="\t")
    
    # Melt and prepare for plotting
    id_vars = ["transcript_id", "mean_exp"]
    value_vars = [col for col in sorted_data.columns if col not in id_vars]
    melted = pd.melt(sorted_data, id_vars=id_vars, value_vars=value_vars, var_name="dataset", value_name="tpm")
    melted[["accession", "tissue", "dbxref"]] = melted['dataset'].str.split("+", n=2, expand=True)
    
    # Format tissue names: replace periods with spaces and capitalize first letter
    melted['tissue'] = (melted['tissue']
                        .str.replace('.', ' ')
                        .str.replace('-positive', '+', regex=False)
                        .str.replace('alpha-beta', 'αβ', regex=False)
                        .str.capitalize())
    
    # melted.to_csv(f"{gene_name}_expression_data_melted.tsv", index=False, header=True, sep="\t")
    
    # Merge with transcript info
    full_data = melted.merge(gene_transcripts, left_on='transcript_id', right_on='transcript_id', how='left')

    full_data.to_csv(f"{gene_name}_{file_name}.tsv", index=False, header=True, sep="\t")
    full_data.rename(columns={'transcript_id_x': 'transcript_id'}, inplace=True)
    
    # Call all plot functions using full_data
    res_plot = create_plot(full_data, gene_name, file_name)
    res_stacked = plot_stacked_bar(full_data, gene_name, file_name)
    res_combined = create_combined_plot(full_data, gene_name, file_name)
    
    # Return full_data instead of just a boolean flag
    return full_data if (res_plot and res_stacked and res_combined) else None

def plot_stacked_bar(full_data, gene_name, file_name):
    # Subset data for the gene of interest
    full_data = full_data[full_data['gene_name'] == gene_name]
    full_data.drop_duplicates(inplace=True)
    # Group by tissue and transcript_type and average TPM values across datasets
    grouped = full_data.groupby(['tissue', 'transcript_id', 'transcript_type'])['tpm'].mean().reset_index()
    
    # Now sum by transcript type within each tissue
    grouped = grouped.groupby(['tissue', 'transcript_type'])['tpm'].sum().unstack(fill_value=0)
    
    if 'protein_coding' not in grouped.columns:
        grouped['protein_coding'] = 0
    
    other_cols = [col for col in grouped.columns if col != 'protein_coding']
    grouped['other'] = grouped[other_cols].sum(axis=1)
    
    grouped['total'] = grouped['protein_coding'] + grouped['other']
    sorted_tissues = grouped.sort_values('total', ascending=True).index
    
    fig = go.Figure()
    
    fig.add_trace(go.Bar(
        name='Protein-coding Transcripts',
        x=list(sorted_tissues),
        y=grouped.loc[sorted_tissues, 'protein_coding'],
        marker_color='cornflowerblue'
    ))
    
    fig.add_trace(go.Bar(
        name='Non-coding Transcripts',
        x=list(sorted_tissues),
        y=grouped.loc[sorted_tissues, 'other'],
        marker_color='grey'
    ))
    
    fig.update_layout(
        barmode='stack',
        title=f"{gene_name} Protein-coding vs Non-coding Transcripts",
        xaxis_title="Tissue / Cell Line",
        yaxis_title="Average Expression (TPM)",
        xaxis_tickangle=45,
        width=1200,
        height=600,
        template=theme,
        legend=dict(
            x=1.05,  # Move legend further to the right
            y=1,     # Keep legend at the top
            xanchor='left',  # Anchor legend to the left
            yanchor='top'    # Anchor legend to the top
        )
    )
    
    fig.write_html(f'{gene_name}_stacked_nc_c.html')
    fig.write_image(f'{gene_name}_stacked_nc_c.png', scale=2)
    fig.write_image(f'{gene_name}_stacked_nc_c.svg', scale=2)  # Added SVG output
    return True

def create_plot(full_data, gene_name, file_name):
    style_map = {
        "protein_coding": "circle-dot",
        "retained_intron": "square-dot",
        "protein_coding_CDS_not_defined": "cross-dot",
        "lncRNA": "diamond-dot",
        "nonsense_mediated_decay": "x-thin-open",
        "processed_transcript": "asterisk-open"
    }
    
    # Identify the major transcript (highest mean expression)
    major_transcript = full_data.loc[full_data['mean_exp'].idxmax()]['transcript_id']
    major_transcript_type = full_data.loc[full_data['mean_exp'].idxmax()]['transcript_type']
    
    fig = go.Figure()
    plotted_transcript_types = set()
    
    # Sort tissues by major transcript expression
    print(gene_name, major_transcript)

    reference_order = (full_data[full_data['transcript_id'] == major_transcript]
                      .groupby('tissue')['tpm']
                      .mean()
                      .sort_values(ascending=True)
                      .index)

    # Keep the original reference order for plotting
    full_data['tissue'] = pd.Categorical(full_data['tissue'], categories=reference_order, ordered=True)

    # Create a mapping for updated x-axis labels
    x_labels = [tissue.replace('-positive', '+').replace('alpha-beta', 'αβ') for tissue in reference_order]

    # Calculate max TPM for y-axis range based on the highest mean TPM across tissues for the major transcript
    major_transcript_group = full_data[full_data['transcript_id'] == major_transcript].groupby("tissue")['tpm'].mean()
    max_tpm = np.ceil(major_transcript_group.max() * 1.1) if not major_transcript_group.empty else 1  # Add 10% buffer

    # Set tick interval to 10 or 1 based on max_tpm
    tick_interval = 10 if max_tpm > 10 else 1

    fig.update_layout(
        width=2000,
        height=800,
        yaxis=dict(
            range=[0, max_tpm],
            title='Expression (TPM)',
            tickmode='linear',  # Use linear ticks
            dtick=tick_interval,  # Set consistent tick intervals
            showgrid=True,  # Show gridlines
            gridcolor='lightgrey',  # Light grey gridlines for better visibility
            title_standoff=10  # Add space between axis title and ticks
        ),
        xaxis=dict(
            title='Tissue / Cell Line',
            tickangle=60,  # Increase angle to prevent overlap
            ticktext=x_labels,  # Use updated labels for display
            tickvals=reference_order  # Use original values for plotting
        ),
        title=f'{gene_name} Isoform-Level Transcript Expression',
        template=theme,
        legend=dict(
            x=1.05,  # Move legend further to the right
            y=1,     # Keep legend at the top
            xanchor='left',  # Anchor legend to the left
            yanchor='top'    # Anchor legend to the top
        )
    )
    
    for (transcript_name, transcript_type), group in full_data.groupby(['transcript_name', 'transcript_type']):
        group = group.sort_values(by='tissue')
        mean_tpm = group.groupby("tissue")['tpm'].mean()
        se_tpm = group.groupby("tissue")['tpm'].sem()
        
        # Get the transcript ID for this group
        transcript_id = group['transcript_id'].iloc[0]
        is_major = transcript_id == major_transcript
        
        # Determine color based on transcript type and whether it's the major transcript
        if transcript_type == "protein_coding":
            color = "darkblue" if is_major else "cornflowerblue"
        elif transcript_type == "retained_intron":
            color = "red" if is_major else "lightcoral"
        else:
            color = "grey"
            
        symbol = style_map.get(transcript_type, 'circle-dot')

        # Only add traces if there is data to plot
        if not mean_tpm.empty:
            plotted_transcript_types.add(transcript_type)
            fig.add_trace(go.Scatter(
                x=mean_tpm.index,
                y=mean_tpm.values,
                error_y=dict(type='data', array=se_tpm.values, visible=True),
                mode='lines+markers',
                marker=dict(symbol=symbol, size=10, color=color),
                line=dict(color=color, width=2),
                name=transcript_name,
                legendgroup='Transcript Names',
                legendgrouptitle_text="Transcript Name",
                showlegend=True,
                meta=[transcript_name, transcript_type],
                hovertemplate="Transcript: %{meta[0]}<br>Tissue: %{x}<br>TPM: %{y:.2f}<extra></extra>",
            ))
    
    # Add transcript type legend items only if they were plotted
    if "protein_coding" in plotted_transcript_types:
        fig.add_trace(go.Scatter(
            x=[None], y=[None], mode='markers',
            marker=dict(symbol=style_map.get("protein_coding"), size=10, color="darkblue"),
            name="Protein-coding (major)",
            legendgroup='Transcript Types', legendgrouptitle=dict(text="Transcript Type"),
            showlegend=True, hoverinfo='skip'
        ))
        fig.add_trace(go.Scatter(
            x=[None], y=[None], mode='markers',
            marker=dict(symbol=style_map.get("protein_coding"), size=10, color="cornflowerblue"),
            name="Protein-coding",
            legendgroup='Transcript Types', legendgrouptitle=dict(text="Transcript Type"),
            showlegend=True, hoverinfo='skip'
        ))
    
    if "retained_intron" in plotted_transcript_types:
        if major_transcript_type == "retained_intron":  # Only add (major) if the major transcript is retained_intron
            fig.add_trace(go.Scatter(
                x=[None], y=[None], mode='markers',
                marker=dict(symbol=style_map.get("retained_intron"), size=10, color="red"),
                name="Retained intron (major)",
                legendgroup='Transcript Types', legendgrouptitle=dict(text="Transcript Type"),
                showlegend=True, hoverinfo='skip'
            ))
        fig.add_trace(go.Scatter(
            x=[None], y=[None], mode='markers',
            marker=dict(symbol=style_map.get("retained_intron"), size=10, color="lightcoral"),
            name="Retained intron",
            legendgroup='Transcript Types', legendgrouptitle=dict(text="Transcript Type"),
            showlegend=True, hoverinfo='skip'
        ))
    
    for ttype in [t for t in style_map.keys() if t not in ["protein_coding", "retained_intron"]]:
        if ttype in plotted_transcript_types:
            symbol = style_map.get(ttype, 'circle-dot')
            formatted_ttype = ttype.replace('_', ' ').replace('cds', 'CDS').capitalize()
            fig.add_trace(go.Scatter(
                x=[None], y=[None], mode='markers',
                marker=dict(symbol=symbol, size=10, color="grey"),
                name=formatted_ttype, legendgroup='Transcript Types',
                legendgrouptitle=dict(text="Transcript Type"),
                showlegend=True, hoverinfo='skip'
            ))
    
    # Print unique tissues and cell types
    unique_tissues = full_data['tissue'].unique()
    print(f"Unique tissues and cell types for {gene_name}:")
    print(unique_tissues)
    print(f"Total unique tissues and cell types: {len(unique_tissues)}")

    # Calculate mean proportion of SRA1 expression by all non-protein-coding transcripts
    if gene_name == "SRA1":
        non_protein_coding = full_data[full_data['transcript_type'] != 'protein_coding']
        total_tpm = full_data['tpm'].sum()
        non_protein_coding_tpm = non_protein_coding['tpm'].sum()
        mean_proportion = (non_protein_coding_tpm / total_tpm) if total_tpm > 0 else 0
        print(f"Mean proportion of SRA1 expression by non-protein-coding transcripts: {mean_proportion:.2%}")

    fig.write_html(f'{gene_name}_{file_name}_plot.html')
    fig.write_image(f'{gene_name}_{file_name}_plot.png', scale=2)
    fig.write_image(f'{gene_name}_{file_name}_plot.svg', scale=2)  # Added SVG output
    return True

def create_combined_plot(full_data, gene_name, file_name):
    style_map = {
        "protein_coding": "circle-dot",
        "retained_intron": "square-dot",
        "protein_coding_CDS_not_defined": "cross-dot",
        "lncRNA": "cross-dot",
        "nonsense_mediated_decay": "cross-dot",
        "processed_transcript": "cross-dot"
    }

    # Identify the major transcript (highest mean expression)
    major_transcript = full_data.loc[full_data['mean_exp'].idxmax()]['transcript_id']
    major_transcript_name = full_data.loc[full_data['mean_exp'].idxmax()]['transcript_name']

    # Sort tissues by major transcript expression
    reference_order = (full_data[full_data['transcript_id'] == major_transcript]
                      .groupby('tissue')['tpm']
                      .mean()
                      .sort_values(ascending=True)
                      .index)

    full_data['tissue'] = pd.Categorical(full_data['tissue'], categories=reference_order, ordered=True)

    # Create a mapping for updated x-axis labels
    x_labels = [tissue.replace('-positive', '+').replace('alpha-beta', 'αβ') for tissue in reference_order]

    # Prepare data for the stacked bar plot
    grouped = full_data.groupby(['tissue', 'transcript_id', 'transcript_type'])['tpm'].mean().reset_index()
    grouped = grouped.groupby(['tissue', 'transcript_type'])['tpm'].sum().unstack(fill_value=0)

    if 'protein_coding' not in grouped.columns:
        grouped['protein_coding'] = 0

    other_cols = [col for col in grouped.columns if col != 'protein_coding']
    grouped['other'] = grouped[other_cols].sum(axis=1)

    # Calculate proportions
    grouped['total'] = grouped['protein_coding'] + grouped['other']
    grouped['protein_coding_prop'] = grouped['protein_coding'] / grouped['total']
    grouped['other_prop'] = grouped['other'] / grouped['total']

    # Create the combined plot
    fig = go.Figure()

    # Add the stacked bar plot (proportions)
    fig.add_trace(go.Bar(
        name='Proportion of Protein-coding Transcripts',
        x=list(reference_order),
        y=grouped.loc[reference_order, 'protein_coding_prop'],
        marker_color='royalblue',
        opacity=0.5,
        offsetgroup=1,
        xaxis='x',
        yaxis='y1'
    ))

    fig.add_trace(go.Bar(
        name='Proportion of Non-coding Transcripts',
        x=list(reference_order),
        y=grouped.loc[reference_order, 'other_prop'],
        marker_color='black',
        opacity=0.5,
        offsetgroup=1,
        xaxis='x',
        yaxis='y1'
    ))

    # Add the expression plot (log scale, plotted on the right y-axis)
    for (transcript_name, transcript_type), group in full_data.groupby(['transcript_name', 'transcript_type']):
        group = group.sort_values(by='tissue')
        mean_tpm = group.groupby("tissue")['tpm'].mean()
        se_tpm = group.groupby("tissue")['tpm'].sem()

        # Calculate log_max based on the maximum TPM value
        max_tpm = np.ceil(group['tpm'].max()) if group['tpm'].max() > 0 else 1

        transcript_id = group['transcript_id'].iloc[0]
        is_major = transcript_id == major_transcript

        if transcript_type == "protein_coding":
            color = "darkblue" if is_major else "royalblue"
        elif transcript_type == "retained_intron":
            color = "red" if is_major else "firebrick"
        else:
            color = "dimgrey"

        symbol = style_map.get(transcript_type, 'circle-dot')

        if not mean_tpm.empty:
            fig.add_trace(go.Scatter(
                x=mean_tpm.index,
                y=mean_tpm.values,
                error_y=dict(type='data', array=se_tpm.values, visible=True),
                mode='lines+markers',
                marker=dict(symbol=symbol, size=8, color=color),
                line=dict(color=color, width=2),
                name=None,  # Do not show transcript names in the legend
                legendgroup='Transcript Types',
                showlegend=False,
                meta=[transcript_name, transcript_type],
                hovertemplate="Transcript: %{meta[0]}<br>Tissue: %{x}<br>TPM: %{y:.2f}<extra></extra>",
                xaxis='x',
                yaxis='y2'
            ))

    # Add legend items for transcript types
    fig.add_trace(go.Scatter(
        x=[None], y=[None], mode='markers',
        marker=dict(symbol=style_map.get("protein_coding"), size=10, color="darkblue"),
        name=f"Major Transcript ({major_transcript_name})",
        legendgroup='Transcript Types',
        showlegend=True, hoverinfo='skip'
    ))
    fig.add_trace(go.Scatter(
        x=[None], y=[None], mode='markers',
        marker=dict(symbol=style_map.get("protein_coding"), size=10, color="royalblue"),
        name="Protein-coding",
        legendgroup='Transcript Types',
        showlegend=True, hoverinfo='skip'
    ))
    fig.add_trace(go.Scatter(
        x=[None], y=[None], mode='markers',
        marker=dict(symbol=style_map.get("retained_intron"), size=10, color="firebrick"),
        name="Retained intron",
        legendgroup='Transcript Types',
        showlegend=True, hoverinfo='skip'
    ))
    fig.add_trace(go.Scatter(
        x=[None], y=[None], mode='markers',
        marker=dict(symbol=style_map.get("lncRNA"), size=10, color="dimgrey"),
        name="Other non-coding",
        legendgroup='Transcript Types',
        showlegend=True, hoverinfo='skip'
    ))

    # Update layout for combined plot
    fig.update_layout(
        barmode='stack',
        width=1600,
        height=650,
        title=f'{gene_name} Combined Isoform-Level Transcript Expression and Stacked Bar Plot',
        xaxis=dict(
            title='Tissue / Cell Line',
            tickangle=90,
            ticktext=x_labels,
            tickvals=reference_order,
            domain=[0, 1],
            title_standoff=20
        ),
        yaxis=dict(
            title='Proportion',
            range=[0, 1],
            side='right',
            titlefont=dict(color='black'),
            tickfont=dict(color='black'),
            anchor='x'
        ),
        yaxis2=dict(
            title='Isoform Expression (TPM)',
            range=[0, max_tpm],
            side='left',
            overlaying='y',
            titlefont=dict(color='grey'),
            tickfont=dict(color='grey'),
            anchor='x',
            minor_showgrid=False,
            dtick=2
        ),
        template=theme,
        legend=dict(
            x=1.05,
            y=1,
            xanchor='left',
            yanchor='top'
        )
    )

    fig.write_html(f'{gene_name}_combined_plot.html')
    fig.write_image(f'{gene_name}_combined_plot.png', scale=2)
    fig.write_image(f'{gene_name}_combined_plot.svg', scale=2)
    return True

def plot_genes_comparison(all_genes_data, output_name="genes_comparison"):
    """
    Plot protein coding vs non-coding TPM for each gene (averaged across all datasets).
    Each x-axis position represents a different gene, arranged by ratio of non-coding to protein-coding TPM.
    
    Args:
        all_genes_data: Dictionary of DataFrames with gene_name as key and full_data as value
        output_name: Base name for output files
    """
    gene_data = []
    
    for gene_name, full_data in all_genes_data.items():
        if full_data is None or full_data.empty:
            continue
            
        # Get dataset columns
        datasets = full_data['dataset'].unique()
        
        # Group first by transcript type and dataset, summing across all transcripts of that type
        grouped = full_data.groupby(['transcript_type', 'dataset'])['tpm'].sum().reset_index()
        
        # Now compute the mean across all datasets for each transcript type
        type_means = grouped.groupby('transcript_type')['tpm'].mean()
        
        # Get protein coding TPM (sum of all protein coding transcripts, averaged across datasets)
        protein_coding_tpm = type_means.get('protein_coding', 0)
        
        # Sum all non-protein-coding transcript types
        non_coding_types = [t for t in type_means.index if t != 'protein_coding']
        non_coding_tpm = sum(type_means.get(t, 0) for t in non_coding_types)
        
        # Calculate ratio (handle division by zero)
        if protein_coding_tpm == 0:
            ratio = float('inf')  # Infinity for genes with no protein-coding expression
        else:
            ratio = non_coding_tpm / protein_coding_tpm
        
        # Calculate proportions
        total_tpm = protein_coding_tpm + non_coding_tpm
        protein_coding_prop = protein_coding_tpm / total_tpm if total_tpm > 0 else 0
        non_coding_prop = non_coding_tpm / total_tpm if total_tpm > 0 else 0
        
        gene_data.append((gene_name, protein_coding_tpm, non_coding_tpm, ratio, protein_coding_prop, non_coding_prop))
    
    if not gene_data:
        print("No valid gene data for comparison plot")
        return False
        
    # Sort genes by ratio of non-coding to protein-coding TPM (descending order)
    gene_data.sort(key=lambda x: x[3], reverse=True)  # Sort by ratio (4th element) in descending order
    
    gene_names = [data[0] for data in gene_data]
    protein_coding_tpm = [data[1] for data in gene_data]
    non_coding_tpm = [data[2] for data in gene_data]
    ratios = [data[3] for data in gene_data]
    protein_coding_props = [data[4] for data in gene_data]
    non_coding_props = [data[5] for data in gene_data]
    
    # Create plotly figure
    from plotly.subplots import make_subplots
    fig = make_subplots(rows=2, cols=1, shared_xaxes=True, 
                       vertical_spacing=0.1,
                       row_heights=[0.7, 0.3])  # Allocate 70% to expression, 30% to proportions
    
    # Add protein coding TPM bars to top subplot
    fig.add_trace(go.Bar(
        x=gene_names,
        y=protein_coding_tpm,
        name='Protein-coding TPM',
        marker_color='cornflowerblue',
        offsetgroup=0
    ), row=1, col=1)
    
    # Add non-coding TPM bars to top subplot
    fig.add_trace(go.Bar(
        x=gene_names,
        y=non_coding_tpm,
        name='Non-coding TPM',
        marker_color='grey',
        offsetgroup=1
    ), row=1, col=1)
    
    # Add stacked proportion bars to bottom subplot
    fig.add_trace(go.Bar(
        x=gene_names,
        y=protein_coding_props,
        name='Proportion Protein-coding',
        marker_color='royalblue',
        width=0.8,  # Make bars wider
        showlegend=True
    ), row=2, col=1)
    
    fig.add_trace(go.Bar(
        x=gene_names,
        y=non_coding_props,
        name='Proportion Non-coding',
        marker_color='lightgrey',
        width=0.8,  # Make bars wider
        showlegend=True
    ), row=2, col=1)
    
    # Update layout
    fig.update_layout(
        title="Protein-coding vs Non-coding Expression by Gene (Sorted by NC:PC Ratio)",
        width=max(1200, len(gene_names) * 100),  # Adjust width based on gene count
        height=800,  # Increased height for two subplots
        template=theme,
        legend=dict(
            x=1.05,
            y=1,
            xanchor='left',
            yanchor='top'
        )
    )
    
    # Update subplot-specific settings
    fig.update_layout(barmode='group')  # For top subplot
    
    # Apply stacked mode specifically to the second row
    fig.update_layout(
        barmode='stack',
        bargap=0.15
    )
    
    # Update axes
    fig.update_yaxes(title="Average Expression (TPM)", row=1, col=1)
    fig.update_yaxes(title="Proportion", range=[0, 1], tickformat=".0%", row=2, col=1)
    fig.update_xaxes(title="Gene (arranged by non-coding : protein-coding ratio)", tickangle=45, row=2, col=1)
    
    # Add annotation about sorting
    fig.add_annotation(
        text="Genes arranged by decreasing ratio of non-coding to protein-coding expression",
        xref="paper", yref="paper",
        x=0.5, y=1.05,
        showarrow=False,
        font=dict(size=12)
    )
    
    # Save figures
    fig.write_html(f'{output_name}.html')
    fig.write_image(f'{output_name}.png', scale=2)
    fig.write_image(f'{output_name}.svg', scale=2)
    
    # Save the ratio data for reference
    ratio_data = pd.DataFrame({
        'gene': gene_names,
        'protein_coding_tpm': protein_coding_tpm,
        'non_coding_tpm': non_coding_tpm,
        'nc_pc_ratio': ratios,
        'protein_coding_prop': protein_coding_props,
        'non_coding_prop': non_coding_props
    })
    # ratio_data.to_csv(f'{output_name}_ratios.tsv', sep='\t', index=False)

    return True

def main():
    # Load gene list
    with open(args.gene_list) as f:
        genes = [line.strip() for line in f if line.strip()]
    
    print(f"Extracting transcript IDs for {len(genes)} genes...")
    transcript_ids = get_transcript_ids_for_genes(args.gtf, genes)
    
    print("Loading TPM data...")
    tpm_data = load_tpm_data(args.input_list, transcript_ids)

    # tpm_data.to_csv("tpm_data.tsv", index=False, header=True, sep="\t")
    
    # Load all transcript info
    print("Loading transcript annotations...")
    transcript_info = load_transcripts_from_gtf(args.gtf)
    # transcript_info.to_csv("transcript_info.tsv", index=False, header=True, sep="\t")
    
    print(f"Processing {len(genes)} genes...")
    
    # Prepare arguments for parallel processing
    process_args = [(gene, tpm_data, transcript_info, args.topN, args.min_exp, args.outfiles) for gene in genes]
    
    # Process genes in parallel and collect all full_data results
    all_genes_data = {}
    with ProcessPoolExecutor() as executor:
        results = list(executor.map(process_gene, process_args))
        
        # Create a dictionary of gene_name -> full_data
        for i, result in enumerate(results):
            if result is not None:
                gene_name = genes[i]
                all_genes_data[gene_name] = result
    
    successful = len(all_genes_data)
    print(f"Successfully processed {successful} out of {len(genes)} genes")
    
    # Create plots for each gene
    for gene_name, full_data in all_genes_data.items():
        print(f"Creating plots for gene: {gene_name}")
        create_plot(full_data, gene_name, args.outfiles)

if __name__ == "__main__":
    main()
