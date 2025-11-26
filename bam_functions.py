# python module that contains helpful functions for working with bam files

import pysam
import pandas as pd
import numpy as np
import scipy
import re
import seaborn as sns
from Bio import SeqIO
from collections import Counter
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.patches as mpatches
from matplotlib.patches import Rectangle
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

def get_reads(sam, contig, start, stop, strand):
    '''Get reads in a certain region from a BAM file.'''

    data = {
        'contig': [],
        'region_start': [],
        'region_stop': [],
        'region_len': [],
        'strand': [],
        'read_number': [],
        'read_name': [],
        'read_start': [],
        'read_stop': [],
        'read_aln_len': [],
        'read_len_frac': [],
        'read_query_len': [],
        'read_query_overlap': [],
        'cigar_tuples': [],
        'read_seq': [],
        'read_mapq': []
    }

    reads = sam.fetch(contig, start, stop)

    try:
        for i, read in enumerate(reads):
            if (read.is_reverse and strand == '+') or (read.is_forward and strand == '-') or (read.is_secondary):
                pass
            else:
                data['contig'].append(contig)
                data['region_start'].append(start),
                data['region_stop'].append(stop),
                data['region_len'].append(stop-start),
                data['strand'].append(strand),
                data['read_number'].append(i),
                data['read_name'].append(read.query_name),
                data['read_start'].append(read.reference_start),
                data['read_stop'].append(read.reference_end),
                data['read_aln_len'].append(read.reference_length),
                data['read_len_frac'].append(read.reference_length / (stop - start)),
                data['read_query_len'].append(read.infer_query_length()),
                data['read_query_overlap'].append(read.get_overlap(start, stop)),
                data['cigar_tuples'].append(read.cigartuples),
                data['read_seq'].append(read.get_forward_sequence()),
                data['read_mapq'].append(read.mapping_quality)

            
        df = pd.DataFrame(data)
        return df

    except:
        raise

def filter_by_coverage(df, num):
    '''Filter regions in a BAM file by minimum number of reads mapped to that region (not necessarily overlapping).'''
    # print(df[['contig', 'region_start', 'region_stop']].to_string())
    return df[df.groupby(['contig', 'region_start', 'region_stop'])['contig'].transform('size') >= num]

def filter_by_quality(df, score):
    '''Filter reads by mapping quality.'''
    return df.loc[df['read_mapq'] > score]

def filter_by_overlaps(df, min_overlaps):
    '''Filter by number of bases in read aligning the given reference region.'''
    return df.loc[df['read_query_overlap'] >= min_overlaps]

def get_average_error(sam):
    '''Estimate average error of the long-read dataset.'''
    mgs = [] # per-base gap-compressed sequence divergence
    for read in sam.fetch():
        mgs.append(read.get_tag('mg'))

    return 100 - np.mean(mgs)    

def normalize_ends(row):
    '''Get read start and stop relative to the reference region, and get percent of read overlapping the region.'''
    if row['strand'] == '+':
        return ((row['read_start'] - row['region_start']) / (row['region_len']),
                (row['read_stop'] - row['region_start']) / (row['region_len']),
                (row['read_query_overlap'] / row['read_query_len']))
    elif row['strand'] == '-':
        return ((row['region_stop'] - row['read_stop']) / (row['region_len']), 
                (row['region_stop'] - row['read_start']) / (row['region_len']),
                (row['read_query_overlap'] / row['read_query_len']))

def get_mets(seq):
    '''Find all ATGs in sequence regardless of reading frame or codon boxes.'''
    met_inds = [m.start() for m in re.finditer('(?=ATG)', seq)]

    return met_inds

def get_coverage(sam, contig, start, stop):
    '''Get the per-position read coverage of the specified region.'''
    coverage_dict = {key: 0 for key in range(start, stop)}
    cov = np.array(sam.count_coverage(contig, start, stop)) # per base / per position coverage

    for i in range(cov.shape[1]):
        coverage_dict[i + start] = np.sum([arr[i] for arr in cov], axis=0)

    return coverage_dict

def get_splice_data(row):
    '''Get the splice junctions and introns for BAM read data.'''
    start, stop, cigar, strand = row[['read_start', 'read_stop', 'cigar_tuples', 'strand']]
    trans_struct = None
    splice_juncs = None
    introns = []
    pos = None

    pos = start
    splice = []

    for tup in cigar:
        if tup[0] not in [3, 1, 4]:
            pos += tup[1]
        elif tup[0] == 3:
            introns.append((pos, pos + tup[1]))

            splice.append(pos)
            pos += tup[1]
            splice.append(pos-1)

    splice_juncs = tuple(splice)

    trans_struct = tuple([start, *splice, stop])

    return introns, trans_struct, splice_juncs

def structure_abundance(df, normalize=False):
    '''Return the frequency (fractional or count) of the transcript structure (start, splices, stop).'''
    value_count = pd.DataFrame()

    if normalize:
        value_count = df.groupby(['contig', 'region_start', 'region_stop', 'strand', 'internal_splice_juncs']).internal_splice_juncs.agg(lambda x: x.count().div(len(x)))
    else:
        value_count = df.groupby(['contig', 'region_start', 'region_stop', 'strand', 'internal_splice_juncs']).internal_splice_juncs.agg(count=(lambda x: x.count()))

    return value_count

def sashimi_junc(ax, introns):
    '''Render sashimi plot of read coverage / splice junctions.'''
    Path = mpath.Path
    introns = Counter([i for sublist in introns for i in sublist])
    for intron, count in introns.items():
        pt1 = (intron[0], count)
        midpoint = ((intron[0] + intron[1])/2, 2*count)
        pt2 = (intron[1], count)

        if count > 1:
            pp1 = mpatches.PathPatch(
                Path([pt1, midpoint, pt2], [Path.MOVETO, Path.CURVE3, Path.CURVE3]),
                fc="none", transform=ax.transData)
            ax.text((intron[0] + intron[1])/2, 1.5*count, str(count), 
                fontsize=5, va='center', ha='center', backgroundcolor='w',
                bbox={'facecolor':'white', 'edgecolor': 'white', 'pad':0.2}, clip_on=True)

        else:
            midpoint = ((intron[0] + intron[1])/2, max(introns.values()))
            Path = mpath.Path
            pp1 = mpatches.PathPatch(
                Path([pt1, midpoint, pt2], [Path.MOVETO, Path.CURVE3, Path.CURVE3]),
                fc="none", transform=ax.transData, ls=':')

        ax.add_patch(pp1)

def get_reliability_score(df, thresh=0.2):
    '''Calculate reliability score to classify long-read datasets as suitable or not for downstream analysis.'''

    bounded_fl_annotated = df[df['norm_read_starts'].between(-0.1, 0.1) & \
        df['norm_read_stops'].between(0.9, 1.0)]
    bounded_fl_annotated_perc = len(bounded_fl_annotated) / len(df)

    svals, bs = np.histogram(df['norm_read_starts'], bins=np.arange(-0.25, 1.25, 0.1), density=True)
    stvals, bs = np.histogram(df['norm_read_stops'], bins=np.arange(-0.25, 1.25, 0.1), density=True)

    unexpected_peak = False

    unannotated_starts = sum([x for i, x in enumerate(svals) if x!=2]) 
    unannotated_stops = sum([x for i, x in enumerate(stvals) if x!=12])

    if unannotated_starts > thresh or unannotated_stops > thresh:
        unexpected_peak = True


    return bounded_fl_annotated_perc, unexpected_peak

def plot_reliability(df, key, axes):
    '''Plot comprehensive reliability plots showing start and end site frequencies, and query coverage.'''
    data = df[key]
    flat_overlaps = df['norm_overlaps'].to_list()

    ### UNBOUNDED START / END SITE HISTOGRAM ###

    flat_ends = data.to_list()

    stats, edges, nums = scipy.stats.binned_statistic(x=flat_ends, values=flat_overlaps, statistic='mean', bins=40)

    n, bins, patches = axes[0].hist(flat_ends, log=True, bins=40, edgecolor='black', linewidth=1)
    cm = plt.cm.binary

    for i, p in enumerate(patches):
        plt.setp(p, 'facecolor', cm(stats[i]))

    norm = mpl.colors.Normalize(vmin=0, vmax=1)
    sm = plt.cm.ScalarMappable(cmap=cm, norm=norm)
    sm.set_array([])
    clb = plt.colorbar(sm, ax=axes[0], label='Average Fraction of\nRead Overlapping Gene')

    clb.ax.set_label('Average Fraction of Read Overlapping Gene')

    ### BOUNDED START / END SITE COUNTS ###

    flat_lens = df['read_len_frac'].to_list()

    vals, bs = np.histogram(flat_ends, bins=np.arange(-0.225, 1.225, 0.05))
    stats, edges, nums = scipy.stats.binned_statistic(x=flat_ends, values=flat_lens, statistic=lambda x: len([i for i in x if i > 0.9 and i < 1.1]), bins=np.arange(-0.225, 1.225, 0.05))

    axes[1].hist(bs[:-1], bs, weights=stats, linewidth=1, facecolor='lightgrey', edgecolor='black', label='"full length" reads')
    # print(bs)
    axes[1].set_xticks([x + 0.025 for x in bs[:-1:2]])
    axes[1].set_xticklabels(['{:.1f}'.format(x + 0.025) for x in bs[:-1:2]])

    axes[1].plot([x + 0.025 for x in bs[:-1]], vals, color='red', label='all reads')
    axes[1].legend()

    ### START / END SITE POSITION HEAT MAP ###

    hist_array = []

    for gene, gene_vals in df.sort_values(by=['region_len'], ascending=True).groupby(['contig', 'region_start', 'region_stop'], sort=False):
        gene_end = gene_vals[key].to_list()
        if key == 'norm_read_starts':
            vals, bs = np.histogram(gene_end, bins=np.arange(-1.05, 1.05, 0.10))
            hist_array.append(vals)
        elif key == 'norm_read_stops':
            vals, bs = np.histogram(gene_end, bins=np.arange(-0.05, 2.05, 0.10))
            hist_array.append(vals)

    hist_array = np.float_(hist_array)

    tots = np.sum(hist_array, axis=1, keepdims=True)
    norm_array = np.divide(hist_array, tots, out=np.zeros_like(hist_array), where=tots!=0)
    sns.heatmap(norm_array, ax=axes[2], xticklabels='auto')
    axes[2].set_xticklabels(['{:.2f}'.format(x + 0.05) for x in bs[0:-1:2]], rotation=45)
    axes[2].set_yticks(np.arange(0, len(norm_array), 100))
    axes[2].set_yticklabels(np.arange(0, len(norm_array), 100))


def find_in_ref(ref, contig, start, stop, seq, surround=0):
    '''Get occurences of a given sequence in a larger fasta file (intended for reference of BAM file).'''
    chrm = None
    for record in SeqIO.parse(ref, "fasta"):
        if record.id == contig:
            chrm = record

    region = str(chrm.seq[start: stop]).upper()
    
    pattern = re.compile(seq)
    
    matches = [m for m in re.finditer(seq, region)]
    
    seqs = []

    for match in matches:
        idx, end = match.start(), match.end()
        length = end - idx
        rresult = region[idx - surround: idx + length + surround]
        cresult = chrm.seq.upper()[start + idx - surround: start + idx + length + surround]
        
        seqs.append(cresult)
        
    return [x.start() + start for x in matches], seqs

def plot_ends_hist(sam, contig, start, stop, strand):
    '''Deprecated function to plot coverage of starts and ends of reads.'''
    reads, stops, starts, lens = get_reads(sam, contig, start, stop, strand)

    sns.set_style('darkgrid')
    fig, ax = plt.subplots(1, 1, figsize=(12, 2))

    ax.hist(starts, bins=range(start, stop + 1, 3), log=True, color='green', edgecolor='green')
    ax.hist(stops, bins=range(start, stop + 1, 3), log=True, color='red', edgecolor='red')

    ax.set_xlim(start - 1000, stop)
    ax.set_xticks(np.linspace(start, stop - 2000, 6))

    plt.show()
