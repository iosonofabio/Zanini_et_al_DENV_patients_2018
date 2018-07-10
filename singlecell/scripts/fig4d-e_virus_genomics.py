# vim: fdm=indent
'''
author:     Fabio Zanini
date:       04/06/18
content:    Try and explore a bit the virus genomics.
'''
import os
import sys
import argparse
import numpy as np
import pandas as pd
import pysam
from Bio import SeqIO
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

from singlecell.filenames import experiments_foldername_small as experiments_foldername
from singlecell.util import pair_generator


def trim_short_cigars(read, min_match_length=25):
    cigar = read.cigar

    # Check there is at least a match
    if not any(ct == 0 for (ct, cl) in cigar):
        raise ValueError('read has no matches')

    # Find longest match
    lmax = max(cl for (ct, cl) in cigar if ct == 0)
    if lmax < min_match_length:
        raise ValueError('read has short matches only')

    # Find matches that are long enough
    is_good_match = np.array([(ct == 0) and (cl >= min_match_length) for (ct, cl) in cigar])
    tmp = np.nonzero(is_good_match)[0]
    first_good, last_good = tmp[0], tmp[-1]
    if first_good != last_good:
        is_good_match[first_good: last_good] = True

    # Trim the read
    posr = read.reference_start
    seq = []
    qual = []
    posq = 0
    is_first_match = True
    for ic, (ct, cl) in enumerate(cigar):
        if is_good_match[ic]:
            if ct == 0:
                seq.append(read.seq[posq: posq + cl])
                qual.append(read.query_qualities[posq: posq + cl])
                posq += cl
                is_first_match = False
            elif ct == 1:
                seq.append(read.seq[posq: posq + cl])
                qual.append(read.query_qualities[posq: posq + cl])
                posq += cl
            elif ct == 2:
                continue
            else:
                ValueError('CIGAR not in (0, 1, 2)')
        else:
            if ct in (0, 1):
                posq += cl
            elif ct == 2:
                pass
            else:
                ValueError('CIGAR not in (0, 1, 2)')
            if (ct in (0, 2)) and is_first_match:
                posr += cl
    seq = ''.join(seq)
    qual = read.query_qualities.__class__('B', np.concatenate(qual))

    # Set the read
    read.pos = posr
    read.seq = seq
    read.query_qualities = qual
    read.cigar = [t for igc, t in zip(is_good_match, cigar) if igc]
    return read


def get_allele_counts(fn_bam, ref):
    alphal = ['A', 'C', 'G', 'T', '-', 'N']
    alpha = np.array(alphal)
    ac = np.zeros((len(alpha), len(ref)), int)
    with pysam.AlignmentFile(fn_bam, 'rb') as bamfile:
        for ir, (read1, read2) in enumerate(pair_generator(bamfile)):
            seq = np.array(['X' for i in ref])
            try:
                (read1, read2) = (trim_short_cigars(read1), trim_short_cigars(read2))
            except ValueError:
                continue

            for read in (read1, read2):
                for (posq, posr) in read.get_aligned_pairs():
                    if posr is None:
                        continue
                    # If the position is unset, set it
                    if seq[posr] == 'X':
                        if posq is None:
                            seq[posr] = '-'
                        else:
                            seq[posr] = read.seq[posq]

                    else:
                        if (posq is None) and (seq[posr] == '-'):
                            continue
                        # If the two reads do not agree, unset it
                        if (posq is None) or (read.seq[posq] != seq[posr]):
                            seq[posr] = 'X'

            # Fill the allele count matrix
            for inu, nuc in enumerate(alpha):
                ac[inu, seq == nuc] += 1

    return {
        'allele_counts': ac,
        'alpha': alpha,
        }



if __name__ == '__main__':

    pa = argparse.ArgumentParser(description='Explore virus genomics data')
    pa.add_argument('--experiment', choices=['10017016', '10017022'],
                    default='10017022',
                    help='What experiment to analyze')
    args = pa.parse_args()

    print('Get hybrid reference')
    fn_ref = '../../data/dengue_reference_hybrid_{:}.fasta'.format(args.experiment)
    ref = SeqIO.read(fn_ref, 'fasta')

    print('Get allele counts')
    fn_bam = '../../data/dengue_remapped_{:}.bam'.format(args.experiment)
    tmp = get_allele_counts(fn_bam, ref)
    ac = tmp['allele_counts']
    alpha = tmp['alpha']

    print('Get minor allele frequencies')
    af_min = np.ma.masked_all(len(ref), float)
    for pos, a in enumerate(ac.T):
        if a.sum() >= 10:
            i_min = np.argsort(a)[-2]
            af_min[pos] = 1.0 * a[i_min] / a.sum()

    # FIG 4G
    print('Plot coverage and allele frequencies')
    x = np.arange(len(ref)) + 1
    fig, ax = plt.subplots(1, 1, figsize=(6.27, 2.1))
    ax.plot(x, 0.1 + ac.sum(axis=0), lw=2, color='darkred', label='Coverage', zorder=5, alpha=0.7)
    ax.set_xlim(0, len(ref) + 1)
    ax.set_ylim(ymin=0.09)
    ax.set_yscale('log')
    ax.set_xlabel('Position in DENV genome')
    ax2 = ax.twinx()
    ax2.plot(x, af_min, lw=2, color='steelblue', label='Minor allele frequency', zorder=4, alpha=0.6)
    ax2.set_xlim(0, len(ref) + 1)
    ax2.set_ylim(0.009, 1.1)
    ax2.set_yscale('log')

    ax.set_ylabel('Coverage (red)')
    ax2.set_ylabel('MAF (blue)')
    fig.text(0.01, 0.98, 'G', ha='left', va='top', fontsize=14)
    fig.tight_layout(rect=(0.01, 0, 0.99, 1))
    fig.savefig('../../figures/fig4G.png')
    fig.savefig('../../figures/fig4G.svg')

    ## Try the same with 3 plots
    ## FIG 4G
    #x = np.arange(len(ref)) + 1
    ##start = 5280
    ##end = 5640
    #start = 7100
    #end = 7500
    #start2 = 10200
    #end2 = 10600

    #from matplotlib import gridspec
    #fig = plt.figure(figsize=(8.27, 4.5))
    #gs = gridspec.GridSpec(2, 2)
    #axs = [fig.add_subplot(gs[0, :])]
    #axs.append(fig.add_subplot(gs[2], sharey=axs[0]))
    #axs.append(fig.add_subplot(gs[3], sharey=axs[0]))

    #ax = axs[0]
    #ax.plot(x, 0.1 + ac.sum(axis=0), lw=2, color='darkred', label='Coverage', zorder=5, alpha=0.7)
    #ax.set_xlim(0, len(ref) + 1)
    #ax.set_ylim(ymin=0.09)
    #ax.set_yscale('log')
    #ax2 = ax.twinx()
    #ax2.plot(x, af_min, lw=2, color='steelblue', label='Minor allele frequency', zorder=4, alpha=0.6)
    #ax2.set_xlim(0, len(ref) + 1)
    #ax2.set_ylim(0.009, 1.1)
    #ax2.set_yscale('log')
    #ax2.plot([start] * 2, [1e-3, 1.1], lw=1.5, color='k', zorder=8)
    #ax2.plot([end] * 2, [1e-3, 1.1], lw=1.5, color='k', zorder=8)
    #ax2.plot([start2] * 2, [1e-3, 1.1], lw=1.5, color='k', zorder=8)
    #ax2.plot([end2] * 2, [1e-3, 1.1], lw=1.5, color='k', zorder=8)

    #ax = axs[1]
    #ax.plot(x[start: end], 0.1 + ac.sum(axis=0)[start: end], lw=2, color='darkred', label='Coverage', zorder=5, alpha=0.7)
    #ax.set_xlim(start, end + 1)
    #ax.set_ylim(ymin=0.09)
    #ax.set_yscale('log')
    #ax2 = ax.twinx()
    #ax2.plot(x[start: end], af_min[start: end], lw=2, color='steelblue', label='Minor allele frequency', zorder=4, alpha=0.6)
    #ax.set_xlim(start, end + 1)
    #ax2.set_ylim(0.009, 1.1)
    #ax2.set_yscale('log')

    #plt.setp(ax2.get_yticklabels(), visible=False)

    #ax = axs[2]
    #ax.plot(x[start2: end2], 0.1 + ac.sum(axis=0)[start2: end2], lw=2, color='darkred', label='Coverage', zorder=5, alpha=0.7)
    #ax.set_xlim(start, end + 1)
    #ax.set_ylim(ymin=0.09)
    #ax.set_yscale('log')
    #ax2 = ax.twinx()
    #ax2.plot(x[start2: end2], af_min[start2: end2], lw=2, color='steelblue', label='Minor allele frequency', zorder=4, alpha=0.6)
    #ax.set_xlim(start2, end2 + 1)
    #ax2.set_ylim(0.009, 1.1)
    #ax2.set_yscale('log')
    #plt.setp(axs[2].get_yticklabels(), visible=False)

    #fig.text(0.03, 0.52, 'Coverage (red)', ha='center', va='center', rotation=90)
    #fig.text(0.965, 0.52, 'MAF (blue)', ha='center', va='center', rotation=90)
    #fig.text(0.01, 0.98, 'G', ha='left', va='top', fontsize=14)
    #fig.tight_layout(rect=(0.03, 0, 0.97, 1), w_pad=0.1)
    ##fig.savefig('../../figures/fig4G.png')
    ##fig.savefig('../../figures/fig4G.svg')

    print('Load cross-sectional allele frequencies')
    from singlecell.filenames import support_foldername
    cs_fn = support_foldername+'dengue/alignments/derived/allele_frequencies_dengue3.npz'
    tmp = np.load(cs_fn)
    alpha = tmp['alpha']
    afs = tmp['allele_frequencies']
    cons = ''.join(alpha[np.argmax(afs, axis=0)])
    ind_nogaps = np.array([c != '-' for c in cons])
    indi_nogaps = np.nonzero(ind_nogaps)[0]
    cons_nogaps = ''.join(np.array(list(cons))[ind_nogaps])
    afs_nogaps = afs[:, ind_nogaps]

    print('Align hybrid reference and cross sectional allele frequencies')
    from seqanpy import align_global, align_overlap
    score, ali1, ali2 = align_overlap(ref, cons_nogaps, score_gapopen=-20)
    ali_fn = support_foldername+'sequencing/experiments/{:}/ali_to_consensus_nongap.fasta'.format(
            args.experiment)
    #with open(ali_fn, 'wt') as f:
    #    f.write('> hybridReference_{:}\n'.format(args.experiment))
    #    f.write(ali1+'\n')
    #    f.write('> serotypeConsNoGaps\n')
    #    f.write(ali2+'\n')

    # No gaps in the largest region 89: 10280
    # NOTE: annotate the UTRs
    if args.experiment == '10017022':
        pos_cs = np.arange(89, 10280)
    else:
        raise ValueError('Range not yet defined for this experiment')
    afs_cs_ali = afs_nogaps[:, pos_cs]
    start_hyb = len(ali2) - len(ali2.lstrip('-'))
    pos_hyb = pos_cs - start_hyb
    afs_pat = 1.0 * ac / ac.sum(axis=0)
    afs_pat_ali = afs_pat[:, pos_hyb]

    print('Calculate entropies')
    def get_site_entropy(nu, alpha=['A', 'C', 'G', 'T', '-', 'N'], gaps_max=0.2):
        '''Get site entropy from allele frequencies'''
        s = np.ma.array(-(nu * np.log2(nu + 1e-10)).sum(axis=0),
                        shrink=False)
        s[s < 0] = 0
        # Exclude sites of many gaps
        s.mask = np.zeros(len(s), bool)
        s.mask[nu[-2] > gaps_max] = True
        return s
    ent_cs_ali = get_site_entropy(afs_cs_ali)
    ent_pat_ali = get_site_entropy(afs_pat_ali)

    print('Restrict to decent intrapatient coverage')
    cov_min = 200
    cov_ali = ac.sum(axis=0)[pos_hyb]
    ent_cs_cov = ent_cs_ali[cov_ali >= cov_min]
    ent_pat_cov = ent_pat_ali[cov_ali >= cov_min]

    # FIG 4H
    print('Scatter entropies')
    from scipy.stats import pearsonr, spearmanr
    r = pearsonr(ent_cs_cov, ent_pat_cov)[0]

    fig, ax = plt.subplots(figsize=(2, 2.1))
    ax.scatter(
            ent_cs_cov, ent_pat_cov, s=20, alpha=0.5,
            label='$r={:.2f}$'.format(r),
            )
    ax.set_xlabel('Serotype entropy [bits]', fontsize=9)
    ax.set_ylabel('Patient entropy [bits]', fontsize=9)
    ax.set_xlim(xmin=1e-2)
    ax.set_ylim(ymin=1e-2)
    ax.set_xscale('log')
    ax.set_yscale('log')
    #ax.legend(loc='lower right')
    ax.text(0.99, 0.01, '$r={:.2f}$'.format(r), ha='right', va='bottom',
            transform=ax.transAxes)
    ax.grid(True)
    fig.text(0.01, 0.98, 'H', ha='left', va='top', fontsize=16)
    plt.tight_layout(rect=(0, 0, 1, 0.9))
    fig.savefig('../../figures/fig4H.png')
    fig.savefig('../../figures/fig4H.svg')

    plt.ion()
    plt.show()
