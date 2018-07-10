# vim: fdm=indent
'''
author:     Fabio Zanini
date:       13/06/18
content:    Look at gene expression of bystander B cells in patients with
            infected cells.
'''
# Modules
import os
import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

os.environ['SINGLET_CONFIG_FILENAME'] = 'singlet.yml'
sys.path.append('/home/fabio/university/postdoc/singlet')
from singlet.dataset import Dataset, CountsTable, FeatureSheet




# Script
if __name__ == '__main__':

    print('Load data (L1 normalized)')
    ds = Dataset(
            samplesheet='dengue_patients',
            counts_table='dengue_patients',
            featuresheet='dengue_patients',
            )

    print('Filter low quality cells')
    ds.samplesheet['coverage'] = ds.samplesheet['coverage'].astype(int)
    ds.query_samples_by_metadata('cellType not in ("unknown", "PBMC", "low-quality", "ambiguous")', inplace=True)

    print('Translate gene names')
    ds.featuresheet['EnsemblID'] = ds.featurenames
    ds.rename(axis='features', column='GeneName', inplace=True)

    print('Annotate dengue severity')
    ds.samplesheet['dengue_severity'] = 0
    # FIXME FIXME FIXME: 10017014 is SD, 10017015 is DF!!!
    ds.samplesheet.loc[
            ds.samplesheet['experiment'].isin(['10017013', '10017015']),
            'dengue_severity'] = 1
    ds.samplesheet.loc[
            ds.samplesheet['experiment'].isin(['10017014', '10017016', '10017021', '10017022']),
            'dengue_severity'] = 2
    ds.samplesheet['severe_dengue'] = ds.samplesheet['dengue_severity'] == 2

    print('Ignore X/Y-linked, HLA types and TCR and BCR variable regions')
    # Keep only autosomal genes
    genes_good = ds.featuresheet['Chromosome/scaffold name'].isin([str(i+1) for i in range(23)])
    # Discard HLA
    #genes_good &= ~ds.featurenames.str.startswith('HLA')
    genes_good &= ~ds.featurenames.str.startswith('HLA-A')
    genes_good &= ~ds.featurenames.str.startswith('HLA-B')
    genes_good &= ~ds.featurenames.str.startswith('HLA-C')
    genes_good &= ~ds.featurenames.str.startswith('HLA-E')
    genes_good &= ~ds.featurenames.str.startswith('HLA-F')
    genes_good &= ~ds.featurenames.str.startswith('HLA-H')
    genes_good &= ~ds.featurenames.str.startswith('HLA-L')
    # Discard IGHV/IGKV/IGLV
    genes_good &= ~ds.featurenames.str.startswith('IGHV')
    genes_good &= ~ds.featurenames.str.startswith('IGKV')
    genes_good &= ~ds.featurenames.str.startswith('IGLV')
    # Discard TRBV
    genes_good &= ~ds.featurenames.str.startswith('TRAV')
    genes_good &= ~ds.featurenames.str.startswith('TRAJ')
    genes_good &= ~ds.featurenames.str.startswith('TRBV')
    genes_good &= ~ds.featurenames.str.startswith('TRBJ')
    # Discard mitochondrial RNA
    genes_good &= ~ds.featurenames.str.startswith('MTR')
    # Discard ribosomal
    genes_good &= ~ds.featurenames.str.startswith('RPS')
    genes_good &= ~ds.featurenames.str.startswith('RPL')
    genes_good &= ~ds.featurenames.str.startswith('RP3')
    genes_good &= ~ds.featurenames.str.startswith('RP5')
    genes_good &= ~ds.featurenames.str.startswith('RP11')
    ds.counts = ds.counts.loc[genes_good]

    print('Set normalized dengue reads')
    n = ds.samplesheet['numberDengueReads'].astype(float)
    cov = ds.samplesheet['coverage'].astype(float)
    ds.samplesheet['virus_reads_per_million'] = 1e6 * n / (n + cov)
    ds.samplesheet['log_virus_reads_per_million'] = (np.log10(0.1 + ds.samplesheet['virus_reads_per_million'])).fillna(0)

    print('Focus on B cells')
    dsB = ds.query_samples_by_metadata('cellType == "B cell"')

    print('Restrict to patients 10017016 and 10017022, the ones with infected cells')
    dss = dsB.query_samples_by_metadata('experiment in ("10017016", "10017022")')

    print('Divide by infected, bystanders, and controls')
    dsi = dss.query_samples_by_metadata('virus_reads_per_million >= 30')
    dsb = dss.query_samples_by_metadata('virus_reads_per_million == 0')
    # NOTE: think whether this control is appropriate
    dsc = dsB.query_samples_by_metadata('dengue_severity == 0')

    print('Differential expression between bystanders and controls')
    comp = dsb.compare(dsc)

    print('Plot distributions in bystanders vs controls')
    # NOTE: these genes come out of the smallest P-value for 'comp' above
    #genes_diff = ['IFI6', 'IFI44L', 'IFIT3', 'ISG15', 'MX1']

    # FIG 4G
    genes_diff = ['IFI6', 'IFI44L', 'IFIT3']
    fig, axs = plt.subplots(
            1, len(genes_diff), sharex=True, sharey=True,
            figsize=(3.6, 2.3))
    axs = axs.ravel()
    colors = [sns.color_palette(n_colors=7)[i] for i in [2, 0]]
    for gname, ax in zip(genes_diff, axs):
        dfi = comp.loc[gname]
        datab = np.log10(0.1 + dsb.counts.loc[[gname]].T)
        datab['is_bystander'] = 'yes'
        datac = np.log10(0.1 + dsc.counts.loc[[gname]].T)
        datac['is_bystander'] = 'no'
        data = pd.concat([datab, datac], axis=0)
        with sns.axes_style('whitegrid'):
            sns.violinplot(
                    data=data,
                    x='is_bystander',
                    y=gname,
                    inner='quartile',
                    ax=ax,
                    order=['no', 'yes'],
                    zorder=10,
                    scale='width',
                    palette=colors,
                    )
            plt.setp(ax.collections, alpha=.6)
            sns.swarmplot(
                    data=data,
                    x='is_bystander',
                    y=gname,
                    ax=ax,
                    order=['no', 'yes'],
                    zorder=10,
                    s=1.5,
                    alpha=1.0,
                    palette=colors,
                    )
        ax.set_title(gname)
        ax.set_ylabel('')
        ax.set_xlabel('')
        ax.text(0.06, 0.96, 'P = {:.0e}'.format(dfi['P-value']),
                ha='left',
                va='top',
                transform=ax.transAxes,
                bbox={'edgecolor': 'grey', 'facecolor': 'white', 'alpha': 0.9, 'lw': 1, 'pad': 3},
                )
        for y in [0, 1, 2, 3, 4, 5]:
            ax.plot([-1, 2], [y] * 2, color='grey', alpha=0.5, zorder=0, lw=1)

    ax.set_ylim(-1.1, 6.3)
    ax.set_yticks([-1, 0, 1, 2, 3, 4, 5])
    ax.set_yticklabels(['$0$', '$1$', '$10$', '$10^2$', '$10^3$', '$10^4$', '$10^5$'])
    fig.text(0.55, 0.08,
             'is the B cell a bystander?',
             ha='center',
             )
    fig.text(0.02, 0.65, 'counts per million', rotation=90, ha='center')
    plt.tight_layout(rect=(0.03, 0.03, 1, 1), w_pad=0)
    t = fig.text(0.01, 0.99, 'F', ha='left', va='top', fontsize=16)
    fig.savefig('../../figures/fig4G.png')
    fig.savefig('../../figures/fig4G.svg')

    plt.ion()
    plt.show()
