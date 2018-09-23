# vim: fdm=indent
'''
author:     Fabio Zanini
date:       27/05/18
content:    Supplementary Fig 13 shows the level of cross-talk in virus reads.
'''
# Modules
import os
import sys
import argparse
import numpy as np
import pandas as pd
import igraph as ig
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

os.environ['SINGLET_CONFIG_FILENAME'] = 'singlet.yml'
sys.path.append('/home/fabio/university/postdoc/singlet')
from singlet.dataset import Dataset




# Script
if __name__ == '__main__':

    pa = argparse.ArgumentParser(description='Find viral reads in healthy samples.')
    args = pa.parse_args()

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

    print('Set normalized dengue reads')
    n = ds.samplesheet['numberDengueReads'].astype(float)
    cov = ds.samplesheet['coverage'].astype(float)
    ds.samplesheet['virus_reads_per_million'] = 1e6 * n / (n + cov)
    ds.samplesheet['log_virus_reads_per_million'] = (np.log10(0.1 + ds.samplesheet['virus_reads_per_million'])).fillna(0)

    print('Annotate dengue severity')
    ds.samplesheet['dengue_severity'] = 0
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
    # Add back CD123 which is on X
    genes_good['IL3RA'] = True
    # Discard HLA
    genes_good &= ~ds.featurenames.str.startswith('HLA')
    # Discard IGHV/IGKV/IGLV
    genes_good &= ~ds.featurenames.str.startswith('IGHV')
    genes_good &= ~ds.featurenames.str.startswith('IGKV')
    genes_good &= ~ds.featurenames.str.startswith('IGLV')
    # Discard TRBV
    genes_good &= ~ds.featurenames.str.startswith('TRBV')
    # Discard mitochondrial RNA
    genes_good &= ~ds.featurenames.str.startswith('MTR')
    # Discard ribosomal
    genes_good &= ~ds.featurenames.str.startswith('RPS')
    genes_good &= ~ds.featurenames.str.startswith('RPL')
    genes_good &= ~ds.featurenames.str.startswith('RP3')
    genes_good &= ~ds.featurenames.str.startswith('RP5')
    genes_good &= ~ds.featurenames.str.startswith('RP11')
    ds.counts = ds.counts.loc[genes_good]

    # Check background DENV reads in healthy controls
    ns = {}
    for key, datum in ds.samplesheet.loc[:, ['dengue_severity', 'numberDengueReads']].groupby('dengue_severity'):
        ns[key] = datum['numberDengueReads'].values
    labels = {0: 'healthy', 1: 'dengue', 2: 'severe dengue'}
    colors = sns.color_palette(n_colors=3)
    fig, ax = plt.subplots(1, 1, figsize=(4.25, 3.2))
    for key, datum in ns.items():
        x = 0.1 + np.sort(datum)
        y = 1.0 - np.linspace(0, 1, len(x)) + 1e-4
        ax.plot(x, y, label=labels[key], color=colors[key])
        ax.axhline(1.0 / len(x), color=colors[key], ls='--')
    ax.plot([30, 30], [0.9e-4, 1.5], lw=1.5, color='grey', ls='--')
    ax.legend(loc='upper right', fontsize=8, ncol=3)
    ax.set_xlim(0.09, 60)
    ax.set_xscale('log')
    ax.set_ylim(0.9e-4, 1.2)
    ax.set_yscale('log')
    ax.set_xlabel('# of DENV reads')
    ax.set_ylabel('Fraction of cells with >= x DENV reads')
    plt.tight_layout()
    fig.savefig('../../data/supplementary_fig7.png', dpi=600)
    fig.savefig('../../data/supplementary_fig7.svg')

    plt.ion()
    plt.show()
