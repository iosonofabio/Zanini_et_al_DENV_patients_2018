# vim: fdm=indent
'''
author:     Fabio Zanini
date:       17/03/18
content:    Quality controls by patient and cell type.
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
    ds.samplesheet['infected'] = ds.samplesheet['virus_reads_per_million'] >= 30
    ds.samplesheet['numberGenes'] = (ds.counts >= 1).sum(axis=0)

    # Plot housekeeping by sample
    pats = ['3-013-1', '3-027-1', '3-018-1', '3-006-1',
            '1-008-1', '1-020-1',
            '1-013-1', '1-026-1', '1-010-1', '1-036-1']
    genes_house = ['ACTB', 'ACTG1', 'B2M', 'PTPRC', 'UBC', 'UBB', 'PSMB4', 'GAPDH']
    fig, axs = plt.subplots(4, 2, figsize=(5, 7), sharex=True, sharey=True)
    axs = axs.ravel()
    for gene, ax in zip(genes_house, axs):
        df = ds.samplesheet[['patientSample']].copy()
        df[gene] = np.log10(0.1 + ds.counts.loc[gene])

        sns.violinplot(
                x='patientSample',
                y=gene,
                data=df,
                inner='quartile',
                ax=ax,
                order=pats,
                scale='width',
                palette=sns.color_palette(n_colors=len(pats)),
                )
        plt.setp(ax.collections, alpha=.8)
        ax.set_xlabel('')
        for tk in ax.get_xticklabels():
            tk.set_rotation(90)

    plt.tight_layout(w_pad=0.5, h_pad=0.1)
    fig.savefig('../../figures/qc_housekeeping_bysample.png', dpi=600)
    fig.savefig('../../figures/qc_housekeeping_bysample.svg')

    # Plot housekeeping by cell type
    cell_types = ['T cell', 'NK cell', 'NKT cell', 'B cell', 'monocyte', 'cDC', 'pDC']
    fig, axs = plt.subplots(4, 2, figsize=(5, 7), sharex=True, sharey=True)
    axs = axs.ravel()
    for gene, ax in zip(genes_house, axs):
        df = ds.samplesheet[['patientSample', 'cellType', 'infected']].copy()
        df[gene] = np.log10(0.1 + ds.counts.loc[gene])

        sns.violinplot(
                x='cellType',
                y=gene,
                data=df,
                inner='quartile',
                ax=ax,
                order=cell_types,
                scale='width',
                palette=sns.color_palette(n_colors=len(cell_types)),
                )
        plt.setp(ax.collections, alpha=.8)
        ax.set_xlabel('')
        for tk in ax.get_xticklabels():
            tk.set_rotation(90)

    plt.tight_layout(w_pad=0.5, h_pad=0.1)
    fig.savefig('../../figures/qc_housekeeping_bycelltype.png', dpi=600)
    fig.savefig('../../figures/qc_housekeeping_bycelltype.svg')

    # Plot housekeeping by association with DENV
    dsip = ds.query_samples_by_metadata('(patient in ("1-026", "1-036")) & (cellType == "B cell")')
    fig, axs = plt.subplots(4, 2, figsize=(2.5, 7), sharex=True, sharey=True)
    axs = axs.ravel()
    for gene, ax in zip(genes_house, axs):
        df = dsip.samplesheet[['patientSample', 'cellType', 'infected']].copy()
        df.loc[df['infected'] == True, 'infected'] = 'with DENV RNA'
        df.loc[df['infected'] == False, 'infected'] = 'no DENV RNA'
        df[gene] = np.log10(0.1 + dsip.counts.loc[gene])

        sns.violinplot(
                x='infected',
                y=gene,
                data=df,
                inner='quartile',
                ax=ax,
                order=['no DENV RNA', 'with DENV RNA'],
                scale='width',
                palette=sns.color_palette(n_colors=2),
                )
        plt.setp(ax.collections, alpha=.8)
        ax.set_xlabel('')
        for tk in ax.get_xticklabels():
            tk.set_rotation(90)

    plt.tight_layout(w_pad=0.5, h_pad=0.1)
    fig.savefig('../../figures/qc_housekeeping_byDENV_Bcells2636.png', dpi=600)
    fig.savefig('../../figures/qc_housekeeping_byDENV_Bcells2636.svg')

    # Plot number of genes detected by sample
    pats = ['3-013-1', '3-027-1', '3-018-1', '3-006-1',
            '1-008-1', '1-020-1',
            '1-013-1', '1-026-1', '1-010-1', '1-036-1']
    fig, ax = plt.subplots(1, 1, figsize=(2.5, 2.5))
    df = ds.samplesheet[['patientSample', 'numberGenes']].copy()
    sns.violinplot(
            x='patientSample',
            y='numberGenes',
            data=df,
            inner='quartile',
            ax=ax,
            order=pats,
            scale='width',
            palette=sns.color_palette(n_colors=len(pats)),
            )
    plt.setp(ax.collections, alpha=.8)
    ax.set_xlabel('')
    ax.set_ylabel('number of genes >= 1 cpm')
    for tk in ax.get_xticklabels():
        tk.set_rotation(90)
    ax.set_ylim(0, 4000)

    plt.tight_layout(w_pad=0.5, h_pad=0.1)
    fig.savefig('../../figures/qc_number_genes_bysample.png', dpi=600)
    fig.savefig('../../figures/qc_number_genes_bysample.svg')

    # Plot number of genes detected by cell type
    cell_types = ['T cell', 'NK cell', 'NKT cell', 'B cell', 'monocyte', 'cDC', 'pDC']
    fig, ax = plt.subplots(1, 1, figsize=(2.5, 2.5))
    df = ds.samplesheet[['cellType', 'numberGenes']].copy()
    sns.violinplot(
            x='cellType',
            y='numberGenes',
            data=df,
            inner='quartile',
            ax=ax,
            order=cell_types,
            scale='width',
            palette=sns.color_palette(n_colors=len(cell_types)),
            )
    plt.setp(ax.collections, alpha=.8)
    ax.set_xlabel('')
    ax.set_ylabel('number of genes >= 1 cpm')
    for tk in ax.get_xticklabels():
        tk.set_rotation(90)
    ax.set_ylim(0, 4000)

    plt.tight_layout(w_pad=0.5, h_pad=0.1)
    fig.savefig('../../figures/qc_number_genes_bycelltype.png', dpi=600)
    fig.savefig('../../figures/qc_number_genes_bycelltype.svg')

    # Plot number of genes detected by DENV RNA
    dsip = ds.query_samples_by_metadata('(patient in ("1-026", "1-036")) & (cellType == "B cell")')
    fig, ax = plt.subplots(1, 1, figsize=(2.5, 2.5))
    df = dsip.samplesheet[['patientSample', 'cellType', 'infected', 'numberGenes']].copy()
    df.loc[df['infected'] == True, 'infected'] = 'with DENV RNA'
    df.loc[df['infected'] == False, 'infected'] = 'no DENV RNA'

    sns.violinplot(
            x='infected',
            y='numberGenes',
            data=df,
            inner='quartile',
            ax=ax,
            order=['no DENV RNA', 'with DENV RNA'],
            scale='width',
            palette=sns.color_palette(n_colors=len(cell_types)),
            )
    plt.setp(ax.collections, alpha=.8)
    ax.set_xlabel('')
    ax.set_ylabel('number of genes')
    for tk in ax.get_xticklabels():
        tk.set_rotation(90)
    ax.set_ylim(0, 4000)

    plt.tight_layout(w_pad=0.5, h_pad=0.1)
    fig.savefig('../../figures/qc_number_genes_byDENV_Bcells2636.png', dpi=600)
    fig.savefig('../../figures/qc_number_genes_byDENV_Bcells2636.svg')


    plt.ion()
    plt.show()
