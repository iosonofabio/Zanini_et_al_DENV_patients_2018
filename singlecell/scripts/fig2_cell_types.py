# vim: fdm=indent
'''
author:     Fabio Zanini
date:       27/05/18
content:    Figure 2 shows that we capture and identify many cell types from
            human blood. Minor points:
            - little batch effects (or can be controlled)
            - enrichment factor (look at FACS data)

            Note that this figure should be rather compact but can have
            supplementary figures for support.
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

    pa = argparse.ArgumentParser(description='Scatter patient PBMCs.')
    pa.add_argument('--regenerate_layout', action='store_true',
                    help='Regenerate graph layout')
    pa.add_argument('--no-markers', action='store_true',
                    help='Use only rough markers')
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

    # FIG 2A-F
    layout_name = 'lgl'
    k = 5
    k_slice = 6
    k_threshold = 0.2
    graph_layout_fn = '../../data/fig2_graph_layout_{:}_k{:}_slice{:}_threshold{:.2f}.tsv'.format(
            layout_name,
            k, k_slice, k_threshold)
    if (not os.path.isfile(graph_layout_fn)) or args.regenerate_layout:
        print('Make knn graph')
        knn = ds.graph.lshknn(axis='samples', n_neighbors=k, slice_length=k_slice, threshold=k_threshold)
        edges = list(zip(knn.row, knn.col))
        G = ig.Graph(edges=edges, edge_attrs={'weight': knn.data})

        print('Get graph layout')
        if layout_name == 'lgl':
            layout = G.layout_lgl(maxiter=500)
        elif layout_name == 'fr':
            layout = G.layout_fruchterman_reingold(maxiter=500)
        else:
            raise ValueError('Graph layout name not implemented: {:}'.format(layout_name))
        layout_df = pd.DataFrame(
                data=np.array(layout),
                index=ds.samplenames,
                columns=['dim1', 'dim2'],
                )
        layout_df.to_csv(graph_layout_fn, sep='\t', index=True)

    else:
        print('Load existing layout')
        layout_df = pd.read_csv(graph_layout_fn, sep='\t', index_col=0)

    print('Plot graph layout')
    genes_plot = [
            ('PTPRC', 'CD45'),
            ('TRAC', 'TRAC'),
            ('GZMA', 'GZMA'),
            ('MS4A1', 'CD20'),
            ('CD14', 'CD14'),
            ('log_virus_reads_per_million', 'virus abundance'),
            ]
    cmap = 'viridis'
    x = layout_df.values[:, 0]
    y = layout_df.values[:, 1]
    fig, axs = plt.subplots(2, 3, figsize=(7, 8))
    axs = axs.ravel()
    for ig, ((gname, glabel), ax) in enumerate(zip(genes_plot, axs)):
        if gname in ds.counts.index:
            v = np.log10(0.1 + ds.counts.loc[gname].fillna(0).values)
        else:
            v = ds.samplesheet[gname].values
        v = (v - v.min()) / (v.max() - v.min())
        c = mpl.cm.get_cmap(cmap)(v)

        # Plot higher dots on top
        bins = [0, 0.2, 0.4, 0.6, 0.8, 1.1]
        for ib in range(len(bins) - 1):
            bmin = bins[ib]
            bmax = bins[ib + 1]
            ind = (v >= bmin) & (v < bmax)
            ax.scatter(
                    x[ind], y[ind], c=c[ind],
                    s=20, alpha=0.1 + (0.02 + 0.04 * ('virus' in gname)) * ib,
                    zorder=5 + ib,
                    )

        ax.set_title(glabel)
        ax.set_axis_off()
        ax.text(0.01, 0.98, chr(65 + ig),
                transform=ax.transAxes,
                ha='left',
                va='top',
                fontsize=14)

    ax = fig.add_axes([0.09, 0.348, 0.88, 0.02])
    norm = mpl.colors.Normalize(vmin=0, vmax=1)
    cb = mpl.colorbar.ColorbarBase(
        ax=ax, cmap=cmap, norm=norm,
        orientation="horizontal")
    cb.set_label('Gene/virus expression relative to highest expressing cell')

    # FIG 2G
    print('Plot number of cells for each cell type in various patients')
    n_cell_types = ds.samplesheet.groupby(['experiment', 'cellType']).size().unstack().fillna(0).stack().to_frame()
    n_cell_types.columns = ['n']
    n_cell_types['experiment'] = n_cell_types.index.get_level_values(0)
    n_cell_types['cellType'] = n_cell_types.index.get_level_values(1)
    n_cell_types_plot = n_cell_types.copy()
    n_cell_types_plot['n'] += 0.1
    ax = fig.add_axes([0.16, 0.07, 0.81, 0.18])
    sns.swarmplot(
            data=n_cell_types_plot,
            y='n',
            x='cellType',
            order=['T cell', 'NK cell', 'NKT cell', 'B cell', 'monocyte', 'cDC', 'pDC'],
            ax=ax,
            )
    ax.set_ylabel('number of cells')
    ax.set_ylim(0.09, 1100)
    ax.set_yscale('log')
    ax.grid(True, axis='y')
    ax.set_xlabel('cell type')
    ax.set_yticks([0.1, 1, 10, 100, 1000])
    ax.set_yticklabels(['$0$', '$1$', '$10$', '$10^2$', '$10^3$'])
    fig.text(0.085, 0.275, chr(65 + ig + 1),
            transform=fig.transFigure,
            ha='left',
            va='top',
            fontsize=14)

    plt.tight_layout(rect=(0, 0.34, 1, 1), w_pad=0.1, h_pad=0.1)
    #fig.savefig('../../figures/fig2A-G.png')
    #fig.savefig('../../figures/fig2A-G.svg')

    # SUPPLEMENTARY FIG 10
    print('T cells')
    dsT = ds.query_samples_by_metadata('cellType == "T cell"')
    markers_T = [
            'CD4',
            'CD8A', 'CD8B',
            'CXCR4',
            'CCR4',
            'CCR7',
            'CD27',
            'CD28',
            'IL7R',  # CD127
            'FAS',  # CD95
            'B3GAT1',  # CD57
            'CD44',
            'SELL',  # CD62L
            'IL2RA',  # CD25
            'IL2RB',  # CD122
            'SPN',  # CD43
            'KLRG1',
            'PRF1',
            'GZMA',
            'GZMB',
            'UNC13D',
            'FOXP3',
            'CTLA4',
            'TNFRSF18',
            'CD69',
            'IL2',
            'IL4',
            'IL10',
            'IL15',
            'IFNG',
            'TGFB1',
            'CLC',
            'CD40LG',
            'TRAC',  # alpha receptor constant
            'TRBC1',  # beta receptor constant 1
            'TRBC2',  # beta receptor constant 2
            'TRGC1',  # gamma receptor constant 1
            'TRGC2',  # gamma receptor constant 2
            'TRDC',  # delta receptor
            ]
    dsTpl = dsT.query_features_by_name(markers_T, ignore_missing=True)
    vs = pd.read_csv('../../data/scatter_subtypes_Tcells.tsv', index_col=0, sep='\t')

    genes_select = [None, 'CD4', 'CD8A', 'CXCR4', 'CCR7', 'IL7R', 'SELL', 'PRF1']
    fig, axs = plt.subplots(2, 4, figsize=(7, 4), sharex=True, sharey=True)
    axs = axs.ravel()
    for gene, ax in zip(genes_select, axs):
        kwargs = {}
        if gene is None:
            kwargs['color'] = 'grey'

        dsTpl.plot.scatter_reduced_samples(
                vs,
                color_by=gene,
                color_log=True,
                cmap='viridis',
                ax=ax,
                tight_layout=False,
                alpha=0.3,
                s=10,
                **kwargs,
                )
        ax.set_xlabel('')
        ax.set_ylabel('')
        if gene is not None:
            ax.set_title(gene)

    fig.text(0.5, 0.02, 'dimension 1', ha='center')
    fig.text(0.03, 0.5, 'dimension 2', ha='center', va='center', rotation=90)
    fig.suptitle('T cells')
    plt.tight_layout(rect=[0.025, 0.04, 0.88, 0.96])
    ax = fig.add_axes([0.89, 0.1, 0.015, 0.8])
    norm = mpl.colors.Normalize(vmin=0, vmax=1)
    cb = mpl.colorbar.ColorbarBase(
        ax=ax, cmap=cmap, norm=norm,
        orientation="vertical")
    cb.set_label('relative expression')

    fig.savefig('../../figures/supplementary_fig10.svg')
    fig.savefig('../../figures/supplementary_fig10.png')

    # SUPPLEMENTARY FIG 11
    print('NK cells')
    dsNK = ds.query_samples_by_metadata('cellType == "NK cell"')
    markers_NK = [
            'TRAC',  # needed to show it's not T cells
            'NCAM1',  # CD56
            'FCGR3A',  # CD16
            'SELL',  # CD62L
            'B3GAT1',  # CD57
            'IL2RA',  # CD25
            'IL3RA',
            'IL12RB1',
            'IL12RB2',
            'IL15RA',
            'IL18R1',
            'GZMA',
            'GZMB',
            'SPON2',
            'KLRD1',  # CD94
            'KLRC1',  # CD159a
            'KLRC2',
            'KLRB1',
            'CLEC12A',
            'KLRC3',
            'IL15',
            'KLRK1',
            'KLRC4',
            'CLEC4A',
            'CLEC4D',
            'CLEC1B',
            'KIR2DL1',  # CD158a
            'KIR2DL3',  # CD158b
            'CD69',
            'CXCR6',
            ]
    dsNKpl = dsNK.query_features_by_name(markers_NK, ignore_missing=True)
    vs = pd.read_csv('../../data/scatter_subtypes_NKcells.tsv', index_col=0, sep='\t')

    genes_select = [None, 'NCAM1', 'FCGR3A', 'SELL', 'GZMA', 'KLRC1', 'KLRB1', 'KIR2DL3']
    fig, axs = plt.subplots(2, 4, figsize=(7, 4), sharex=True, sharey=True)
    axs = axs.ravel()
    for gene, ax in zip(genes_select, axs):
        kwargs = {}
        if gene is None:
            kwargs['color'] = 'grey'

        dsNKpl.plot.scatter_reduced_samples(
                vs,
                color_by=gene,
                color_log=True,
                cmap='viridis',
                ax=ax,
                tight_layout=False,
                alpha=0.3,
                s=10,
                **kwargs,
                )
        ax.set_xlabel('')
        ax.set_ylabel('')
        if gene is not None:
            ax.set_title(gene)

    fig.text(0.5, 0.02, 'dimension 1', ha='center')
    fig.text(0.03, 0.5, 'dimension 2', ha='center', va='center', rotation=90)
    fig.suptitle('NK cells')
    plt.tight_layout(rect=[0.025, 0.04, 0.88, 0.96])
    ax = fig.add_axes([0.89, 0.1, 0.015, 0.8])
    norm = mpl.colors.Normalize(vmin=0, vmax=1)
    cb = mpl.colorbar.ColorbarBase(
        ax=ax, cmap=cmap, norm=norm,
        orientation="vertical")
    cb.set_label('relative expression')

    fig.savefig('../../figures/supplementary_fig11.svg')
    fig.savefig('../../figures/supplementary_fig11.png')

    # SUPPLEMENTARY FIG 12
    print('B cells')
    dsB = ds.query_samples_by_metadata('cellType == "B cell"')
    markers_B = [
            'CD19', 'MS4A1',  # B cells
            'TCL1A', 'SYK',
            'SDC1',  # CD138 (memory)
            'PRDM1', 'JCHAIN', 'IRF4', 'XBP1', 'CD38', 'CD27',  # plasmablast markers (thx Derek)
            'HLA-DRA', 'HLA-DRB1',  # MHC II
            'CD69',
            'IGHM', 'IGHD', 'IGHG1', 'IGHG2', 'IGHG3', 'IGHG4', 'IGHE', 'IGHA1', 'IGHA2',  # isotypes
            'AICDA',  # AID (NOTE: not expressed??)
            'MME',  # CD10 (early)
            'CD34', 'CD24',  # transitional?
            'SPN',  # CD43
            'CD40',
            'CD93',
            'AKT1',  # partner of TCL1 and one of Shirit's kinases
            'IL4R', 'AIM2',
            'AQP3',  # activated B cells (and dendritic cells)
            'TYROBP',  # late only
            'PAX5',  # early stages only
            'BANK1',  # BCR-induced mobilization of intracellular Ca2+
            'FCMR',  # antiapoptotic receptor for IGHM
            'MZB1',  # promotes IGHM secretion
            'LTB',  # involved in spleen architechture for immune cells
            'CD22',  # maturation marker, aka siglec2
            'LAT',
            'LAT2',  # NTAL
            'CD80',
            ]
    dsBpl = dsB.query_features_by_name(markers_B, ignore_missing=True)
    vs = pd.read_csv('../../data/scatter_subtypes_Bcells.tsv', index_col=0, sep='\t')

    genes_select = [None, 'MS4A1', 'TCL1A', 'CD40', 'CD69', 'PRDM1', 'TYROBP', 'IGHM', 'IGHG1', 'IGHA1']
    fig, axs = plt.subplots(2, 5, figsize=(9, 4), sharex=True, sharey=True)
    axs = axs.ravel()
    for gene, ax in zip(genes_select, axs):
        kwargs = {}
        if gene is None:
            kwargs['color'] = 'grey'

        dsBpl.plot.scatter_reduced_samples(
                vs,
                color_by=gene,
                color_log=True,
                cmap='viridis',
                ax=ax,
                tight_layout=False,
                alpha=0.3,
                s=10,
                **kwargs,
                )
        ax.set_xlabel('')
        ax.set_ylabel('')
        if gene is not None:
            ax.set_title(gene)

    fig.text(0.5, 0.02, 'dimension 1', ha='center')
    fig.text(0.02, 0.5, 'dimension 2', ha='center', va='center', rotation=90)
    fig.suptitle('B cells')
    plt.tight_layout(rect=[0.025, 0.04, 0.9, 0.96])
    ax = fig.add_axes([0.91, 0.1, 0.015, 0.8])
    norm = mpl.colors.Normalize(vmin=0, vmax=1)
    cb = mpl.colorbar.ColorbarBase(
        ax=ax, cmap=cmap, norm=norm,
        orientation="vertical")
    cb.set_label('relative expression')

    fig.savefig('../../figures/supplementary_fig12.svg')
    fig.savefig('../../figures/supplementary_fig12.png')

    # SUPPLEMENTARY FIG 13
    print('Monocytes')
    dsM = ds.query_samples_by_metadata('cellType == "monocyte"')
    markers_M = [
            'CD14',
            'FCGR3A',  # CD16
            'CX3CR1',
            'ITGAM',
            'CSF1R',
            'CCR2',
            'CCR5',
            'CXCR4',
            'VCAM1',
            'ICAM1',
            'HLA-DRA', 'HLA-DRB1',  # MHC II
            'CD68',
            'MRC1',
            'ADGRE1',  # maybe??
            'IFNGR1',
            'IFNGR2',
            'IL17RA',
            'APOER',
            'TLR1',
            'TLR2',
            'TLR4',
            'TLR6',
            'CD36',
            'CCL5',
            'CXCL1',
            'CXCL10',
            'TNF',
            'IL6',
            'IL10',
            'CCL17',
            'CCL18',
            'CCL22',
            'CCL24',
            ]
    dsMpl = dsM.query_features_by_name(markers_M, ignore_missing=True)
    vs = pd.read_csv('../../data/scatter_subtypes_monocytes.tsv', index_col=0, sep='\t')

    genes_select = [None, 'CD14', 'FCGR3A', 'ICAM1', 'ITGAM', 'CX3CR1', 'CCL5', 'CXCL10', 'IL6', 'TNF']
    fig, axs = plt.subplots(2, 5, figsize=(9, 4), sharex=True, sharey=True)
    axs = axs.ravel()
    for gene, ax in zip(genes_select, axs):
        kwargs = {}
        if gene is None:
            kwargs['color'] = 'grey'

        dsMpl.plot.scatter_reduced_samples(
                vs,
                color_by=gene,
                color_log=True,
                cmap='viridis',
                ax=ax,
                tight_layout=False,
                alpha=0.3,
                s=10,
                **kwargs,
                )
        ax.set_xlabel('')
        ax.set_ylabel('')
        if gene is not None:
            ax.set_title(gene)

    fig.text(0.5, 0.02, 'dimension 1', ha='center')
    fig.text(0.02, 0.5, 'dimension 2', ha='center', va='center', rotation=90)
    fig.suptitle('monocytes')
    plt.tight_layout(rect=[0.025, 0.04, 0.9, 0.96])
    ax = fig.add_axes([0.91, 0.1, 0.015, 0.8])
    norm = mpl.colors.Normalize(vmin=0, vmax=1)
    cb = mpl.colorbar.ColorbarBase(
        ax=ax, cmap=cmap, norm=norm,
        orientation="vertical")
    cb.set_label('relative expression')

    fig.savefig('../../figures/supplementary_fig13.svg')
    fig.savefig('../../figures/supplementary_fig13.png')

    # FIG 2H
    print('All subtype tSNEs together')
    vs_by_celltype = {}
    vs_by_celltype['T cell'] = pd.read_csv('../../data/scatter_subtypes_Tcells.tsv', index_col=0, sep='\t')
    vs_by_celltype['NK cell'] = pd.read_csv('../../data/scatter_subtypes_NKcells.tsv', index_col=0, sep='\t')
    vs_by_celltype['B cell'] = pd.read_csv('../../data/scatter_subtypes_Bcells.tsv', index_col=0, sep='\t')
    vs_by_celltype['monocyte'] = pd.read_csv('../../data/scatter_subtypes_monocytes.tsv', index_col=0, sep='\t')

    fig, axs = plt.subplots(
            1, 4,
            figsize=(7, 2.2),
            gridspec_kw={'width_ratios': [1, 1, 1.5, 1]})
    axs = axs.ravel()
    for ct, ax in zip(['T cell', 'NK cell', 'B cell', 'monocyte'], axs):
        vs = vs_by_celltype[ct]
        x, y = vs.values.T
        ax.scatter(x, y, s=15, color='grey', alpha=0.2)
        ax.set_title(ct+'s', pad=20)
        ax.set_axis_off()
    fig.text(0.5, 0.03, 'dimension 1', ha='center')
    fig.text(0.024, 0.5, 'dimension 2', ha='center', va='center', rotation=90)
    plt.tight_layout(rect=[0.005, 0.04, 1, 0.96])
    #fig.savefig('../../figures/fig2H.png')
    #fig.savefig('../../figures/fig2H.svg')

    plt.ion()
    plt.show()
