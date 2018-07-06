# vim: fdm=indent
'''
author:     Fabio Zanini
date:       17/03/18
content:    Look into infected B cells for a figure for the paper.
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
            samplesheet='virus',
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

    print('Restrict to patients 10017016 and 10017022, the ones with infected cells')
    dss = ds.query_samples_by_metadata('experiment in ("10017016", "10017022")')

    # FIG 4A
    print('Plot fraction of infected cells by cell type')
    nvr = 30
    ctypes_order = ['T cell', 'NK cell', 'NKT cell', 'B cell', 'monocyte', 'cDC', 'pDC']
    dss.samplesheet['infected'] = dss.samplesheet['virus_reads_per_million'] >= nvr
    gby = dss.samplesheet[['infected', 'cellType']].groupby('cellType')
    inf_by = gby.mean()['infected']
    inf_by = inf_by[ctypes_order]
    fig, axs = plt.subplots(1, 2, figsize=(6.3, 3), gridspec_kw={'width_ratios': [1.8, 1]})
    ax = axs[0]
    colors = sns.color_palette(n_colors=7)
    inf_by.plot.barh(ax=ax, zorder=10)
    for iy, ct in enumerate(inf_by.index):
        n = gby.sum()['infected'][ct]
        ax.text(
            inf_by[ct] + 0.02 + (0.03 * (inf_by[ct] > 0)), iy,
            '{:}/{:}'.format(str(int(n)), str(int(gby.count().loc[ct]))),
            ha='left', va='center',
            bbox={'edgecolor': colors[iy], 'facecolor': 'white', 'alpha': 0.7, 'lw': 3, 'pad': 5},
            )
    ax.set_ylabel('')
    ax.set_xlim(0, 0.57)
    ax.grid(True)
    ax.set_xticks([0, 0.1, 0.2, 0.3, 0.4, 0.5])
    ax.set_xticklabels(['0', '10',  '20', '30', '40', '50'])
    ax.set_xlabel('% cells infected'.format(nvr))
    ax.set_ylim(6.5, -0.5)

    ax = axs[1]
    data = dss.query_samples_by_metadata('cellType != "unknown"').samplesheet[['cellType', 'log_virus_reads_per_million']]
    sns.violinplot(
            x='log_virus_reads_per_million',
            y='cellType',
            data=data,
            inner=None,
            ax=ax,
            order=ctypes_order,
            scale='width',
            palette=sns.color_palette(n_colors=7),
            )
    plt.setp(ax.collections, alpha=.4)
    sns.swarmplot(
            x='log_virus_reads_per_million',
            y='cellType',
            data=data,
            ax=ax,
            order=ctypes_order,
            palette=sns.color_palette(n_colors=7),
            alpha=0.5,
            s=2,
            )
    ax.set_xlim(-1.1, 5.1)
    ax.set_xticks([-1, 0, 1, 2, 3, 4, 5])
    ax.set_xticklabels(['$0$', '$1$', '$10$', '$10^2$', '$10^3$', '$10^4$', '$10^5$'], va='bottom')
    ax.tick_params(axis='x', which='major', pad=14)
    ax.set_xlabel('Virus reads per million RNA')
    ax.set_ylabel('')
    ax.set_yticklabels([])
    for yt in np.arange(-1, 6):
        ax.plot([yt] * 2, [-1, 8], lw=1, color='grey', alpha=0.5, zorder=0.5)
    fig.text(0.01, 0.98, 'A', ha='left', va='top', fontsize=14)
    plt.tight_layout(w_pad=0.1)

    if False:
        print('Differential expression in infected B cells')
        dsb = dss.query_samples_by_metadata('cellType == "B cell"')
        dsb.samplesheet.loc[:, 'infected'] = 'ambiguous'
        dsb.samplesheet.loc[dss.samplesheet['virus_reads_per_million'] >= nvr, 'infected'] = 'yes'
        dsb.samplesheet.loc[dss.samplesheet['virus_reads_per_million'] == 0, 'infected'] = 'no'

        dsb.counts.log(inplace=True)
        dsi = dsb.split(phenotypes=['infected'])
        comp = dsi['yes'].compare(dsi['no'])
        comp['avg_inf'] = dsi['yes'].counts.mean(axis=1)
        comp['avg_uninf'] = dsi['no'].counts.mean(axis=1)
        comp['fold-change'] = comp['avg_inf'] - comp['avg_uninf']

        df = comp.nsmallest(21, columns=['P-value']).sort_values('fold-change', ascending=False)
        dsbp = dsb.query_samples_by_metadata('infected != "ambiguous"')
        from matplotlib.patches import Rectangle
        fig, axs = plt.subplots(3, 7, sharex=True, sharey=True, figsize=(10, 7))
        axs = axs.ravel()
        for (gname, dfi), ax in zip(df.iterrows(), axs):
            data = dsbp.counts.loc[[gname]].T
            data['infected'] = dsbp.samplesheet['infected']
            with sns.axes_style('whitegrid'):
                sns.violinplot(
                        data=data,
                        x='infected',
                        y=gname,
                        inner='quartile',
                        ax=ax,
                        order=['no', 'yes'],
                        zorder=10,
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
            if dfi['fold-change'] > 0:
                rcol = 'darkred'
            else:
                rcol = 'steelblue'

            ax.add_patch(Rectangle(
                    (0, 0),
                    1, 1,
                    transform=ax.transAxes,
                    lw=5,
                    facecolor='none', edgecolor=rcol,
                    ))

        ax.set_ylim(-1.1, 6.3)
        ax.set_yticks([-1, 0, 1, 2, 3, 4, 5])
        ax.set_yticklabels(['$0$', '$1$', '$10$', '$10^2$', '$10^3$', '$10^4$', '$10^5$'])
        fig.text(0.5, 0.02,
                 'is the cell infected?'.format(nvr),
                 ha='center',
                 )
        fig.text(0.02, 0.58, 'counts per million', rotation=90, ha='center')
        plt.tight_layout(rect=(0.03, 0.03, 1, 1))

    # FIG 4B
    print('Differential expression in infected B cells, selected genes')
    dsb = dss.query_samples_by_metadata('cellType == "B cell"')
    dsb.samplesheet.loc[:, 'infected'] = 'ambiguous'
    dsb.samplesheet.loc[dss.samplesheet['virus_reads_per_million'] >= nvr, 'infected'] = 'yes'
    dsb.samplesheet.loc[dss.samplesheet['virus_reads_per_million'] == 0, 'infected'] = 'no'

    dsb.counts.log(inplace=True)
    dsi = dsb.split(phenotypes=['infected'])
    comp = dsi['yes'].compare(dsi['no'])
    comp['avg_inf'] = dsi['yes'].counts.mean(axis=1)
    comp['avg_uninf'] = dsi['no'].counts.mean(axis=1)
    comp['fold-change'] = comp['avg_inf'] - comp['avg_uninf']

    genes = ['IGHM', 'IGHD', 'TCL1A', 'CXCR4', 'CD69', 'IRF1', 'FCRL1', 'TXNIP']
    dsbp = dsb.query_samples_by_metadata('infected != "ambiguous"')
    from matplotlib.patches import Rectangle
    fig, axs = plt.subplots(2, 4, sharex=True, sharey=True, figsize=(6.3, 4.5))
    axs = axs.ravel()
    for gname, ax in zip(genes, axs):
        dfi = comp.loc[gname]
        data = dsbp.counts.loc[[gname]].T
        data['infected'] = dsbp.samplesheet['infected']
        with sns.axes_style('whitegrid'):
            sns.violinplot(
                    data=data,
                    x='infected',
                    y=gname,
                    inner='quartile',
                    ax=ax,
                    order=['no', 'yes'],
                    zorder=10,
                    )
            plt.setp(ax.collections, alpha=.6)
            sns.swarmplot(
                    data=data,
                    x='infected',
                    y=gname,
                    ax=ax,
                    order=['no', 'yes'],
                    zorder=10,
                    s=1.5,
                    alpha=1.0,
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
        for y in np.arange(-1, 6):
            ax.plot([-1, 3], [y] * 2, lw=1, color='grey', alpha=0.5, zorder=0)
    ax.set_ylim(-1.1, 6.8)
    ax.set_yticks([-1, 0, 1, 2, 3, 4, 5])
    ax.set_yticklabels(['$0$', '$1$', '$10$', '$10^2$', '$10^3$', '$10^4$', '$10^5$'])
    fig.text(0.5, 0.02,
             'is the cell infected?'.format(nvr),
             ha='center',
             )
    fig.text(0.027, 0.63, 'counts per million', rotation=90, ha='center')
    fig.text(0.01, 0.98, 'B', ha='left', va='top', fontsize=14)
    plt.tight_layout(rect=(0.03, 0.03, 1, 1))


    if False:
        print('Differential expression stratified by patient')
        genes_diff = df.index
        fig, axs = plt.subplots(3, 7, sharex=True, sharey=True, figsize=(13, 7))
        axs = axs.ravel()
        for gene, ax in zip(genes_diff, axs):
            dtmp = dsbp.counts.loc[[gene]].T
            dtmp['experiment'] = dsbp.samplesheet['experiment']
            dtmp['infected'] = dsbp.samplesheet['infected']
            sns.violinplot(
                    data=dtmp,
                    y=gene,
                    x='experiment',
                    hue='infected',
                    hue_order=('no', 'yes'),
                    ax=ax,
                    )

            ax.set_ylim(-1, 5.1)
            ax.set_yticks([-1, 0, 1, 2, 3, 4, 5])
            ax.set_yticklabels(['$0$', '$1$', '$10$', '$10^2$', '$10^3$', '$10^4$', '$10^5$'])
            if ax != axs[0]:
                ax.legend_.remove()
            ax.set_title(gene)
            ax.set_ylabel('')
            ax.set_xlabel('')

        fig.text(0.5, 0.02, 'experiment', ha='center')
        plt.tight_layout()


    print('Differential expression of genes that are not patient specific')
    genes_diff = ('CXCR4', 'IGHD', 'IGHM', 'CD69', 'TXNIP', 'VPREB3', 'IRF1',
                  'ZFP36L1', 'NFKBIA', 'VIM', 'ZEB2', 'GPR183')
    fig, axs = plt.subplots(2, 6, sharex=True, sharey=True, figsize=(10, 6))
    axs = axs.ravel()
    for gname, ax in zip(genes_diff, axs):
        dfi = comp.loc[gname]
        data = dsbp.counts.loc[[gname]].T
        data['infected'] = dsbp.samplesheet['infected']
        with sns.axes_style('whitegrid'):
            sns.violinplot(
                    data=data,
                    x='infected',
                    y=gname,
                    inner='quartile',
                    ax=ax,
                    order=['no', 'yes'],
                    zorder=10,
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

        if dfi['fold-change'] > 0:
            rcol = 'darkred'
        else:
            rcol = 'steelblue'

        ax.add_patch(Rectangle(
                (0, 0),
                1, 1,
                transform=ax.transAxes,
                lw=5,
                facecolor='none', edgecolor=rcol,
                ))

    ax.set_ylim(-1.1, 6.3)
    ax.set_yticks([-1, 0, 1, 2, 3, 4, 5])
    ax.set_yticklabels(['$0$', '$1$', '$10$', '$10^2$', '$10^3$', '$10^4$', '$10^5$'])
    fig.text(0.5, 0.02,
             'is the cell infected?'.format(nvr),
             ha='center',
             )
    fig.text(0.02, 0.58, 'counts per million', rotation=90, ha='center')

    plt.tight_layout(rect=(0.03, 0.03, 1, 1))


    if False:
        print('Plot dimensionality reduction of B cells on differentially expressed genes')
        dfdr = comp.nsmallest(100, columns=['P-value'])
        dsdr = dsb.copy()
        dsdr.counts = dsdr.counts.loc[dfdr.index]
        # Exclude ribo/mito and batch effect genes
        good_genes = [fn for fn in dsdr.featurenames if (fn not in ('COL6A3', 'TNR', 'ACACB')) and (not fn.startswith('MT')) and (not fn.startswith('RN'))]
        dsdr.counts = dsdr.counts.loc[good_genes]


        #vs = dsdr.dimensionality.pca(transform=None)
        vs = dsdr.dimensionality.tsne(transform=None, perplexity=20)

        fig, axs = plt.subplots(2, 3, figsize=(8, 6), sharex=True, sharey=True)
        axs = axs.ravel()
        x = vs.iloc[:, 0]
        y = vs.iloc[:, 1]
        for ax, gname in zip(axs, ['log_virus_reads_per_million', 'IGHD', 'experiment', 'SYK', 'TCL1A', 'none']):
            if gname in dsdr.counts.index:
                v = dsdr.counts.loc[gname].values
                ax.set_title(gname)
            elif gname in dsdr.samplesheet.columns:
                v = dsdr.samplesheet[gname]
                if 'virus' in gname:
                    ax.set_title('Virus abundance')
                elif 'exp' in gname:
                    ax.set_title('Patient')
            else:
                v = 'grey'

            if isinstance(v, str):
                color = v
            elif isinstance(v[0], str):
                v_uniques = list(np.unique(v))
                cols = sns.color_palette(n_colors=len(v_uniques))
                color = [cols[v_uniques.index(vi)] for vi in v.values]
            else:
                vmin = v.min()
                vmax = v.max()
                vnorm = (v - vmin) / (vmax - vmin)
                color = mpl.cm.get_cmap('viridis')(vnorm)
            ax.scatter(x, y, s=15, c=color, alpha=0.6)
            #ax.grid(True)
            ax.set_axis_off()
        ax.set_xlabel('PC1')
        ax.set_ylabel('PC2')
        plt.tight_layout()


    print('Plot dimensionality reduction on B cell maturation genes')
    genes_Bmat = [
        'PTPRC',  # CD45
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
        ]

    print('Look at genes that correlate with maturation markers')
    # NOTE: done by hand already
    #corr = dsbp.correlation.correlate_features_features('all', genes_Bmat, method='spearman').fillna(0)
    genes_Bmat_corr = [
            'IL4R', 'AIM2',
            'AQP3',  # activated B cells (and dendritic cells)
            'TYROBP',  # late only
            'PAX5',  # early stages only
            'BANK1',  # BCR-induced mobilization of intracellular Ca2+
            'FCMR',  # antiapoptotic receptor for IGHM
            'MZB1',  # promotes IGHM secretion
            'LTB',  # involved in spleen architechture for immune cells
            'CD22',  # maturation marker, aka siglec2
            ]

    #genes_Bmat = [g for g in genes_Bmat if g in dsb.counts.index]
    genes_Bmat = [g for g in genes_Bmat + genes_Bmat_corr if g in dsb.counts.index]
    dsdr = dsb.copy()
    dsdr.counts = dsdr.counts.loc[genes_Bmat]
    # Exclude ribo/mito and batch effect genes
    good_genes = [fn for fn in dsdr.featurenames if (fn not in ('COL6A3', 'TNR', 'ACACB')) and (not fn.startswith('MT')) and (not fn.startswith('RN'))]
    dsdr.counts = dsdr.counts.loc[good_genes]

    #vs = dsdr.dimensionality.pca(transform=None)
    vs = dsdr.dimensionality.tsne(transform=None, perplexity=20)

    if False:
        #fig, axs = plt.subplots(4, 8, figsize=(15, 8), sharex=True, sharey=True)
        fig, axs = plt.subplots(3, 4, figsize=(9, 7), sharex=True, sharey=True)
        axs = axs.ravel()
        x = vs.iloc[:, 0]
        y = vs.iloc[:, 1]
        feas_plot = [
                #'PTPRC',
                'HLA-DRA',
                #'HLA-DRB1',
                'IGHM',
                'IGHD',
                'IGHG1',
                'MS4A1',
                'SYK',
                'TCL1A',
                'CD27',
                'CD38',
                'log_virus_reads_per_million',
                'experiment',
                 'none',
                ]
        for ax, gname in zip(axs, feas_plot):
            if gname in dsdr.counts.index:
                v = dsdr.counts.loc[gname].values
                ax.set_title(gname)
            elif gname in dsdr.samplesheet.columns:
                v = dsdr.samplesheet[gname]
                if 'virus' in gname:
                    ax.set_title('Virus abundance')
                elif 'exp' in gname:
                    ax.set_title('Patient')
            else:
                v = 'grey'

            if isinstance(v, str):
                color = v
            elif isinstance(v[0], str):
                v_uniques = list(np.unique(v))
                cols = sns.color_palette(n_colors=len(v_uniques))
                color = [cols[v_uniques.index(vi)] for vi in v.values]
            else:
                vmin = v.min()
                vmax = v.max()
                vnorm = (v - vmin) / (vmax - vmin)
                color = mpl.cm.get_cmap('viridis')(vnorm)

            if 'exp' in gname:
                alpha = 0.2
            else:
                alpha = 0.6

            ax.scatter(x, y, s=15, c=color, alpha=alpha)
            #ax.grid(True)
            ax.set_axis_off()
        ax.set_xlabel('dimension 1')
        ax.set_ylabel('dimension 2')
        plt.tight_layout()

    # FIG 4C
    fig, axs = plt.subplots(6, 1, figsize=(2, 9.5), sharex=True, sharey=True)
    axs = axs.ravel()
    x = vs.iloc[:, 0]
    y = vs.iloc[:, 1]
    feas_plot = [
            'log_virus_reads_per_million',
            'MS4A1',
            'JCHAIN',
            'IGHM',
            'TCL1A',
            'TYROBP',
            ]
    for ax, gname in zip(axs, feas_plot):
        if gname in dsdr.counts.index:
            v = dsdr.counts.loc[gname].values
            ax.set_title(gname)
        elif gname in dsdr.samplesheet.columns:
            v = dsdr.samplesheet[gname]
            if 'virus' in gname:
                ax.set_title('Virus abundance')
            elif 'exp' in gname:
                ax.set_title('Patient')
        else:
            v = 'grey'

        if isinstance(v, str):
            color = v
        elif isinstance(v[0], str):
            v_uniques = list(np.unique(v))
            cols = sns.color_palette(n_colors=len(v_uniques))
            color = [cols[v_uniques.index(vi)] for vi in v.values]
        else:
            vmin = v.min()
            vmax = v.max()
            vnorm = (v - vmin) / (vmax - vmin)
            color = mpl.cm.get_cmap('viridis')(vnorm)

        if 'exp' in gname:
            alpha = 0.2
        else:
            alpha = 0.6

        ax.scatter(x, y, s=15, c=color, alpha=alpha)
        #ax.grid(True)
        ax.set_axis_off()
    ax.set_xlabel('dimension 1')
    ax.set_ylabel('dimension 2')
    fig.text(0.01, 0.98, 'C', ha='left', va='top', fontsize=14)
    plt.tight_layout(rect=[-0.15, 0.04, 1, 1])
    ax = fig.add_axes([0.15, 0.07, 0.7, 0.015])
    norm = mpl.colors.Normalize(vmin=0, vmax=1)
    cb = mpl.colorbar.ColorbarBase(
        ax=ax, cmap='viridis', norm=norm,
        orientation="horizontal")
    cb.set_label('Gene/virus expression\n(relative to highest cell)')
    #fig.savefig('../../../papers/dengue_patients/draft_20180527/figures/fig4C.svg')
    #fig.savefig('../../../papers/dengue_patients/draft_20180527/figures/fig4C.png')

    print('Plot many more genes, no good for paper')
    dssp = ds.query_samples_by_metadata('experiment in ("10017016", "10017022")')
    dssp.counts = dssp.counts.loc[:, dsdr.counts.columns]
    dssp.counts.log(inplace=True)
    x = vs.iloc[:, 0]
    y = vs.iloc[:, 1]
    fig, axs = plt.subplots(4, 8, figsize=(15, 8), sharex=True, sharey=True)
    axs = axs.ravel()
    feas_plot = [
            'PTPRC',
            #'TRAC', 'GZMA', 'CD14',
            'HLA-DRA',
            'HLA-DRB1',
            'IGHM',
            'IGHD',
            'IGHG1',
            'IGHG2',
            'IGHG3',
            'CD19',
            'MS4A1',
            'SYK',
            'TCL1A',
            'AKT1',
            'PRDM1',
            'JCHAIN',
            'CD27',
            'CD38',
            'SDC1',
            'IL4R', 'AIM2',
            'AQP3',  # activated B cells (and dendritic cells)
            'TYROBP',  # late only
            'PAX5',  # early stages only
            'BANK1',  # BCR-induced mobilization of intracellular Ca2+
            'FCMR',  # antiapoptotic receptor for IGHM
            'MZB1',  # promotes IGHM secretion
            'LTB',  # involved in spleen architechture for immune cells
            'CD22',  # maturation marker, aka siglec2
            'AXL',
            'CD1C',
            'log_virus_reads_per_million',
            'experiment',
             #'none',
            ]
    for ax, gname in zip(axs, feas_plot):
        if gname in dssp.counts.index:
            v = dssp.counts.loc[gname].values
            ax.set_title(gname)
        elif gname in dsdr.samplesheet.columns:
            v = dssp.samplesheet[gname]
            if 'virus' in gname:
                ax.set_title('Virus abundance')
            elif 'exp' in gname:
                ax.set_title('Patient')
        else:
            v = 'grey'

        if isinstance(v, str):
            color = v
        elif isinstance(v[0], str):
            v_uniques = list(np.unique(v))
            cols = sns.color_palette(n_colors=len(v_uniques))
            color = [cols[v_uniques.index(vi)] for vi in v.values]
        else:
            vmin = v.min()
            vmax = v.max()
            vnorm = (v - vmin) / (vmax - vmin)
            color = mpl.cm.get_cmap('viridis')(vnorm)

        if 'exp' in gname:
            alpha = 0.2
        else:
            alpha = 0.6

        ax.scatter(x, y, s=15, c=color, alpha=alpha)
        #ax.grid(True)
        ax.set_axis_off()
    ax.set_xlabel('dimension 1')
    ax.set_ylabel('dimension 2')
    plt.tight_layout()
    #fig.savefig('../../../papers/dengue_patients/draft_20180527/figures/infected_B_cells_tsne_supplementary.svg')
    #fig.savefig('../../../papers/dengue_patients/draft_20180527/figures/infected_B_cells_tsne_supplementary.png')

    print('Expression of in vitro huh7 candidates on infected B cells')
    from scipy.stats import spearmanr
    genes_vitro = [
            'SEC61A1',
            'SEC61A2',
            'SEC61B',
            'SEC61G',
            'SSR1',
            'SSR2',
            'SSR3',
            #'SSR4',
            'SRP9',
            'SRP14',
            'SEC11A',
            'SEC11C',
            'SPCS1',
            'SPCS2',
            'SPCS3',
            'STT3A',
            'STT3B',
            'OSTC',
            'OST4',
            'DDOST',
            'RPN1',
            'RPN2',
            'DAD1',
            #'MAGT1',
            'KRTCAP2',
            #'TUSC3',
            'TMEM258',
            'DDIT3',
            'TRIB3',
            'HM13',
            'TRAM1',
            'TMED2',
            'RRBP1',
            'COPE',
            'ID2',
            ]
    x = dssp.samplesheet['log_virus_reads_per_million']
    fig, axs = plt.subplots(8, 4, sharex=True, sharey=True, figsize=(7, 11))
    axs = axs.ravel()
    for gname, ax in zip(genes_vitro, axs):
        y = dssp.counts.loc[gname]
        ax.scatter(x, y, s=15, alpha=0.4, zorder=4)
        ax.set_xlim(-1.1, 5.1)
        ax.set_xticks([-1, 0, 1, 2, 3, 4, 5])
        ax.set_xticklabels(['$0$', '$1$', '$10$', '$10^2$', '$10^3$', '$10^4$', '$10^5$'])
        ax.set_ylim(-1.1, 5.1)
        ax.set_yticks([-1, 0, 1, 2, 3, 4, 5])
        ax.set_yticklabels(['$0$', '$1$', '$10$', '$10^2$', '$10^3$', '$10^4$', '$10^5$'])
        for tk in range(-1, 6):
            ax.plot([-2, 6], [tk] * 2, color='grey', lw=1, alpha=0.3, zorder=1)
        for tk in range(-1, 6):
            ax.plot([tk] * 2, [-2, 6], color='grey', lw=1, alpha=0.3, zorder=1)
        ax.set_title(gname)
        rho, pval = spearmanr(x.values, y.values)
        ax.text(0.98, 0.98, '$\\rho = {:.2f}$'.format(rho),
                transform=ax.transAxes,
                ha='right',
                va='top',
                )
    fig.text(0.5, 0.015, 'Virus reads per million', ha='center', va='bottom')
    fig.text(0.025, 0.5, 'Gene counts per million', ha='center', va='center', rotation=90)
    plt.tight_layout(rect=[0.03, 0.025, 1, 1], h_pad=0, w_pad=0)
    #fig.savefig('../../../papers/dengue_patients/draft_20180527/figures/correlation_infected_B_cells_invitro_genes_supplementary.svg')
    #fig.savefig('../../../papers/dengue_patients/draft_20180527/figures/correlation_infected_B_cells_invitro_genes_supplementary.png')

    plt.ion()
    plt.show()
