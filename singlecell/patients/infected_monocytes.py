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

    print('Differential expression in infected monocytes')
    nvr = 10
    dsb = dss.query_samples_by_metadata('cellType == "monocyte"')
    dsb.samplesheet.loc[:, 'infected'] = 'ambiguous'
    dsb.samplesheet.loc[dss.samplesheet['virus_reads_per_million'] >= nvr, 'infected'] = 'yes'
    dsb.samplesheet.loc[dss.samplesheet['virus_reads_per_million'] == 0, 'infected'] = 'no'

    print('Print number of infected monocytes for each patient')
    print(dsb.samplesheet[['patient', 'infected', 'experiment']].groupby(['patient', 'infected']).count())

    print('Annotate monocyte subtype')
    cell_subtypes =  {
        'classical': '(CD14 >= 100) & (FCGR3A < 100)',
        'nonclassical': '(CD14 < 100) & (FCGR3A >= 100)',
        'double_positive': '(CD14 >= 100) & (FCGR3A >= 100)',
        'all': '',
        }
    for subtype, query in cell_subtypes.items():
        if query:
            snames = dsb.query_samples_by_counts(query).samplenames
            dsb.samplesheet.loc[snames, 'cellSubtype'] = subtype

    print('Print numbers and fractions')
    table = dsb.samplesheet[['infected', 'cellSubtype', 'experiment']].groupby(['infected', 'cellSubtype']).count()['experiment'].unstack().fillna(0).astype(int).loc[['no', 'yes'], ['classical', 'double_positive', 'nonclassical']]
    print(table)

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
                    scale='width',
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

    ## FIG 4B
    #print('Differential expression in infected B cells, selected genes')
    #dsb = dss.query_samples_by_metadata('cellType == "B cell"')
    #dsb.samplesheet.loc[:, 'infected'] = 'ambiguous'
    #dsb.samplesheet.loc[dss.samplesheet['virus_reads_per_million'] >= nvr, 'infected'] = 'yes'
    #dsb.samplesheet.loc[dss.samplesheet['virus_reads_per_million'] == 0, 'infected'] = 'no'

    #dsb.counts.log(inplace=True)
    #dsi = dsb.split(phenotypes=['infected'])
    #comp = dsi['yes'].compare(dsi['no'])
    #comp['avg_inf'] = dsi['yes'].counts.mean(axis=1)
    #comp['avg_uninf'] = dsi['no'].counts.mean(axis=1)
    #comp['fold-change'] = comp['avg_inf'] - comp['avg_uninf']

    #genes = ['IGHM', 'IGHD', 'TCL1A', 'CXCR4', 'CD69', 'IRF1', 'FCRL1', 'TXNIP']
    #dsbp = dsb.query_samples_by_metadata('infected != "ambiguous"')
    #from matplotlib.patches import Rectangle
    #fig, axs = plt.subplots(2, 4, sharex=True, sharey=True, figsize=(6.3, 4.5))
    #axs = axs.ravel()
    #for gname, ax in zip(genes, axs):
    #    dfi = comp.loc[gname]
    #    data = dsbp.counts.loc[[gname]].T
    #    data['infected'] = dsbp.samplesheet['infected']
    #    with sns.axes_style('whitegrid'):
    #        sns.violinplot(
    #                data=data,
    #                x='infected',
    #                y=gname,
    #                inner='quartile',
    #                ax=ax,
    #                order=['no', 'yes'],
    #                zorder=10,
    #                )
    #        plt.setp(ax.collections, alpha=.6)
    #        sns.swarmplot(
    #                data=data,
    #                x='infected',
    #                y=gname,
    #                ax=ax,
    #                order=['no', 'yes'],
    #                zorder=10,
    #                s=1.5,
    #                alpha=1.0,
    #                )
    #    ax.set_title(gname)
    #    ax.set_ylabel('')
    #    ax.set_xlabel('')
    #    ax.text(0.06, 0.96, 'P = {:.0e}'.format(dfi['P-value']),
    #            ha='left',
    #            va='top',
    #            transform=ax.transAxes,
    #            bbox={'edgecolor': 'grey', 'facecolor': 'white', 'alpha': 0.9, 'lw': 1, 'pad': 3},
    #            )
    #    for y in np.arange(-1, 6):
    #        ax.plot([-1, 3], [y] * 2, lw=1, color='grey', alpha=0.5, zorder=0)
    #ax.set_ylim(-1.1, 6.8)
    #ax.set_yticks([-1, 0, 1, 2, 3, 4, 5])
    #ax.set_yticklabels(['$0$', '$1$', '$10$', '$10^2$', '$10^3$', '$10^4$', '$10^5$'])
    #fig.text(0.5, 0.02,
    #         'is the cell infected?'.format(nvr),
    #         ha='center',
    #         )
    #fig.text(0.027, 0.63, 'counts per million', rotation=90, ha='center')
    #fig.text(0.01, 0.98, 'B', ha='left', va='top', fontsize=14)
    #plt.tight_layout(rect=(0.03, 0.03, 1, 1))

    if False:
        print('Plot dimensionality reduction of monocytes on differentially expressed genes')
        dfdr = comp.nsmallest(100, columns=['P-value'])
        dsdr = dsb.copy()
        dsdr.counts = dsdr.counts.loc[dfdr.index]
        # Exclude ribo/mito and batch effect genes
        good_genes = [fn for fn in dsdr.featurenames if (fn not in ('COL6A3', 'TNR', 'ACACB')) and (not fn.startswith('MT')) and (not fn.startswith('RN'))]
        dsdr.counts = dsdr.counts.loc[good_genes]


        #vs = dsdr.dimensionality.pca(transform=None)
        vs = dsdr.dimensionality.tsne(transform=None, perplexity=20)

        genes_plot = ['log_virus_reads_per_million', 'IGHD', 'experiment', 'SYK', 'TCL1A', 'none']
        fig, axs = plt.subplots(2, 3, figsize=(8, 6), sharex=True, sharey=True)
        axs = axs.ravel()
        x = vs.iloc[:, 0]
        y = vs.iloc[:, 1]
        for ax, gname in zip(axs, genes_plot):
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
        ax.set_xlabel('dimension 1')
        ax.set_ylabel('dimension 2')
        plt.tight_layout()


    plt.ion()
    plt.show()
