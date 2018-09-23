# vim: fdm=indent
'''
author:     Fabio Zanini
date:       27/05/18
content:    Try to predict single cell origin (healthy/dengue/severe dengue).
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


cell_types = ['T cell', 'B cell', 'NK cell', 'NKT cell', 'monocyte', 'pDC', 'cDC']


def regenerate_datasets(ds):
    '''Regenerate cell type specific datasets for faster prototyping'''
    for ct in cell_types:
        #FIXME
        if ct != 'B cell':
            continue
        print(ct)
        dsc = ds.query_samples_by_metadata('cellType == @ct', local_dict=locals())
        dsc.counts.to_csv(
            '../../data/counts_10_10_unique_L1_{:}s.tsv.gz'.format(ct.replace(' ', '')),
            sep='\t',
            index=True,
            compression='gzip',
            )
        dsc.samplesheet['name'] = dsc.samplesheet.index
        dsc.samplesheet.to_csv(
            '../../data/samplesheet_10_10_unique_L1_{:}s.tsv'.format(ct.replace(' ', '')),
            sep='\t',
            index=False,
            )


def max_distance_cumulative(data1, data2):
    """
    Computes the maximal vertical distance between cumulative distributions
    (this is the statistic for KS tests). Code mostly copied from
    scipy.stats.ks_twosamp

    Parameters
    ----------
    data1 : array_like
        First data set
    data2 : array_like
        Second data set
    Returns
    -------
    d : float
        Max distance, i.e. value of the Kolmogorov Smirnov test. Sign is + if
        the cumulative of data1 < the one of data2 at that location, else -.
    x : float
        Value of x where maximal distance d is reached.
    """
    from numpy import ma

    (data1, data2) = (ma.asarray(data1), ma.asarray(data2))
    (n1, n2) = (data1.count(), data2.count())
    mix = ma.concatenate((data1.compressed(), data2.compressed()))
    mixsort = mix.argsort(kind='mergesort')
    csum = np.where(mixsort < n1, 1./n1, -1./n2).cumsum()

    # Check for ties
    if len(np.unique(mix)) < (n1+n2):
        ind = np.r_[np.diff(mix[mixsort]).nonzero()[0], -1]
        csum = csum[ind]
        mixsort = mixsort[ind]

    csumabs = ma.abs(csum)
    i = csumabs.argmax()

    d = csum[i]
    # mixsort[i] contains the index of mix with the max distance
    x = mix[mixsort[i]]

    return (d, x)


def plot_KS(dss, top_genes, distances, n_plots=None):
    if n_plots:
        fig, axs = plt.subplots(n_plots, 1, figsize=(2.1, 0.9 + 1.2 * n_plots), sharex=True, sharey=True)
    else:
        fig, axs = plt.subplots(3, 5, figsize=(8, 4), sharex=True, sharey=True)
        axs = axs.ravel()
    for gene, ax in zip(top_genes, axs):
        x2 = np.sort(dss[True].counts.loc[gene].values)
        y2 = 1.0 - np.linspace(0, 1, len(x2))
        ax.plot(x2, y2, color='darkred', label='s')
        x1 = 0.1 + np.sort(dss[False].counts.loc[gene].values)
        y1 = 1.0 - np.linspace(0, 1, len(x1))
        ax.plot(x1, y1, color='steelblue', label='ns')
        ax.grid(True)
        ax.set_title(gene)

        # Max distance
        xmax = distances[gene][1]
        ax.plot([xmax] * 2, [0, 1.1], color='grey', alpha=0.8, lw=1.5, ls='--')

    axs[0].legend(loc='upper right', fontsize=8)
    axs[0].set_xlim(0.9e-1, 1e5)
    axs[0].set_ylim(-0.02, 1.02)
    axs[0].set_xscale('log')
    if n_plots:
        fig.text(0.06, 0.52, 'fraction of cells > x', ha='center', va='center', rotation=90)
        fig.text(0.52, 0.02, 'counts per million', ha='center')
        plt.tight_layout(rect=(0.06, 0.02, 1, 1), h_pad=0.1)
    else:
        fig.text(0.02, 0.52, 'fraction of cells > x', ha='center', va='center', rotation=90)
        fig.text(0.52, 0.02, 'cpm', ha='center')
        plt.tight_layout(rect=(0.02, 0.02, 1, 1))

    return {
        'fig': fig,
        'axs': axs,
        }


def build_test_model(
        dstrain,
        ct,
        ntops=(5, 9, 15, 21, 26, 31, 41, 51, 61, 81, 101),
        plotKS=False,
        dstest=None,
        ):

    def predict(dst, distances):
        # Let's try to combine these based on the KS thresholds (argmax distance)
        # Consensus voting should be enough
        top_genes = list(distances.keys())
        dstop = dst.query_features_by_name(top_genes)
        vote = []
        for gene in top_genes:
            xmax = distances[gene][1]
            # check over versus underexpression
            sign = bool((1 + np.sign(distances[gene][0])) / 2)

            # Because xmax can be 0, we need strictly larger if the severe is overexpressed
            # The same could be achieved in practice by setting xmax = max(0.01, xmax)
            tally = (dstop.counts.loc[gene] > xmax) ^ (not sign)
            vote.append(tally)
        vote = pd.DataFrame(
                data=np.vstack(vote),
                index=top_genes,
                columns=dstop.samplenames,
                )
        consensus = vote.mean(axis=0) > 0.5
        consensus.name = 'consensus'
        return consensus

    # Test on training data by default
    if dstest is None:
        test_is_train = True
        dstest = dstrain
    else:
        test_is_train = False

    # KS
    dss = dstrain.split(phenotypes='severe_dengue')
    comp = dss[True].compare(dss[False])
    # Note: KS works on maximal cumulative distance, so the fraction of misclassified
    # cells is directly related to the statistics; because the number of observations
    # (cells) is constant across genes, this directly reflects the P value
    #top_genes = comp['P-value'].nsmallest(15).index


    # Try different numbers of genes
    panelA = False
    results = []
    for ntop in ntops:

        #######################
        # Make the model here
        #######################
        top_genes = comp['P-value'].nsmallest(ntop).index

        # Try to get the max distances and locations
        distances = []
        for gene in top_genes:
            x1 = dss[False].counts.loc[gene].values
            x2 = dss[True].counts.loc[gene].values
            (d, xmax) = max_distance_cumulative(x1, x2)
            distances.append((gene, d, xmax))
        distances = {x[0]: (x[1], x[2]) for x in distances}
        #######################

        # Make cute plots
        if plotKS and (not panelA):
            print(ct)
            d = plot_KS(dss, top_genes, distances, n_plots=3)
            fig = d['fig']
            #fig.savefig('../../figures/supplementary_fig5A.png', dpi=600)
            #fig.savefig('../../figures/supplementary_fig5A.svg')
            panelA = True

        #######################
        # Predict here
        #######################
        consensus = predict(dstest, distances)
        #######################

        #######################
        # Test the model here
        #######################
        identity = dstest.samplesheet['severe_dengue'].copy()
        identity.name = 'identity'
        tab = pd.concat([consensus, identity], axis=1)
        tab['tmp'] = 1
        tab = tab.groupby(['consensus', 'identity']).count()['tmp'].unstack().fillna(0)

        # Predictor characteristics
        tpr = 1.0 * tab.loc[True, True] / tab.loc[:, True].sum()
        fpr = 1.0 * tab.loc[True, False] / tab.loc[:, False].sum()

        # Accumulate results
        results.append({
            'cellType': ct,
            'nGenes': ntop,
            'genes': top_genes.tolist(),
            'FPR': fpr,
            'TPR': tpr,
            'nCellsTrain': dstrain.n_samples,
            'nCellsTest': dstest.n_samples,
            'testIsTrain': test_is_train,
            })
        #######################

    results = pd.DataFrame(results).set_index('nGenes', drop=False)
    return results


def build_test_model_allcelltypes(ds, celltypes, **kwargs):
    results_alltypes = []
    for ct in celltypes:
        if len(celltypes) != 1:
            dst = ds.query_samples_by_metadata('cellType == @ct', local_dict=locals())
        else:
            dst = ds

        results = build_test_model(dst, ct, **kwargs)
        results_alltypes.append(results.set_index(['cellType', 'nGenes'], drop=False))
    results_alltypes = pd.concat(results_alltypes, axis=0)
    return results_alltypes


def plot_single_ROC(results, ct):
    results = results.query('cellType == @ct')
    fig, axs = plt.subplots(1, 2, figsize=(5.5, 2.6))
    for iax, ax in enumerate(axs):
        ax.plot(
            results['FPR'],
            results['TPR'],
            'o-',
            color='steelblue',
            markersize=2+3*iax,
            lw=1.5+iax,
            )
        ax.scatter(
            results['FPR'].iloc[[0]],
            results['TPR'].iloc[[0]],
            edgecolor='darkred',
            facecolor='none',
            s=(4+2*iax)**2,
            zorder=5,
            lw=2,
            )
        ax.set_xlabel('FPR')
        ax.grid(True)
    axs[0].plot([0, 1], [0, 1], lw=1.5, color='grey', ls='-')
    axs[0].set_ylabel('TPR')
    axs[0].set_xlim(-0.02, 1.02)
    axs[0].set_ylim(-0.02, 1.02)
    fig.suptitle(ct+'s')
    plt.tight_layout(w_pad=0.05, rect=(0, 0, 1, 0.95))


def find_highest_Lehmann(results_alltypes):
    # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4590288/

    celltypes = np.unique(results_alltypes['cellType'])

    amax = []
    for ct in celltypes:
        results = results_alltypes.query('cellType == @ct')
        x = results['FPR'].values
        y = results['TPR'].values
        # Within the Ansatz y := x^(1.0/a), we get
        # 1.0 / a = log_x (y) = log(y) / log(x)
        # which means:
        # a = log(x) / log(y)
        # but of course we need pseudocounts
        a = np.log(x + 1e-3) / np.log(y + 1e-3)
        i = np.argmax(a)
        am = a[i]
        ntopm = results['nGenes'].values[i]
        amax.append({
            'cellType': ct,
            'amax': am,
            'nGenes': ntopm,
            'FPR': x[i],
            'TPR': y[i],
            })
    amax = pd.DataFrame(amax).set_index('cellType', drop=False)
    return amax


def split_train_test(dst, frac_train=0.90):
    n_train = int(frac_train * dst.n_samples)
    ind = np.arange(dst.n_samples)
    np.random.shuffle(ind)
    ind_test = ind[n_train:]
    dst.samplesheet['is_train'] = True
    dst.samplesheet.loc[dst.samplesheet.index[ind_test], 'is_train'] = False
    dstt = dst.split(phenotypes=['is_train'])
    return dstt


def plot_supplementary_fig5(results_alltypes):
    fig, ax = plt.subplots(1, 1, figsize=(3.8, 3.65))
    colors = sns.color_palette(n_colors=len(args.celltypes))
    for ir, ct in enumerate(celltypes_plot):
        results = results_alltypes.query('cellType == @ct')
        if results['FPR'].values.max() == 0:
            continue
        color = colors[ir]
        ax.plot(
            results['FPR'],
            results['TPR'],
            'o-',
            color=color,
            markersize=2+3*iax,
            lw=1.5+iax,
            label=ct,
            zorder=4,
            )
        ax.scatter(
            results['FPR'].iloc[[0]],
            results['TPR'].iloc[[0]],
            edgecolor='darkred',
            facecolor='none',
            s=(4+2*iax)**2,
            zorder=5,
            lw=2,
            label='',
            )
        ax.plot([0, 1], [0, 1], lw=1.5, color='grey', ls='-', label='')
        xfit = np.linspace(0, 1, 1000)
        for a in [2, 5, 10, 18]:
            yfit = xfit**(1.0 / a)
            ax.plot(xfit, yfit, lw=1.2, color='k', alpha=0.1, label='', zorder=2)

    ax.grid(True)
    ax.set_xlabel('False positive rate')
    ax.set_ylabel('True positive rate')
    ax.legend(loc='lower right', fontsize=9, ncol=2)
    plt.tight_layout()
    return {
        'fig': fig,
        'ax': ax,
        }


# Script
if __name__ == '__main__':

    pa = argparse.ArgumentParser(description='Find viral reads in healthy samples.')
    pa.add_argument(
            '--celltypes', choices=cell_types, nargs='+', default=cell_types,
            help='Restrict to only a cell type (faster)',
            )
    pa.add_argument(
            '--regenerate-datasets',
            action='store_true',
            help='Regenerate cell type specific datasets for convenience')
    args = pa.parse_args()

    print('Load data (L1 normalized)')
    ctn = 'dengue_patients'
    if len(args.celltypes) == 1:
        ctn += '_'+args.celltypes[0].replace(' ', '')+'s'
    ds = Dataset(
            counts_table=ctn,
            samplesheet='dengue_patients',
            featuresheet='dengue_patients',
            )

    if args.regenerate_datasets:
        print('Regenerate datasets for single cell types')
        regenerate_datasets(ds)
        print('Finished regenerating datasets')

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

    # So we want to figure out what genes are best predictors of single cells
    # belonging to severe dengue. This can be done either by distribution level
    # statistics e.g. KS or U or by proper cross-validation, or both.
    results_alltypes = build_test_model_allcelltypes(ds, args.celltypes, plotKS=True)
    for ct in args.celltypes:
        plot_single_ROC(results_alltypes, ct)

    print('Find the highest point in Lehmann curves')
    amax = find_highest_Lehmann(results_alltypes)

    print('Print ROC for each cell type in subplots')
    fig, axs = plt.subplots(
        len(args.celltypes), 2,
        figsize=(4.8, 1.3 + 1.3 * len(args.celltypes)),
        squeeze=False)
    for ir, (ct, axs_row) in enumerate(zip(args.celltypes, axs)):
        results = results_alltypes.query('cellType == @ct')
        for iax, ax in enumerate(axs_row):
            ax.plot(
                results['FPR'],
                results['TPR'],
                'o-',
                color='steelblue',
                markersize=2+3*iax,
                lw=1.5+iax,
                )
            ax.scatter(
                results['FPR'].iloc[[0]],
                results['TPR'].iloc[[0]],
                edgecolor='darkred',
                facecolor='none',
                s=(4+2*iax)**2,
                zorder=5,
                lw=2,
                )
            ax.grid(True)
        axs_row[0].plot([0, 1], [0, 1], lw=1.5, color='grey', ls='-')
        afit = amax.loc[ct, 'amax']
        xfit = np.linspace(0, 1, 1000)
        yfit = xfit**(1.0 / afit)
        axs_row[0].plot(xfit, yfit, lw=1.5, color='black', alpha=0.6, ls='-', zorder=2.5)
        axs_row[0].text(
            0.92, 0.08, '$\\alpha = {:.0f}$'.format(afit), ha='right', va='bottom',
            bbox={'facecolor': 'white', 'boxstyle': 'square'})
        axs_row[0].set_xlim(-0.02, 1.02)
        axs_row[0].set_ylim(-0.02, 1.02)
        axs_row[0].set_ylabel(ct, rotation=0, labelpad=15, ha='right')
        if ir != len(axs) - 1:
            axs_row[0].set_xticklabels([])
        xfit = np.linspace(results['FPR'].values.min(), results['FPR'].values.max(), 1000)
        yfit = xfit**(1.0 / afit)
        axs_row[1].plot(xfit, yfit, lw=1.5, color='black', alpha=0.6, ls='-', zorder=2.5)

    fig.text(0.63, 0.01, 'FPR', ha='center', fontsize=14)
    fig.text(0.04, 0.52, 'TPR', ha='center', va='center', rotation=90, fontsize=14)
    plt.tight_layout(w_pad=0.5, h_pad=0, rect=(0.05, 0.015, 1, 1))

    celltypes_plot = ['T cell', 'NK cell', 'NKT cell', 'B cell', 'monocyte', 'pDC', 'cDC']
    if all(x in args.celltypes for x in celltypes_plot):
        print('Plot everything in one axes, supplementary fig 5')
        d = plot_supplementary_fig5(results_alltypes)
        fig = d['fig']
        #fig.savefig('../../figures/supplementary_fig5C.png', dpi=600)
        #fig.savefig('../../figures/supplementary_fig5C.svg')

    print('Cross validation')
    celltypes = args.celltypes
    cross_validation = []
    for ct in celltypes:
        if len(celltypes) != 1:
            dst = ds.query_samples_by_metadata('cellType == @ct', local_dict=locals())
        else:
            dst = ds

        for ii in range(5):
            print('cross-valdation for cell type {:} #{:}'.format(ct, ii+1))
            # Split train/test
            dstt = split_train_test(dst, frac_train=0.90)

            # Build model
            r = build_test_model(dstt[True], ct, dstest=dstt[False])
            r['nCrossVal'] = ii
            r.set_index(['cellType', 'nGenes', 'nCrossVal'], drop=False, inplace=True)
            cross_validation.append(r)
    cross_validation = pd.concat(cross_validation, axis=0)
    results_alltypes_cv = (cross_validation[['FPR', 'TPR', 'cellType', 'nGenes']]
        .groupby(['cellType', 'nGenes'])
        .mean())

    celltypes_plot = ['T cell', 'NK cell', 'NKT cell', 'B cell', 'monocyte', 'pDC', 'cDC']
    if all(x in args.celltypes for x in celltypes_plot):
        print('Plot everything in one axes, supplementary fig 5')
        d = plot_supplementary_fig5(results_alltypes_cv)
        fig = d['fig']
        fig.savefig('../../figures/supplementary_fig5C.png', dpi=600)
        fig.savefig('../../figures/supplementary_fig5C.svg')


    plt.ion()
    plt.show()
