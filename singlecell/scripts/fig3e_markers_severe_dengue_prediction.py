# vim: fdm=indent
'''
author:     Fabio Zanini
date:       25/05/18
content:    Find markers that are expressed consistently in severe dengue
            patients but not in dengue fever or controls.
'''
# Modules
import os
import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import seaborn as sns

os.environ['SINGLET_CONFIG_FILENAME'] = 'singlet.yml'
sys.path.append('/home/fabio/university/postdoc/singlet')
from singlet.dataset import Dataset, CountsTable, FeatureSheet

# NOTE: to run this script, add the repo folder to your PYTHONPATH
from singlecell.modules.cell_subtypes import cell_subtypes



# Functions
def split_cells_by_subtype_and_patient(ds, cell_subtypes):
    '''Split Dataset of cells by subtype'''
    n_subtype = []
    ds_sub_pats = {}
    for ct, subtypes in cell_subtypes.items():
        if ct != 'all':
            dsct = ds.query_samples_by_metadata('cellType == @ct', local_dict=locals())
        else:
            dsct = ds
        for cst, query in subtypes.items():
            print(ct, cst)
            # FIXME: refine
            if 'isotype' in query:
                dscst = dsct.query_samples_by_metadata(query)
            elif cst == 'all':
                dscst = dsct
            else:
                dscst = dsct.query_samples_by_counts(query)
            print(dscst)
            dspat = dscst.split(phenotypes='experiment')
            for expname, dstmp in dspat.items():
                ds_sub_pats[(ct, cst, expname)] = dstmp
                n_subtype.append({
                    'cellType': ct,
                    'subtype': cst,
                    'experiment': expname,
                    'n': dstmp.n_samples,
                    })
    n_subtype = (
        pd.DataFrame(n_subtype)
          .set_index(['cellType', 'subtype', 'experiment'])
          .unstack()
          .fillna(0)
          .loc[:, 'n']
          .T)

    return {
        'n': n_subtype,
        'datasets': ds_sub_pats,
        }


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
    # Discard HLA
    genes_good &= ~ds.featurenames.str.startswith('HLA')
    # Discard IGHV/IGKV/IGLV
    genes_good &= ~ds.featurenames.str.startswith('IGHV')
    genes_good &= ~ds.featurenames.str.startswith('IGKV')
    genes_good &= ~ds.featurenames.str.startswith('IGLV')
    # Discard TRBV/TRBJ
    genes_good &= ~ds.featurenames.str.startswith('TRBV')
    genes_good &= ~ds.featurenames.str.startswith('TRBJ')
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

    print('Split cells by subtype and patient')
    tmp = split_cells_by_subtype_and_patient(ds, cell_subtypes)
    n_subtype = tmp['n']
    #n_subtype.to_csv('../../data/sequencing/datasets/dengue_patients/n_cells_subtypes_and_patient.tsv', sep='\t', index=True)
    ds_sub_pats = tmp['datasets']

    print('Get patient averages in arbitrary subpopulations and calculare AUCs')
    print('Compare means')
    from scipy.stats import ttest_ind
    idx = []
    pvals = []
    logfold = []
    for ct, subtypes in cell_subtypes.items():
        if ct != 'all':
            dsct = ds.query_samples_by_metadata('cellType == @ct', local_dict=locals())
        else:
            dsct = ds
        for cst, query in subtypes.items():
            print(ct, cst)
            # FIXME: refine
            if 'isotype' in query:
                dscst = dsct.query_samples_by_metadata(query)
            elif cst == 'all':
                dscst = dsct
            else:
                dscst = dsct.query_samples_by_counts(query)

            ds_sub_sd = dscst.split(phenotypes='severe_dengue')
            # Compare means
            exp_neg = ds_sub_sd[False].counts.log().values
            exp_pos = ds_sub_sd[True].counts.log().values
            sta, pval = ttest_ind(exp_neg, exp_pos, axis=1, equal_var=False)
            pval = pd.Series(pval, index=dscst.counts.index)
            logfol = pd.Series(
                    exp_pos.mean(axis=1) - exp_neg.mean(axis=1),
                    index=dscst.counts.index)
            idx.append((ct, cst))
            pvals.append(pval)
            logfold.append(logfol)
    idx = pd.MultiIndex.from_tuples(idx, names=['cellType', 'subtype'])
    pvals = pd.DataFrame(pvals, index=idx).fillna(1).T
    logfold = pd.DataFrame(logfold, index=idx).fillna(0).T

    print('Show best discriminating genes for each subtype')
    def get_log2fold_change(ct, cst, genes):
        return logfold.loc[genes, (ct, cst)] / np.log10(2)
    n_subtype_sd = n_subtype.copy()
    n_subtype_sd['severe_dengue'] = ds.samplesheet[['experiment', 'severe_dengue']].drop_duplicates().set_index('experiment').loc[n_subtype.index, 'severe_dengue']
    #df = pd.DataFrame([pvals.T.stack().values, logfold.T.stack().values], index=['pval', 'logfold'], columns=pvals.T.stack().index).T
    for ic, (ct, cst) in enumerate(pvals.columns):
        pval = pvals[(ct, cst)]
        ptop = pval.nsmallest(20)
        ptop = pd.DataFrame([ptop, logfold.loc[ptop.index, (ct, cst)] / np.log10(2)], index=['P value', 'log2fold'], columns=ptop.index).T
        # Filter only decently overexpressed
        ptop = ptop.loc[ptop['log2fold'] >= 4].iloc[:10]

        nsd = n_subtype_sd.groupby('severe_dengue').sum()[(ct, cst)]
        print('{:} {:}: {:}, {:}'.format(ct, cst, nsd.loc[False], nsd.loc[True]))
        print(ptop)
        print()

    pats_sd = ds.samplesheet[['experiment', 'severe_dengue']].drop_duplicates().set_index('experiment')['severe_dengue']
    sds = pats_sd.values
    pnames = pats_sd.index
    def plot_roc(ct, cst, genes, ax=None, legend=False, labels=None, alpha=0.3):
        from sklearn.metrics import roc_curve
        if ax is None:
            new_ax = True
            fig, ax = plt.subplots(1, 1, figsize=(3, 2))
        else:
            new_ax = False
        cos = []
        for pname in pats_sd.index:
            co = np.log10(0.1 + ds_sub_pats[(ct, cst, pname)].counts.loc[genes]).mean(axis=1)
            cos.append(co)
        cos = np.vstack(cos)

        thresholds = np.linspace(-1, 5, 100)
        colors = sns.color_palette('husl', n_colors=len(genes))
        for ig, gene in enumerate(genes):
            # False positive rate
            xs = [0]
            # True positive rate
            ys = [0]
            for th in thresholds[::-1]:
                n_sd = sds.sum()
                n_healthy = (~sds).sum()
                y = 1.0 * ((cos[:, ig] >= th) & (sds == True)).sum() / n_sd
                x = 1.0 * ((cos[:, ig] >= th) & (sds == False)).sum() / n_healthy
                ys.append(y)
                xs.append(x)
            if labels is None:
                label = gene
            else:
                label = labels[ig]
            xs = np.array(xs)
            ys = np.array(ys)
            # Calculate AUC via trapezoid rule
            auc = (0.5 * np.diff(xs) * (ys[1:] + ys[:-1])).sum()
            print(ct, cst, gene, 'AUC', '{:.2%}'.format(auc))

            # Slight shift to show overlapping lines
            ax.plot(xs + 0.01 * ig, ys, lw=2, label=label, color=colors[ig], alpha=alpha)
        ax.plot(np.linspace(0, 1, 100), np.linspace(0, 1, 100), color='grey')
        if legend:
            ax.legend(loc='lower right', fontsize=8)
        if new_ax:
            ax.set_xlabel('False positive rate')
            ax.set_ylabel('True positive rate')
            ax.set_xlim(-0.05, 1.05)
            ax.set_ylim(-0.05, 1.05)
            plt.tight_layout()
        return ax

    # FIG3E
    print('Plot ROC curves for some more cell types and genes')
    fig, axs = plt.subplots(2, 4, figsize=(8.3, 4.5), sharex=True, sharey=True)
    axs = axs.ravel()
    ctst_pairs = [
            (('T cell', 'all'), ('IFI6', 'ISG15', 'IFI44L')),
            (('NK cell', 'all'), ('IFI6', 'MX1', 'IFI44L')),
            (('B cell', 'all'), ('IFI6', 'IFI44L', 'IFIT3')),
            (('monocyte', 'all'), ('MX1', 'IFIT3', 'RSAD2')),
            (('T cell', 'helper'), ('ISG15', 'IFIT3', 'IFITM3')),
            (('NK cell', 'CD16+'), ('IFI6', 'IFIT3', 'MX1')),
            (('B cell', 'naive'), ('IFITM1', 'MX2', 'UBE2L6')),
            (('monocyte', 'double_positive'), ('CD163', 'IFI27', 'IFIT1')),
            ]
    frame_colors = [sns.color_palette(n_colors=7)[i] for i in [0, 1, 3, 4, 0, 1, 3, 4]]
    for iax, (((ct, cst), genes), ax, frame_color) in enumerate(zip(ctst_pairs, axs, frame_colors)):
        pval = pvals[(ct, cst)]
        log2fold_changes = get_log2fold_change(ct, cst, list(genes))
        labels = ['{:}: {:1.1f}'.format(g, x) for g, x in log2fold_changes.items()]
        plot_roc(ct, cst, list(genes), ax=ax, legend=False, labels=labels, alpha=0.5)
        ax.legend(loc='lower right', fontsize=9)
        if cst == 'all':
            ax.set_title('{:}s'.format(ct), fontsize=11)
        else:
            if cst == 'double_positive':
                ax.set_title('{:}s, CD14+ CD16+'.format(ct), fontsize=9)
            else:
                ax.set_title('{:}s, {:}'.format(ct, cst), fontsize=9)
        ax.set_xlim(-0.05, 1.05)
        ax.set_ylim(-0.05, 1.05)
        ax.set_xticks([0, 0.5, 1])
        ax.set_yticks([0, 0.5, 1])
        rect = Rectangle(
                (0.078 + 0.228 * (iax % 4), 0.53 - 0.46 * (iax // 4)),
                0.218, 0.445,
                transform=fig.transFigure,
                lw=3,
                edgecolor=frame_color,
                facecolor='none',
                zorder=10,
                clip_on=False)
        ax.add_patch(rect)
    fig.text(0.5, 0.02, 'False positive rate', ha='center')
    fig.text(0.02, 0.5, 'True positive rate', ha='center', va='center', rotation=90)
    fig.text(0.01, 0.99, 'E', ha='left', va='top', fontsize=16)
    plt.tight_layout(rect=[0.03, 0.042, 1, 1], w_pad=0.5, h_pad=2.4)
    #fig.savefig('../../figures/fig3E.png', dpi=600)
    #fig.savefig('../../figures/fig3E.svg')

    plt.ion()
    plt.show()
