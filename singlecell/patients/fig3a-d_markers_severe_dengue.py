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
from singlecell.patients.cell_subtypes import cell_subtypes



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

def split_cells_by_subtype(ds, cell_subtypes):
    '''Split Dataset of cells by subtype'''
    n_subtype = []
    ds_sub = {}
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
            ds_sub[(ct, cst)] = dscst
            n_subtype.append({
                'cellType': ct,
                'subtype': cst,
                'n': dscst.n_samples,
                })
    n_subtype = (
        pd.DataFrame(n_subtype)
          .set_index(['cellType', 'subtype'])
          .fillna(0)
          .loc[:, 'n']
          .T)

    return {
        'n': n_subtype,
        'datasets': ds_sub,
        }




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

    print('Split cells by subtype')
    tmp = split_cells_by_subtype(ds, cell_subtypes)
    n_subtype = tmp['n']
    ds_sub = tmp['datasets']
    #n_subtype.to_csv('../../data/sequencing/datasets/dengue_patients/n_cells_subtypes.tsv', sep='\t', index=True)

    print('Get global picture of overexpression in SD')
    genes_over = []
    comps = []
    for (ct, cst), dsi in ds_sub.items():
        # Skip major cell types
        if cst == 'all':
            continue
        print(ct, cst)
        dsid = dsi.split('severe_dengue')
        # P value of comparison
        comp = dsid[True].compare(dsid[False])
        # Fold change of geometric mean
        exp_neg = dsid[False].counts.log(base=10).values
        exp_pos = dsid[True].counts.log(base=10).values
        logfol = (exp_pos.mean(axis=1) - exp_neg.mean(axis=1)) / np.log10(2)
        comp['log2foldchange'] = logfol
        for col in comp.columns:
            comps.append(((ct, cst, col), comp[col]))
        # Select only overexpressed genes
        genes_over.extend(comp.loc[comp['log2foldchange'] > 5.2].index.tolist())
    genes_overset = sorted(set(genes_over))
    compsset = pd.DataFrame([c[1] for c in comps], index=pd.Index([c[0] for c in comps])).T.loc[genes_overset]
    # Some genes are most likely batch effects
    compsset = compsset.loc[~compsset.index.isin(['TNR', 'KNG1', 'bP-2189O9.2', 'MGC39584'])]

    print('Plot global picture of overexpression in SD')
    # FIG 3A
    from scipy.cluster.hierarchy import linkage, leaves_list
    logfol_set = compsset.loc[:, compsset.columns.get_level_values(2) == 'log2foldchange']
    pval_set = compsset.loc[:, compsset.columns.get_level_values(2) == 'P-value']
    negpval_set = -np.log10(pval_set + 1e-210)
    lnk_genes = linkage(logfol_set.values, method='average', metric='euclidean', optimal_ordering=True)
    lnk_subtypes = linkage(logfol_set.values.T, method='average', metric='euclidean', optimal_ordering=True)
    ll_genes = logfol_set.index[leaves_list(lnk_genes)]
    ll_subtypes = np.array([(c[0], c[1]) for c in logfol_set.columns])[leaves_list(lnk_subtypes)]
    camax = np.abs(logfol_set.values).max()
    cnorm = 1.0 * (logfol_set + camax) / (logfol_set.values.max() + camax)
    snorm = 1.0 * (negpval_set - negpval_set.values.min()) / (negpval_set.values.max() - negpval_set.values.min())
    cmap = sns.diverging_palette(220, 10, sep=80, as_cmap=True)
    data = {'x': [], 'y': [], 's': [], 'c': []}
    for ig, gname in enumerate(ll_genes):
        for ict, (ct, cst) in enumerate(ll_subtypes):
            cnor = cnorm.loc[gname, (ct, cst, 'log2foldchange')]
            snor = snorm.loc[gname, (ct, cst, 'P-value')]
            s = 10 + 50 * snor**2
            c = cmap(cnor)
            data['y'].append(ig)
            data['x'].append(ict)
            data['s'].append(s)
            data['c'].append(c)
    fig, ax = plt.subplots(1, 1, figsize=(5.5, 10.2))
    ax.scatter(data['x'], data['y'], s=data['s'], c=data['c'])
    #ax.set_ylabel('Gene')
    #ax.set_xlabel('Subtype')
    ax.set_yticks(np.arange(len(ll_genes)))
    ax.set_yticklabels(ll_genes, fontsize=8)
    ax.set_xticks(np.arange(len(ll_subtypes)))
    ax.set_xticklabels(
            ['{:}, {:}'.format(ll[0], ll[1]) for ll in ll_subtypes],
            rotation=45, ha='right', fontsize=9)
    ax.set_ylim(-0.5, len(ll_genes) - 0.5)
    ax.set_xlim(-0.5, len(ll_subtypes) - 0.5)
    plt.tight_layout()
    fig.savefig('../../../papers/dengue_patients/draft_20180527/figures/fig3A.svg')
    fig.savefig('../../../papers/dengue_patients/draft_20180527/figures/fig3A.png')

    print('Plot distribution of genes common across cell types')
    # FIG 3B-D
    cell_types = ['T cell', 'NK cell', 'NKT cell', 'B cell', 'monocyte', 'cDC', 'pDC']
    genes = ['IFITM1', 'IFIT3', 'CD163']
    dsg = ds.query_features_by_name(genes)
    fig, axs = plt.subplots(
            len(genes), len(cell_types),
            figsize=(8.3, 5.6),
            sharex=True, sharey=True,
            )
    axs = axs.T
    dsct = dsg.split(phenotypes='cellType')
    for icol, (ct, axcol) in enumerate(zip(cell_types, axs)):
        for irow, (gene, ax) in enumerate(zip(genes, axcol)):
            data = np.log10(0.1 + dsct[ct].counts.loc[[gene]].T)
            data['dengue_severity'] = dsct[ct].samplesheet['dengue_severity']
            sns.violinplot(
                    data=data,
                    y=gene,
                    x='dengue_severity',
                    order=[0, 1, 2],
                    ax=ax,
                    zorder=4,
                    scale='width',
                    palette=["#9b59b6", "#3498db", "#2ecc71"],
                    )
            plt.setp(ax.collections, alpha=.8)
            #sns.swarmplot(
            #        data=data,
            #        y=gene,
            #        x='dengue_severity',
            #        order=[0, 1, 2],
            #        ax=ax,
            #        zorder=4,
            #        s=1.5,
            #        alpha=1.0,
            #        palette=["#9b59b6", "#3498db", "#2ecc71"],
            #        )
            if icol == 0:
                ax.set_ylabel(gene, rotation=0, labelpad=18, fontsize=10)
            else:
                ax.set_ylabel('')
            ax.set_xlabel('')
            if irow == len(genes) - 1:
                ax.set_xticklabels(['CT', 'DF', 'SD'])
            elif irow == 0:
                ax.set_title(ct)
            ax.set_ylim(-1.1, 4.3)
            ax.set_yticks(np.log10(np.array([0.1, 1, 10, 100, 1000, 10000])))
            ax.set_yticklabels(['$0$', '$1$', '$10$', '$10^2$', '$10^3$', '$10^4$'])
            for yi in np.arange(-1, 5):
                ax.plot([-1, 3], [yi, yi], lw=1, color='grey', zorder=0)
        rect = Rectangle(
                (0.1257 + 0.1245 * icol, 0.054),
                0.113, 0.936,
                transform=fig.transFigure,
                lw=3,
                edgecolor=sns.color_palette()[icol],
                facecolor='none',
                zorder=10,
                clip_on=False)
        ax.add_patch(rect)
    fig.text(0.52, 0.02, 'dengue severity (CT/DF/SD)', ha='center')
    fig.text(0.015, 0.5, 'counts per million', ha='center', va='center', rotation=90)
    fig.text(0.01, 0.99, 'B', ha='left', va='top', fontsize=16)
    fig.text(0.01, 0.68, 'C', ha='left', va='top', fontsize=16)
    fig.text(0.01, 0.37, 'D', ha='left', va='top', fontsize=16)
    plt.tight_layout(rect=[0.015, 0.04, 1, 1])
    fig.savefig('../../../papers/dengue_patients/draft_20180527/figures/fig3B-D.svg')
    fig.savefig('../../../papers/dengue_patients/draft_20180527/figures/fig3B-D.png')

    if False:
        print('Plot distribution of genes common across cell types (for a slide)')
        cell_types = ['T cell', 'NK cell', 'NKT cell', 'B cell', 'monocyte', 'cDC', 'pDC']
        genes = ['IFITM1', 'IFIT3', 'CXCL10']
        dsg = ds.query_features_by_name(genes)
        fig, axs = plt.subplots(
                len(genes), len(cell_types),
                figsize=(10, 6),
                sharex=True, sharey=True,
                )
        axs = axs.T
        dsct = dsg.split(phenotypes='cellType')
        for icol, (ct, axcol) in enumerate(zip(cell_types, axs)):
            for irow, (gene, ax) in enumerate(zip(genes, axcol)):
                data = np.log10(0.1 + dsct[ct].counts.loc[[gene]].T)
                data['dengue_severity'] = dsct[ct].samplesheet['dengue_severity']
                sns.violinplot(
                        data=data,
                        y=gene,
                        x='dengue_severity',
                        order=[0, 1, 2],
                        ax=ax,
                        zorder=4,
                        scale='width',
                        palette=["#9b59b6", "#3498db", "#2ecc71"],
                        )
                if icol == 0:
                    ax.set_ylabel(gene, rotation=0, ha='right')
                else:
                    ax.set_ylabel('')
                ax.set_xlabel('')
                if irow == len(genes) - 1:
                    ax.set_xticklabels(['CT', 'DF', 'SD'])
                elif irow == 0:
                    ax.set_title(ct)
                ax.set_ylim(-1.1, 4.3)
                ax.set_yticks(np.log10(np.array([0.1, 1, 10, 100, 1000, 10000])))
                ax.set_yticklabels(['$0$', '$1$', '$10$', '$10^2$', '$10^3$', '$10^4$'])
                for yi in np.arange(-1, 5):
                    ax.plot([-1, 3], [yi, yi], lw=1, color='grey', zorder=0)
            rect = Rectangle(
                    (0.135 + 0.1225 * icol, 0.054),
                    0.117, 0.936,
                    transform=fig.transFigure,
                    lw=3,
                    edgecolor=sns.color_palette()[icol],
                    facecolor='none',
                    zorder=10,
                    clip_on=False)
            ax.add_patch(rect)
        fig.text(0.52, 0.02, 'dengue severity (CT/DF/SD)', ha='center')
        fig.text(0.023, 0.5, 'counts per million', ha='center', va='center', rotation=90)
        plt.tight_layout(rect=[0.03, 0.04, 1, 1])

    if False:
        print('Annotate B cells by isotype')
        isotypes = ('M', 'D', 'G1', 'G2', 'G3', 'G4', 'E', 'A1', 'A2')
        dsB = ds.query_samples_by_metadata('cellType == "B cell"')
        dsB.query_features_by_name(['IGH{:}'.format(i) for i in isotypes], inplace=True)
        isos = [isotypes[row.argmax()] for row in dsB.counts.values.T]
        ds.samplesheet['isotype'] = ''
        ds.samplesheet.loc[dsB.samplenames, 'isotype'] = isos
        from singlecell.googleapi.samplesheet import SampleSheet
        ss = SampleSheet(sandbox=False)
        data = ss.get_data(sheetname='sequenced')
        header = data[0]
        name_col = header.index('name')
        iso_col = header.index('isotype')
        for datum in data[1:]:
            if datum[name_col] in dsB.samplesheet.index:
                if len(datum) <= iso_col:
                    datum.extend(['' for i in range(iso_col + 1 - len(datum))])
                datum[iso_col] = ds.samplesheet.loc[datum[name_col], 'isotype']
        #ss.set_sheet(sheetname='sequenced', values=data)

    plt.ion()
    plt.show()
