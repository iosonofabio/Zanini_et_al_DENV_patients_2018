# vim: fdm=indent
'''
author:     Fabio Zanini
date:       08/03/18
content:    Merge counts for all 6 dengue female patients.
'''
# Modules
import os
import sys
import argparse
import numpy as np
import pandas as pd


# Script
if __name__ == '__main__':

    ##################################
    # STEP 1
    ##################################
    expnames = [
            '10017011',  # male, control
            '10017012',  # female, control
            '10017013',  # female, df
            '10017014',  # female, df
            '10017015',  # female, sd
            '10017016',  # female, sd
            '10017017',  # female, control
            '10017018',  # male, control
            '10017021',  # female, sd
            '10017022',  # female, sd
            ]

    #print('Merge raw counts')
    #print('Loading...')
    #dfs = []
    #n_cells = 0
    #for expname in expnames:
    #    print(expname)
    #    df = pd.read_csv(
    #        '../data/sequencing/experiments/{:}/counts.tsv'.format(expname),
    #        sep='\t',
    #        index_col=0,
    #        )
    #    dfs.append(df)
    #    n_cells += df.shape[1]

    #n_genes = df.shape[0]
    #mat = np.zeros((n_genes, n_cells), np.float32)
    #icol = 0
    #cells = []
    #for df in dfs:
    #    mat[:, icol: icol + df.shape[1]] = df.values
    #    cells.extend(df.columns.tolist())
    #    icol += df.shape[1]
    #df_all = pd.DataFrame(
    #        data=mat,
    #        index=df.index,
    #        columns=cells,
    #        )
    #print('Save merged raw counts to file')
    #df_all.to_csv(
    #    '../data/sequencing/datasets/dengue_patients/counts.tsv.gz',
    #    sep='\t',
    #    index=True,
    #    compression='gzip',
    #    )

    #print('Merge virus counts')
    #denv = pd.Series(np.zeros(n_cells, dtype=np.float64), index=df_all.columns)
    #for expname in expnames:
    #    print(expname)
    #    try:
    #        df = pd.read_csv(
    #            '../data/sequencing/experiments/{:}/dengue_counts.tsv'.format(expname),
    #            sep='\t',
    #            index_col=0,
    #            squeeze=True,
    #            )
    #    # Sometimes control patients are missing the data, it's zero anyway
    #    except FileNotFoundError:
    #        continue
    #    denv.loc[df.index] = df
    #denv.to_csv(
    #    '../data/sequencing/datasets/dengue_patients/dengue_counts.tsv',
    #    sep='\t',
    #    index=True,
    #    )

    #print('Free memory')
    #del dfs, df, mat

    #print('Get coverage')
    #coverage = df_all.iloc[:-102].sum(axis=0).astype(np.float64)
    #coverage.to_csv(
    #    '../data/sequencing/datasets/dengue_patients/coverage.tsv',
    #    sep='\t',
    #    index=True,
    #    )

    #fn = '../data/human_genome/GRCh38_ensembl-ERCC/ensemblId_geneName_table.tsv'
    #fs_all = pd.read_csv(fn, sep='\t', index_col=0)
    #fs_all.to_csv(
    #    '../data/sequencing/datasets/dengue_patients/featuresheet.tsv',
    #    sep='\t',
    #    index=True,
    #    )

    #print('Save coverage and number dengue reads to google sheet')
    #from singlecell.googleapi.samplesheet import SampleSheet
    #ss = SampleSheet(sandbox=False)
    #data = ss.get_data(sheetname='sequenced')
    #header = data[0]
    #name_col = header.index('name')
    #cov_col = header.index('coverage')
    #denv_col = header.index('numberDengueReads')
    #for datum in data[1:]:
    #    if datum[name_col] in coverage.index:
    #        if len(datum) <= cov_col:
    #            datum.extend(['' for i in range(cov_col + 1 - len(datum))])
    #        datum[cov_col] = coverage.loc[datum[name_col]]
    #for datum in data[1:]:
    #    if datum[name_col] in denv.index:
    #        datum[denv_col] = denv.loc[datum[name_col]]
    #ss.set_sheet(sheetname='sequenced', values=data)


    ##################################
    # STEP 2
    ##################################
    ## Make a version with less genes
    #print('Reading cell names from virus file')
    #cells = pd.read_csv(
    #    '../data/sequencing/datasets/dengue_patients/dengue_counts.tsv',
    #    sep='\t',
    #    index_col=0,
    #    ).index
    #dtype_dict = {c: np.float32 for c in cells}
    #dtype_dict['Feature'] = str

    #print('Reading merged counts from from file')
    #df_all = pd.read_csv(
    #    '../data/sequencing/datasets/dengue_patients/counts.tsv.gz',
    #    sep='\t',
    #    index_col=0,
    #    compression='gzip',
    #    dtype=dtype_dict,
    #    )
    #print('Reading annotations')
    #fs_all = pd.read_csv(
    #    '../data/sequencing/datasets/dengue_patients/featuresheet.tsv',
    #    sep='\t',
    #    index_col=0,
    #    )

    #print('Select decently expressed genes')
    #ind = (df_all > 10).sum(axis=1) > 10
    #df_sel = df_all.loc[ind]
    #fs_sel = fs_all.loc[ind]

    #print('Filter genes that have more than one EnsemblID')
    #from collections import Counter
    #gdict = Counter(fs_sel['GeneName'].values)
    #glist = [gname for gname, c in gdict.items() if c == 1]
    #ind = [gid for gid, gname in fs_sel['GeneName'].items() if gname in glist]
    #df_sel = df_sel.loc[ind]
    #fs_sel = fs_sel.loc[ind]

    #print('Saving selected gene counts to file')
    #df_sel.to_csv(
    #    '../data/sequencing/datasets/dengue_patients/counts_10_10_unique.tsv.gz',
    #    sep='\t',
    #    index=True,
    #    compression='gzip',
    #    )
    #fs_sel.to_csv(
    #    '../data/sequencing/datasets/dengue_patients/featuresheet_10_10_unique.tsv',
    #    sep='\t',
    #    index=True,
    #    )

    #del df_all

    ##################################
    # STEP 3
    ##################################
    ## Make a normalized version for day-to-day analyses
    #print('Read coverage')
    #coverage = pd.read_csv(
    #    '../data/sequencing/datasets/dengue_patients/coverage.tsv',
    #    sep='\t',
    #    index_col=0,
    #    )

    #print('Read selected gene counts')
    #fn = '../data/sequencing/datasets/female_dengue_patients/counts_10_10_unique.tsv.gz'
    #df_sel = pd.read_csv(
    #    fn,
    #    sep='\t',
    #    index_col=0,
    #    compression='gzip',
    #    )

    #print('Save spikeins and other feature names')
    #spikeins = [gn for gn in df_sel.index[-150:] if gn.startswith('ERCC')]
    #other = df_sel.index[-5:].tolist()
    #with open('../data/sequencing/datasets/dengue_patients/spikeins_10_10_unique.tsv', 'wt') as f:
    #    f.write('\n'.join([str(x) for x in spikeins]))
    #with open('../data/sequencing/datasets/dengue_patients/other_10_10_unique.tsv', 'wt') as f:
    #    f.write('\n'.join([str(x) for x in other]))
    #n_features = df_sel.shape[0] - 5 - len(spikeins)

    #print('L1 normalization')
    #df_L1 = (1e6 * df_sel.iloc[:n_features] / coverage.astype(np.float32)).fillna(0)

    #print('Saving normalized counts to file')
    #fn = '../data/sequencing/datasets/dengue_patients/counts_10_10_unique_L1.tsv.gz'
    #df_L1.to_csv(
    #        fn,
    #        sep='\t',
    #        index=True,
    #        compression='gzip',
    #        )

    #print('Load selected features')
    #fs_sel = pd.read_csv(
    #    '../data/sequencing/datasets/dengue_patients/featuresheet_10_10_unique.tsv',
    #    sep='\t',
    #    index_col=0,
    #    )
    #spikeins = [gn for gn in fs_sel.index[-150:] if gn.startswith('ERCC')]
    #n_features = fs_sel.shape[0] - 5 - len(spikeins)

    #print('Load biomart')
    #fs_L1 = fs_sel.iloc[:n_features]
    #fn_biomart = '../data/sequencing/datasets/dengue_patients/mart_export.tsv'
    #fs_biomart = pd.read_csv(
    #        fn_biomart,
    #        sep='\t',
    #        )

    #print('Drop biomart dups')
    #fs_biomart.drop_duplicates(subset=['Gene stable ID'], inplace=True)
    #fs_biomart.set_index('Gene stable ID', inplace=True)

    #print('Annotate featuresheet based on biomart')
    #ind = np.intersect1d(fs_biomart.index, fs_L1.index)
    #cdict = {
    #        'Chromosome/scaffold name': str,
    #        'Gene start (bp)': int,
    #        'Gene end (bp)': int,
    #        'Strand': str,
    #        'Transcription start site (TSS)': int,
    #        'Gene description': str,
    #        }
    #defaults = {
    #        str: "unknown",
    #        int: -1,
    #        }
    #for key, typ in cdict.items():
    #    fs_L1.loc[:, key] = defaults[typ]
    #    fs_L1.loc[ind, key] = fs_biomart.loc[ind, key].values

    #print('Save featuresheet')
    #fn = '../data/sequencing/datasets/dengue_patients/featuresheet_10_10_unique_L1.tsv'
    #fs_L1.to_csv(
    #        fn,
    #        sep='\t',
    #        index=True,
    #        )

    ##################################
    # STEP 4: intronic counts
    ##################################
    ## FIXME: still missing exps 10017011 and 10017012, gotta call those introns
    #print('Load and taylor intronic counts')
    #fn = '../data/sequencing/datasets/dengue_patients/featuresheet_10_10_unique_L1.tsv'
    #featurenames = pd.read_csv(
    #        fn,
    #        sep='\t',
    #        index_col=0,
    #        ).index
    #featurenamesl = featurenames.tolist()
    #samplenames = pd.read_csv(
    #    '../data/sequencing/datasets/dengue_patients/coverage.tsv',
    #    sep='\t',
    #    index_col=0,
    #    header=None,
    #    ).index
    #samplenamesl = samplenames.tolist()
    #cin = np.zeros((len(featurenames), len(samplenames)), dtype=np.float32)
    #for expname in expnames:
    #    print(expname)
    #    fn = '../data/sequencing/experiments/{:}/counts_intronic.tsv'.format(expname)
    #    if not os.path.isfile(fn):
    #        continue
    #    df = pd.read_csv(
    #        fn,
    #        sep='\t',
    #        index_col=0,
    #        )
    #    fnint = np.intersect1d(featurenames, df.index)
    #    ind = [featurenamesl.index(fn) for fn in fnint]

    #    for sn in df.columns:
    #        print(sn)
    #        i = samplenamesl.index(sn)
    #        cin[ind, i] = df.loc[fnint, sn].astype(np.float32)

    #df_intro = pd.DataFrame(
    #        data=cin,
    #        index=featurenames,
    #        columns=samplenames,
    #        )

    #print('Saving intronic (non-normalized) counts to file')
    #fn = '../data/sequencing/datasets/dengue_patients/counts_10_10_unique_intronic.tsv.gz'
    #df_intro.to_csv(
    #        fn,
    #        sep='\t',
    #        index=True,
    #        compression='gzip',
    #        )
