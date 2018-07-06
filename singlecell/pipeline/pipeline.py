#!/usr/bin/env python
# vim: fdm=indent
'''
author:     Fabio Zanini
date:       20/06/16
content:    Analyze the reads: map against human transcriptome (and virus)
            and count the number of reads per gene.
'''
# Modules
import os
import sys
import argparse

import singlecell.filenames as sc
from singlecell.samples import Sample
from singlecell.experiment import Experiment


# Classes/Functions
class Pipeline:
    def __init__(self, samplenames, online=True, ignore_missing=False):
        from singlecell.samples import SampleTable

        data_all = SampleTable.load(online=online).set_index('name', drop=False)
        # Cache to file to avoid restrictions of Google API
        if online:
            data_all.save()

        data = data_all.loc[samplenames]
        self.samplenames = list(samplenames)
        self.data = data
        self.data_all = data_all
        self.ignore_missing = ignore_missing

    def input_files_raw_reads(self, wildcards):
        fns = Sample(self.data.loc[wildcards["sample"]]).get_raw_reads_filenames()
        return fns

    def get_htseq_single_filenames(self):
        return [sc.get_htseq_count_output_filename(Sample(row))
                for sn, row in self.data.iterrows()]

    def get_htseq_single_filenames_by_experiment(self, wildcards):
        ind = self.data_all['experiment'] == wildcards['expname']
        samples = self.data_all.loc[ind]

        fns = []
        for sn, row in samples.iterrows():
            sample = Sample(row)
            if self.ignore_missing:
                try:
                    sample.get_raw_reads_filenames()
                except IOError:
                    continue
            fns.append(sc.get_htseq_count_output_filename(sample))
        return fns

    def get_virus_counts_single_filenames_by_experiment(self, wildcards):
        ind = self.data_all['experiment'] == wildcards['expname']
        samples = self.data_all.loc[ind]

        fns = []
        for virus in ['dengue', 'zika']:
            for sn, row in samples.iterrows():
                sample = Sample(row)
                if self.ignore_missing:
                    try:
                        sample.get_raw_reads_filenames()
                    except IOError:
                        continue
                fns.append(sc.get_stampy_output_filename(sample, virus, filtered=True))
        return fns

    def get_experiment_names(self):
        return list(set(self.data['experiment'].tolist()))

    def get_htseq_merge_filenames(self, fmts):
        '''Merge htseq counts by experiment only for complete experiments'''
        fns = []
        for expname in self.get_experiment_names():
            # Only merge experiments with all the samples
            data_all = self.data_all.loc[self.data_all['experiment'] == expname]
            if not data_all.index.isin(self.data.index).all():
                continue

            exp = Experiment(expname, load_tables=False)
            for location in ['local', 'data']:
                for fmt in fmts:
                    fns.append(exp.get_counts_filename(fmt=fmt, location=location))
        return fns

    def get_virus_merge_filenames(self, viruses, fmts):
        '''Merge virus counts by experiment only for complete experiments'''
        fns = []
        for expname in self.get_experiment_names():
            # Only merge experiments with all the samples
            data_all = self.data_all.loc[self.data_all['experiment'] == expname]
            if not data_all.index.isin(self.data.index).all():
                continue

            exp = Experiment(expname, load_tables=False)
            for location in ['local', 'data']:
                for virus in viruses:
                    for fmt in fmts:
                        fns.append(exp.get_virus_counts_filename(
                            virus,
                            fmt=fmt,
                            location=location))
        return fns
    
    def get_virus_reference(self, wildcards):
        sample = Sample(self.data.loc[wildcards["sample"]])
        virus = wildcards['virus']
        if virus == 'dengue':
            from singlecell.dengue.filenames import get_serotype_reference_filename
            serotype = sample['serotype']
            # Control patients have no serotype
            if str(serotype) == 'nan':
                serotype = 2
            return get_serotype_reference_filename(serotype, fmt='fasta')
        elif virus == 'zika':
            from singlecell.zika.filenames import get_zika_MR766_filename
            return get_zika_MR766_filename(fmt='fasta')
        else:
            raise ValueError('Virus not recognized: '+virus)
    
    def align_unmapped_star_reads(self, wildcards, fns, fns_out):
        import gzip
        from Bio.SeqIO.QualityIO import FastqGeneralIterator as FGI
    
        sample = Sample(self.data.loc[wildcards["sample"]])
    
        # Get set of unmapped reads
        rs = [set(), set()]
        for ifn, fn in enumerate(fns):
            with open(fn, 'rt') as f:
                for read in FGI(f):
                    rs[ifn].add(read[0].split('\t')[0])
        read_both_set = rs[0] & rs[1]
    
        # Make list of unmapped read pairs
        read_both = {key: [None, None] for key in read_both_set}
        for ifn, fn in enumerate(fns):
            with open(fn, 'rt') as f:
                for read in FGI(f):
                    read_name = read[0].split('\t')[0]
                    read_both[read_name][ifn] = read[1:]

        # Write to file
        with gzip.open(fns_out[0], 'wt') as f1, gzip.open(fns_out[1], 'wt') as f2:
            for rname, (r1, r2) in read_both.items():
                f1.write('@'+rname+'\n'+r1[0]+'\n+\n'+r1[1]+'\n')
                f2.write('@'+rname+'\n'+r2[0]+'\n+\n'+r2[1]+'\n')

    def get_cell_type(self, wildcards):
        return self.data.loc[wildcards["sample"], 'cellType']

    def is_B_cell(self, wildcards):
        return self.get_cell_type(wildcards) == 'B'

    def get_samplenames_celltype(self, cell_type):
        return self.data.index[self.data.loc[:, 'cellType'] == cell_type].tolist()

    def get_samplenames_found(self):
        if not self.ignore_missing:
            return self.samplenames
        else:
            samplenames = []
            for _, s in self.data.itersamples():
                try:
                    s.get_raw_reads_filenames()
                except IOError:
                    continue
                samplenames.append(s['name'])
            return samplenames

    def filter_igblast(self, fn_raw, fn_filtered):
        import pandas as pd

        # Split input into single assembled sequences
        with open(fn_raw, 'rt') as f:
            hits = f.read().split('# IGBLASTN')

        out = []
        for hit in hits:
            if not hit:
                continue
            lines = hit.split('\n')[1:-1]
            if not lines[0].startswith('# Query: '):
                raise ValueError('Parsing failed')
            gname = lines[0][len('# Query: '):].split(' ')[0]
            for il, line in enumerate(lines):
                if line.startswith('# V-(D)-J rearrangement summary for query sequence'):
                    break
            else:
                continue
            tops = lines[il+1].split('\t')

            # Sometimes the D is missing (light chains, but not only)
            if len(tops) == 7:
                tops.insert(1, 'N/A')

            # Write only top hits for now, disregarding chances of rubbish (TODO!)
            outi = [gname] + tops
            out.append(outi)

        header = ['name', 'V', 'D', 'J', 'chain_type', 'stop_codon', 'vj_frame', 'productive', 'strand']

        # Write output
        with open(fn_filtered, 'wt') as fout:
            fout.write('\t'.join(header)+'\n')
            fout.write('\n'.join(['\t'.join(l) for l in out]))

    def filter_virus_reads(
            self, wildcards, fn_unfiltered, fn_filtered,
            mapper='nextgenmap',
            min_longest_match=25):
        import pysam
        from singlecell.util import pair_generator

        def has_bad_cigar(read_pair, op='OR', min_longest_match=25):
            '''Flag read pair if the CIGAR is badu

            op (string): Must be 'OR' or 'AND'. With 'OR', the pair is
                         discarded if at least one read is bad. With 'AND',
                         the pair is kept if at least one read is good.
            '''
            if op not in ('OR', 'AND'):
                raise ValueError('op must be "OR" or "AND"')
            bad = [False, False]
            for ir, read in enumerate(read_pair):
                match_longest = 0
                for (bt, bl) in read.cigar:
                    if (bt == 0) and (bl > match_longest):
                        match_longest = bl
                if match_longest < min_longest_match:
                    if op == 'OR':
                        return True
                    bad[ir] = True
            if bad[0] and bad[1]:
                return True
            return False
    
        sample = Sample(self.data.loc[wildcards["sample"]])
        virus = wildcards['virus']
    
        n_mapped = 0

        # Sometimes the bamfile is actually just a placeholder
        if os.stat(fn_unfiltered).st_size != 0:
            with pysam.Samfile(fn_unfiltered, 'rb') as bamfile:
                with pysam.Samfile(fn_filtered, 'wb', template=bamfile) as bamfile_f:
                    for read_pair in pair_generator(bamfile):
                        # filter out unmapped and discordant reads
                        if read_pair[0].is_unmapped or read_pair[1].is_unmapped:
                            continue
                        if (not read_pair[0].is_proper_pair) or (not read_pair[1].is_proper_pair):
                            continue
    
                        # filter out if CIGAR is a mess
                        if has_bad_cigar(
                                read_pair,
                                op='AND',
                                min_longest_match=min_longest_match):
                            continue
    
                        bamfile_f.write(read_pair[0])
                        bamfile_f.write(read_pair[1])
    
                        n_mapped += 1
        else:
            import pathlib
            pathlib.Path(fn_filtered).touch()
    
        # Delete old file with number of reads if present
        try:
            if mapper == 'nextgenmap':
                fn_n_reads_old = sc.get_nextgenmap_n_mapped_filename(sample, virus=virus)
            elif mapper == 'stampy':
                fn_n_reads_old = sc.get_stampy_n_mapped_filename(sample, virus=virus)
            else:
                raise ValueError('Mapper not recognized')
            os.remove(fn_n_reads_old)
        except (IOError, OSError):
            pass
    
        # Write new file with number of reads
        if mapper == 'nextgenmap':
            fn_n_reads = sc.get_nextgenmap_n_mapped_filename(sample, virus=virus, n=n_mapped)
        elif mapper == 'stampy':
            fn_n_reads = sc.get_stampy_n_mapped_filename(sample, virus=virus, n=n_mapped)
        else:
            raise ValueError('Mapper not recognized')
        with open(fn_n_reads, 'w'):
            pass

    def add_unmapped_to_htseq_count(self, wildcards, fn):
        sample = Sample(self.data.loc[wildcards["sample"]])
        with open(sc.get_star_unaligned_filenames(sample)[0], 'rt') as f:
            n_um = sum(1 for line in f) // 4
        with open(fn, 'rt') as f:
            output = f.read().rstrip('\n').split('\n')
        output[-2] = output[-2].split('\t')[0]+'\t'+str(n_um)
        output = '\n'.join(output)+'\n'
        with open(fn, 'wt') as f:
            f.write('feature\tcount\n')
            f.write(output)

    def merge_htseq_counts(self, fns_out):
        import numpy as np
        import pandas as pd

        matrix_created = False
        for i, (sn, row) in enumerate(self.data.iterrows()):
            print('Merge gene counts from cell n. {:}: {:}'.format(i+1, sn))
            if self.ignore_missing:
                try:
                    Sample(row).get_raw_reads_filenames()
                except IOError:
                    print('Missing')
                    continue

            fn = sc.get_htseq_count_output_filename(Sample(row))
            c = pd.read_csv(fn, sep='\t', index_col=0, squeeze=True)
            if not matrix_created:
                nrows = c.shape[0]
                ncols = self.data.shape[0]
                counts = pd.DataFrame(
                        np.zeros((nrows, ncols), int),
                        index=c.index,
                        columns=self.data.index)
                matrix_created = True
            counts.iloc[:, i] = c.values

        print('Save to file...')
        for fn_out in fns_out:
            ext = fn_out.split('.')[-1]
            if ext == 'pickle':
                counts.to_pickle(fn_out)
            elif ext == 'tsv':
                counts.to_csv(fn_out, sep='\t', index=True)
            else:
                raise ValueError('Format for merged count file not recognized')

    def merge_virus_counts(self, fns_out):
        import numpy as np
        import pandas as pd

        for virus in ['dengue', 'zika']:
            print(virus)
            matrix_created = False
            for i, (sn, row) in enumerate(self.data.iterrows()):
                print('Merge virus counts from cell n. {:}: {:}'.format(i+1, sn))
                if self.ignore_missing:
                    try:
                        Sample(row).get_raw_reads_filenames()
                    except IOError:
                        print('Missing')
                        continue
                # NOTE: use stampy as mapper, more reliable
                fn = sc.get_stampy_n_mapped_filename(Sample(row), virus)
                c = int(fn.split('_')[-1])
                print(c)
                if not matrix_created:
                    print('creating matrix')
                    nrows = self.data.shape[0]
                    counts = pd.Series(
                            np.zeros(nrows, int),
                            index=self.data.index)
                    matrix_created = True
                    print('matrix created')
                counts.iloc[i] = c

            print('Save to file...')
            for fn_out in fns_out:
                if virus not in fn_out:
                    continue
                ext = fn_out.split('.')[-1]
                if ext == 'pickle':
                    counts.to_pickle(fn_out)
                elif ext == 'tsv':
                    counts.to_csv(fn_out, sep='\t', index=True)
                else:
                    raise ValueError('Format for merged count file not recognized')


# Script
if __name__ == '__main__':

    pa = argparse.ArgumentParser(description='RNA Seq pipeline')
    g = pa.add_mutually_exclusive_group(required=True)
    g.add_argument('--samples', nargs='+', default=None,
                    help='samples to process')
    g.add_argument('--experiments', nargs='+', default=None,
                    help='experiments to process')
    g.add_argument('--sequencing-runs', nargs='+', default=None,
                    help='sequencing runs to process')
    m = pa.add_mutually_exclusive_group(required=True)
    m.add_argument('--cluster', action='store_true',
                    help='fork to cluster')
    m.add_argument('--cluster-force', action='store_true',
                    help='fork to cluster and force execution')
    m.add_argument('--dry', action='store_true',
                    help='dry run')
    m.add_argument('--unlock', action='store_true',
                    help='remove stale locks (this should never happen?)')
    pa.add_argument('--partitions', default=('quake', 'normal'), nargs='+',
                    choices=['quake', 'normal'])
    vdb = pa.add_mutually_exclusive_group(required=False)
    vdb.add_argument('--verbose', action='store_true',
                     help='verbose snakemake')
    vdb.add_argument('--debug', action='store_true',
                     help='debug mode snakemake')
    pa.add_argument('--ignore-missing', action='store_true',
                    help='Ignore missing fastq.gz samples (testing only!)')

    args = pa.parse_args()

    # Check that samples/experiments/sequencing runs exist
    if args.samples is not None:
        table = None
        samplenames = args.samples
    elif args.experiments is not None:
        from singlecell.samples import SampleTable
        table = SampleTable.load(online=True)
        table = table.loc[table['experiment'].isin(args.experiments)]
        samplenames = table['name'].tolist()
    elif args.sequencing_runs is not None:
        from singlecell.samples import SampleTable
        table = SampleTable.load(online=True)
        table = table.loc[table['sequencingRun'].isin(args.sequencing_runs)]
        samplenames = table['name'].tolist()

    if not len(samplenames):
        raise ValueError('No samples selected')

    if args.ignore_missing:
        if table is None:
            table = SampleTable.load(online=True)
        table = table.loc[table['name'].isin(samplenames)]
        sns = []
        for _, sample in table.itersamples():
            sns.append(sample['name'])
            try:
                sample.get_raw_reads_filenames()
            except IOError:
                print('Skipping sample: {:}'.format(sample['name']))
                continue
        samplenames = ['ignore_missing'] + sns
        input("Press Enter to continue...")

    # Inject samplenames into Snakefile via env vars
    SINGLECELL_SNAKEFILE_SAMPLENAMES = ':'.join(samplenames)

    # Various ways of calling snakemake
    import subprocess as sp
    from singlecell.filenames import module_foldername
    snakefilename = module_foldername+'pipeline/Snakefile'
    call = 'snakemake all --snakefile '+snakefilename
    if args.verbose:
        call += ' --verbose'
    if args.debug:
        call += ' --debug'
    if args.unlock:
        call += ' --unlock'
    else:
        call += ' --cluster "sbatch --ntasks=1 --job-name={params.name} --cpus-per-task={threads} --partition='+','.join(args.partitions)+' --mem={params.mem} -o {log.logout} -e {log.logerr}" --keep-target-files -j 100 -w 100 -k -T -p -r --rerun-incomplete'
        if args.dry:
            call += ' -n'
        elif args.cluster_force:
            call += ' --forceall'

        # If it is not a dry run, we should make sure log folders exist for the cluster
        # NOTE: loading the table twice (here and in the Snakefile) is not the
        # most optimized version, but it's the cleanest
        pipe = Pipeline(samplenames)
        expnames = pipe.get_experiment_names()
        for samplename in samplenames:
          os.makedirs(sc.samples_foldername+samplename+'/log', exist_ok=True)
        for expname in expnames:
          os.makedirs(sc.experiments_foldername+samplename+'/log', exist_ok=True)

    print('SINGLECELL_SNAKEFILE_SAMPLENAMES='+SINGLECELL_SNAKEFILE_SAMPLENAMES+' '+call)
    sp.run(
        call,
        shell=True,
        env=dict(
            os.environ,
            SINGLECELL_SNAKEFILE_SAMPLENAMES=SINGLECELL_SNAKEFILE_SAMPLENAMES),
        check=True)

