# vim: fdm=indent
'''
author:     Fabio Zanini
date:       03/06/18
content:    Pipeline for virus mapping within patients AFTER the rough virus
            reads have been identified in the Snakemake pipeline. The thing is
            Snakemake is VERY slow to construct that graph ;-)
'''
import os
import sys
import numpy as np
import pysam
import glob
import subprocess as sp
import shutil
import argparse

from singlecell.filenames import experiments_foldername, get_stampy_exec_filename


def shell(call, env=None):
    if env is None:
        env = os.environ.copy()
    return sp.run(call, check=True, shell=True, env=env)


def pq(query_qualities):
    qstring = ''.join([chr(q + 33) for q in query_qualities])
    return qstring


def rc(seq, qual):
    d = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return (''.join([d[x] for x in seq])[::-1], qual[::-1])


def read_dict(read):
    seq = read.query_sequence
    qual = pq(read.query_qualities)
    # reverse reads in BAM are transformed into positive strand, go back
    if read.is_reverse:
        (seq, qual) = rc(seq, qual)
    return {
        'name': read.qname,
        'seq': seq,
        'qual': qual,
        }



if __name__ == '__main__':

    pa = argparse.ArgumentParser(description='Patient virus mapping pipeline')
    pa.add_argument('--experiments', nargs='+', required=True,
                    help='experiments to process')
    pa.add_argument('--virus', choices=['dengue', 'zika'],
                    default='dengue',
                    help='What virus to remap to')
    args = pa.parse_args()
    virus = args.virus

    for expname in args.experiments:
        print(expname)
        root_fdn = experiments_foldername+expname+'/'
        raw_reads_fn = root_fdn+virus+'_reads.bam'
        raw_reads_fastq_fns = [root_fdn+virus+'_read1.fastq', root_fdn+virus+'_read2.fastq']
        remap_reads_fn = root_fdn+virus+'_remapped.bam'
        reference_fn = root_fdn+virus+'_reference_hybrid.fasta'

        if os.path.isfile(remap_reads_fn):
            print('Remapped already, skip')
            continue

        print('First, make fastqs out of the bam')
        with pysam.AlignmentFile(raw_reads_fn, 'rb') as bamfile,\
             open(raw_reads_fastq_fns[0], 'wt') as fr1,\
             open(raw_reads_fastq_fns[1], 'wt') as fr2:
            fr_out = [fr1, fr2]
            readname = None
            pair = []
            bfs = [[], []]
            for read in bamfile:
                if (read.qname != readname) and (len(pair) == 2):
                    for bf, d in zip(bfs, pair):
                        bf.append('@{:}\n{:}\n+\n{:}\n'.format(
                            d['name'],
                            d['seq'],
                            d['qual']))
                    # Keep buffers from overflowing
                    if len(bfs[0]) > 1000:
                        for bf, fr in zip(bfs, fr_out):
                            fr.write(''.join(bf))
                        bfs = [[], []]
                    pair = [read_dict(read)]
                    readname = read.qname
                elif (read.qname == readname) and (len(pair) == 1):
                    pair.append(read_dict(read))
                    readname = read.qname
                # Special case for the initial line
                elif readname is None:
                    pair.append(read_dict(read))
                    readname = read.qname
                else:
                    raise ValueError('Mwo ya?')
                    
            # Empty buffers
            for bf, fr in zip(bfs, fr_out):
                fr.write(''.join(bf))
            bfs = [[], []]

        print('Remap via stampy')
        output_sam=remap_reads_fn[:-3]+'sam'
        output_index=remap_reads_fn[:-3]+'stidx'
        output_hash=remap_reads_fn[:-3]+'sthash'
        output_prefix_sg='/stampy/'+os.path.basename(output_index[:-6])
        reference_folder=os.path.dirname(reference_fn)
        reference_sg='/stampy_reference/'+os.path.basename(reference_fn)
        input_sg=['/stampy_input/'+os.path.basename(i) for i in raw_reads_fastq_fns]
        output_sam_sg='/stampy/'+os.path.basename(output_sam)
        input_folder=os.path.dirname(raw_reads_fn)
        output_folder=os.path.dirname(output_index)
        stampy=get_stampy_exec_filename()
        stampy_call='singularity run -B '+output_folder+':/stampy -B '+input_folder+':/stampy_input -B '+reference_folder+':/stampy_reference '+stampy
        shell("rm -f {:} {:} {:}".format(output_sam, output_index, output_hash))
        shell(stampy_call+" -G {:} {:}".format(output_prefix_sg, reference_sg))
        shell(stampy_call+" -g {:} -H {:}".format(output_prefix_sg, output_prefix_sg))
        shell(stampy_call+" -g {:} -h {:} -o {:} --inputformat=fastq --substitutionrate=0.05 --sensitive -M {:} {:}".format(output_prefix_sg, output_prefix_sg, output_sam_sg, input_sg[0], input_sg[1]))
        shell("samtools view -bT {:} {:} > {:}".format(reference_fn, output_sam, remap_reads_fn))
        shell("rm {:}".format(output_sam))
