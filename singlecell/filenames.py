# vim: fdm=indent
'''
author:     Fabio Zanini
date:       20/06/16
content:    Filenames support module
'''
# Modules
import os
import platform
import warnings


# Globals
nodename = platform.node()
# Find your location
if nodename.startswith('sh-'):
    runtime_location = 'SHERLOCK-NODE'
elif 'sherlock' in nodename:
    runtime_location = 'SHERLOCK-FRONTEND'
elif nodename.startswith('compute-'):
    runtime_location = 'SINGLECELL-NODE'
elif 'singlecell' in nodename:
    runtime_location = 'SINGLECELL-FRONTEND'
elif nodename == 'X260':
    runtime_location = 'LAPTOP'
else:
    warinings.warn('runtime location not recognized')
    runtime_location = 'GENERIC'

# Data and code folders
if 'VIRUS_SINGLECELL_ROOT_DATA_FOLDER' in os.environ:
    data_foldername = os.getenv('VIRUS_SINGLECELL_ROOT_DATA_FOLDER').rstrip('/')+'/'
elif 'SHERLOCK' in runtime_location:
    data_foldername = '/oak/stanford/groups/quake/fzanini/sequencing_data/virus_singlecell/'
elif 'SINGLECELL' in runtime_location:
    data_foldername = '/datastorese/fzanini/sequencing_data/virus_singlecell/'
elif 'LAPTOP' in runtime_location:
    data_foldername = '/home/fabio/university/postdoc/virus_singlecell/data/sequencing/'
else:
    raise IOError('Runtime location not recognized')

if 'VIRUS_SINGLECELL_SUPPORT_FOLDER' in os.environ:
    support_foldername = os.getenv('VIRUS_SINGLECELL_SUPPORT_FOLDER').rstrip('/')+'/'
elif 'SHERLOCK' in runtime_location:
    support_foldername = '/oak/stanford/groups/quake/fzanini/postdoc/virus_singlecell/data/'
elif 'SINGLECELL' in runtime_location:
    support_foldername = '/local10G/fzanini/postdoc/virus_singlecell/data/'
elif 'LAPTOP' in runtime_location:
    support_foldername = '/home/fabio/university/postdoc/virus_singlecell/data/'
else:
    raise IOError('Runtime location not recognized')

if 'SHERLOCK' in runtime_location:
    module_foldername = '/oak/stanford/groups/quake/fzanini/postdoc/virus_singlecell/singlecell/'
elif 'SINGLECELL' in runtime_location:
    module_foldername = '/local10G/fzanini/postdoc/virus_singlecell/singlecell/'
elif 'LAPTOP' in runtime_location:
    module_foldername = '/home/fabio/university/postdoc/virus_singlecell/singlecell/'
else:
    raise IOError('Runtime location not recognized')

raw_reads_foldername = data_foldername+'raw/'
muxed_reads_foldername = data_foldername+'muxed/'
samples_foldername = data_foldername+'samples/'
experiments_foldername = data_foldername+'experiments/'
experiments_foldername_small = support_foldername+'sequencing/experiments/'
tables_foldername = support_foldername+'tables/'


# Filename factories
def get_sample_table_filename(kind='sequenced', sandbox=False):
    fn = tables_foldername+'samples-'+kind
    if sandbox:
        fn += '_sandbox'
    fn += '.tsv'
    return fn


def get_muxed_foldername(experiment):
    '''Get foldername of a muxed run to demux'''
    return muxed_reads_foldername+experiment.raw_sequencing_folder+'/'


def get_sample_foldername(sample):
    '''Get foldername of a sample'''
    return samples_foldername+sample['name']+'/'


def get_raw_reads_filenames(sample):
    fdn = raw_reads_foldername+sample['sequencingFolder']
    fn1 = None
    fn2 = None
    for fn in os.listdir(fdn):
        if '_R1_' in fn:
            fn1 = fdn+fn
        elif '_R2_' in fn:
            fn2 = fdn+fn
    if fn1 is None:
        raise IOError('Read 1 file not found')
    if fn2 is None:
        raise IOError('Read 2 file not found')
    return (fn1, fn2)


def get_nextgenmap_foldername(sample, virus):
    return get_sample_foldername(sample)+'nextgenmap_'+virus+'/'


def get_nextgenmap_input_filenames(sample, virus):
    fdn = get_nextgenmap_foldername(sample, virus)
    return (fdn+'unmapped1.fastq', fdn+'unmapped2.fastq')


def get_nextgenmap_output_filename(sample, virus, filtered=True, fmt='sam'):
    fn =  get_nextgenmap_foldername(sample, virus)+'mapped_to_virus'
    if not filtered:
        fn += '_unfiltered'
    fn += '.'+fmt
    return fn


def get_nextgenmap_n_mapped_filename(sample, virus, n=None):
    fn = get_nextgenmap_foldername(sample, virus)
    if n is None:
        for fn_tmp in os.listdir(fn):
            if 'n_mapped_reads_' in fn_tmp:
                break
        else:
            raise IOError('File with number of reads not found: '+sample.name)
        fn += fn_tmp
    else:
        fn += 'n_mapped_reads_'+str(n)
    return fn


def get_pipeline_step_done_filename(sample, step):
    if (step == 'stampy') or ('virus' in step):
        step += '_'+sample.virus
    return get_sample_foldername(sample)+step+'.done'


def get_igblastn_exec_filename():
    return '/oak/stanford/groups/quake/fzanini/programs/ncbi-igblast-1.8.0/bin/igblastn'


def get_stampy_exec_filename():
    return '/oak/stanford/groups/quake/shared/singularity-images/stampy/quakelab-singularity-stampy.img'


def get_stampy_foldername(sample, virus):
    return get_sample_foldername(sample)+'stampy_'+virus+'/'


def get_stampy_prefix(sample, virus):
    return get_stampy_foldername(sample, virus)+'reference'


def get_stampy_index_filename(sample, virus):
    return get_stampy_prefix(sample, virus)+'.stidx'


def get_stampy_hash_filename(sample, virus):
    return get_stampy_prefix(sample, virus)+'.sthash'


def get_stampy_input_filenames(sample, virus):
    fdn = get_stampy_foldername(sample, virus)
    return (fdn+'unmapped1.fastq', fdn+'unmapped2.fastq')


def get_stampy_output_filename(sample, virus, filtered=True):
    fn =  get_stampy_foldername(sample, virus)+'mapped_to_virus'
    if not filtered:
        fn += '_unfiltered'
    fn += '.bam'
    return fn


def get_stampy_n_mapped_filename(sample, virus, n=None):
    fn = get_stampy_foldername(sample, virus)
    if n is None:
        for fn_tmp in os.listdir(fn):
            if 'n_mapped_reads_' in fn_tmp:
                break
        else:
            raise IOError('File with number of reads not found: '+sample.name)
        fn += fn_tmp
    else:
        fn += 'n_mapped_reads_'+str(n)
    return fn


def get_smalt_foldername(sample):
    return get_sample_foldername(sample)+'smalt/'


def get_smalt_dengue_16681_index_prefix(sample):
    fn = get_smalt_foldername(sample)
    return fn+'U87411_strain16681'


def get_smalt_output_filename(sample, only_prefix=False, filtered=True):
    fn = get_smalt_foldername(sample)
    if not only_prefix:
        fn += 'mapped_to_virus'
        if filtered:
            fn += '_filtered'
        fn += '.bam'
    return fn


def get_smalt_n_mapped_filename(sample, n=None):
    fn = get_smalt_foldername(sample)
    if n is None:
        for fn_tmp in os.listdir(fn):
            if 'n_mapped_reads_' in fn_tmp:
                break
        else:
            raise IOError('File with number of reads not found: '+sample.name)
        fn += fn_tmp
    else:
        fn += 'n_mapped_reads_'+str(n)
    return fn


def get_vicuna_foldername(sample):
    return get_sample_foldername(sample)+'vicuna/'


def get_vicuna_config_stub_filename():
    return support_foldername+'vicuna_config.txt'


def get_vicuna_config_filename(sample):
    return get_vicuna_foldername(sample)+'vicuna_config.txt'


def get_vicuna_input_foldername(sample):
    return get_vicuna_foldername(sample)+'input/'


def get_vicuna_output_foldername(sample):
    return get_vicuna_foldername(sample)+'output/'


def get_iva_foldername(sample):
    return get_sample_foldername(sample)+'iva/'


def get_iva_input_filename(sample):
    return get_iva_foldername(sample)+'input_interleaved.fastq.gz'


def get_ERCC_table_filename(fmt='tsv'):
    return support_foldername+'ERCC/ERCC.'+fmt


def get_star_foldername(sample):
    return get_sample_foldername(sample)+'star/'


def get_star_output_filename(sample, only_prefix=False):
    fn = get_star_foldername(sample)
    if not only_prefix:
        fn += 'Aligned.out.bam'
    return fn


def get_star_unaligned_filenames(sample):
    fdn = get_star_foldername(sample)
    return (fdn+'Unmapped.out.mate1', fdn+'Unmapped.out.mate2')


def get_human_genome_foldername():
    return support_foldername+'human_genome/GRCh38_ensembl-ERCC/'


def get_star_genomedir_foldername():
    return get_human_genome_foldername()+'STAR_DIR/'


def get_htseq_count_foldername(sample):
    return get_sample_foldername(sample)+'htseq/'


def get_htseq_count_output_filename(sample):
    return get_htseq_count_foldername(sample)+'counts.tsv'


def get_htseq_count_annotation_filename():
    return get_human_genome_foldername()+'Homo_sapiens.GRCh38.82.ERCC.gtf'


def get_cluster_log_error_filename(sample):
    return get_sample_foldername(sample)+'cluster_error.log'


def get_cluster_log_output_filename(sample):
    return get_sample_foldername(sample)+'cluster_output.log'


def get_ensembl_id_gene_table_filename():
    return get_human_genome_foldername()+'ensemblId_geneName_table.tsv'


def get_ERCC_names_filename():
    return get_human_genome_foldername()+'ERCC_names.tsv'


def get_junk_names_filename():
    return get_human_genome_foldername()+'junk_names.tsv'


def get_experiment_foldername(experiment, location='local'):
    if location == 'local':
        root_fdn = support_foldername+'sequencing/experiments/'
    else:
        root_fdn = experiments_foldername
    return root_fdn+experiment.name+'/'


def get_experiment_counts_filename(
        experiment,
        fmt='tsv',
        location='local'):
    return get_experiment_foldername(experiment, location=location)+'counts.'+fmt


def get_experiment_virus_counts_filename(
        experiment,
        virus,
        fmt='tsv',
        location='local'):
    return get_experiment_foldername(experiment, location=location)+virus+'_counts.'+fmt


def get_virus_coverage_filename(experiment, virus, fmt='pickle'):
    return get_experiment_foldername(experiment)+'coverage_'+virus+'.'+fmt


def get_virus_allele_counts_filename(experiment, virus, fmt='pickle'):
    return get_experiment_foldername(experiment)+'allele_counts_'+virus+'.'+fmt


def get_immgen_foldername():
    return support_foldername+'immgen/'


def get_immgen_human_table_filename():
    return get_immgen_foldername()+'immgen_v1_converted_human.tsv'


def get_immgen_informative_genes_filename():
    return get_immgen_foldername()+'immgen_v1_top_200_overdispersed_converted_human.tsv'


def get_picogreen_foldername():
    return support_foldername+'picogreen/'


def get_picogreen_filename(experiment):
    fdn = get_picogreen_foldername()
    for fn in os.listdir(fdn):
        if fn.startswith(experiment.name):
            return fdn+fn
    raise IOError('File not found')


def get_virus_data_foldername(virus='dengue'):
    return support_foldername+virus+'/'


def get_virus_reference_filename(refname, virus, fmt='fasta'):
    import glob
    fns = glob.glob(get_virus_data_foldername(virus=virus)+'*'+refname+'*.'+fmt)
    if len(fns) != 1:
        raise IOError('File not found or too many found')
    return fns[0]


def get_virus_alignments_foldername(virus='dengue'):
    return get_virus_data_foldername(virus=virus)+'alignments/'


def get_virus_reference_alignment_filename(refname, virus):
    return get_virus_alignments_foldername(virus=virus)+'alignment_'+virus+'_to_'+refname+'.fasta'


def get_virus_multiple_alignment_filename(virus):
    return get_virus_alignments_foldername(virus=virus)+virus+'.fasta'
