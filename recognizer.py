#!/usr/bin/env python
"""
reCOGnizer - a tool for functional annotation with COGs

By João Sequeira

Nov 2019
"""

from argparse import ArgumentParser, ArgumentTypeError
from glob import glob
import os
from pathlib import Path
from shutil import which
from subprocess import run, Popen, PIPE, check_output, DEVNULL
import sys
import numpy as np
import pandas as pd
from tqdm import tqdm
from multiprocessing import Pool, cpu_count, Manager
from time import time, gmtime, strftime
from Bio import Entrez, SeqIO
from requests import get as requests_get
import xml.etree.ElementTree as ET
import re

__version__ = '1.6.6'

Entrez.email = "A.N.Other@example.com"

default_print_command = False        # for debugging purposes


def get_arguments():
    parser = ArgumentParser(
        description="reCOGnizer - a tool for domain based annotation with the CDD database",
        epilog="Input file must be specified.")
    parser.add_argument("-f", "--file", help="Fasta file with protein sequences for annotation")
    parser.add_argument(
        "-t", "--threads", type=int, default=cpu_count() - 2,
        help="Number of threads for reCOGnizer to use [max available - 2]")
    parser.add_argument("--evalue", type=float, default=1e-3, help="Maximum e-value to report annotations for [1e-2]")
    parser.add_argument(
        "--pident", type=float, default=0.0, help="[DEPRECATED] Minimum pident to report annotations for [0]")
    parser.add_argument(
        "-o", "--output", help="Output directory [reCOGnizer_results]", default='reCOGnizer_results')
    parser.add_argument(
        "-dr", "--download-resources", default=False, action="store_true",
        help='If resources for reCOGnizer are not available at "resources_directory" [false]')
    parser.add_argument(
        "-rd", "--resources-directory", default=os.path.expanduser('~/recognizer_resources'),
        help="Output directory for storing databases and other resources [~/recognizer_resources]")
    parser.add_argument(
        "-dbs", "--databases", nargs='+',
        choices=["CDD", "Pfam", "NCBIfam", "Protein_Clusters", "Smart", "TIGRFAM", "COG", "KOG"],
        default=["CDD", "Pfam", "NCBIfam", "Protein_Clusters", "Smart", "TIGRFAM", "COG", "KOG"],
        help="Databases to include in functional annotation [all available]")
    parser.add_argument(
        "-db", "--database",
        help="Basename of database for annotation. If multiple databases, use comma separated list (db1,db2,db3)")
    parser.add_argument(
        "--custom-database", action="store_true", default=False, help="If database was NOT produced by reCOGnizer")
    parser.add_argument(
        "-mts", "--max-target-seqs", type=int, default=20, help="Number of maximum identifications for each protein [1]")
    parser.add_argument(
        "--keep-spaces", action="store_true", default=False,
        help="BLAST ignores sequences IDs after the first space. "
             "This option changes all spaces to underscores to keep the full IDs.")
    parser.add_argument(
        "--no-output-sequences", action="store_true", default=False,
        help="Protein sequences from the FASTA input will be stored in their own column.")
    parser.add_argument(
        "--no-blast-info", action="store_true", default=False,
        help="Information from the alignment will be stored in their own columns.")
    parser.add_argument(
        "--quiet", action="store_true", default=False, help="Don't output download information, used mainly for CI.")
    parser.add_argument(
        "-sd", "--skip-downloaded", action="store_true", default=False,
        help="Skip download of resources identified as already downloaded, used mainly for CI.")
    parser.add_argument('-v', '--version', action='version', version=f'reCOGnizer {__version__}')

    taxArguments = parser.add_argument_group('Taxonomy Arguments')
    taxArguments.add_argument(
        "--tax-file", default=None,
        help="File with taxonomic identification of proteins inputted (TSV). "
             "Must have one line per query, query name on first column, taxid on second.")
    taxArguments.add_argument(
        "--protein-id-col", default='qseqid',
        help="Name of column with protein headers as in supplied FASTA file [qseqid]")
    taxArguments.add_argument(
        "--tax-col", default='Taxonomic identifier (SPECIES)',
        help="Name of column with tax IDs of proteins [Taxonomic identifier (SPECIES)]")
    taxArguments.add_argument(
        "--species-taxids", default=False, action='store_true',
        help="If tax col contains Tax IDs of species (required for running COG taxonomic)")

    args = parser.parse_args()

    args.output = args.output.rstrip('/')
    args.resources_directory = args.resources_directory.rstrip('/')

    if hasattr(args, "file"):
        for directory in [f'{args.output}/{folder}' for folder in ['fasta', 'asn', 'blast', 'rpsbproc', 'tmp']] + [
                f'{args.resources_directory}/dbs']:
            if not os.path.isdir(directory):
                Path(directory).mkdir(parents=True, exist_ok=True)
                print(f'Created {directory}')

    return args


def timed_message(message):
    print(f'{strftime("%Y-%m-%d %H:%M:%S", gmtime())}: {message}')


def run_command(bash_command, print_command=default_print_command, stdout=None, stderr=None):
    if print_command:
        print(bash_command)
    run(bash_command.split(), stdout=stdout, stderr=stderr)


def human_time(seconds):
    days = seconds // 86400
    if days > 0:
        return strftime(f"{days}d%Hh%Mm%Ss", gmtime(seconds))
    return strftime("%Hh%Mm%Ss", gmtime(seconds))


def run_pipe_command(bash_command, file='', mode='w', print_command=default_print_command):
    if print_command:
        print(f'{bash_command}{f" > {file}" if file != "" else ""}')
    if file == '':
        Popen(bash_command, stdin=PIPE, shell=True).communicate()
    elif file == 'PIPE':
        return Popen(bash_command, stdin=PIPE, shell=True, stdout=PIPE).communicate()[0].decode('utf8')
    else:
        with open(file, mode) as output_file:
            Popen(bash_command, stdin=PIPE, shell=True, stdout=output_file).communicate()


def get_tabular_taxonomy(output):
    res = requests_get('https://ftp.uniprot.org/pub/databases/uniprot/current_release/rdf/taxonomy.rdf.xz')
    with open('taxonomy.rdf.xz', 'wb') as f:
        f.write(res.content)
    run_command(f'unxz taxonomy.rdf.xz')
    print('Reading RDF taxonomy')
    root = ET.parse('taxonomy.rdf').getroot()
    elems = root.findall('{http://www.w3.org/1999/02/22-rdf-syntax-ns#}Description')
    with open(output, 'w') as f:
        written = f.write('\t'.join(
            ['taxid', 'name', 'rank', 'parent_taxid']) + '\n')  # assignment to "written" stops output to console
        for elem in tqdm(elems, desc='Converting XML taxonomy.rdf to TSV format'):
            info = [elem.get('{http://www.w3.org/1999/02/22-rdf-syntax-ns#}about').split('/')[-1]]
            scientific_name = elem.find('{http://purl.uniprot.org/core/}scientificName')
            info.append(scientific_name.text if scientific_name is not None else '')
            rank_elem = elem.find('{http://purl.uniprot.org/core/}rank')
            info.append(rank_elem.get('{http://www.w3.org/1999/02/22-rdf-syntax-ns#}resource').split('/')[-1]
                        if rank_elem is not None else '')
            upper_taxon = elem.find('{http://www.w3.org/2000/01/rdf-schema#}subClassOf')
            info.append(upper_taxon.get('{http://www.w3.org/1999/02/22-rdf-syntax-ns#}resource').split('/')[-1]
                        if upper_taxon is not None else '')
            written = f.write('\t'.join(info) + '\n')


def download_resources(directory, quiet=False, skip_downloaded=False):
    if not os.path.isdir(f'{directory}/smps'):
        Path(f'{directory}/smps').mkdir(parents=True, exist_ok=True)
        print(f'Created {directory}/smps')
    for location in [
        # Download CDD
        'ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/cdd.tar.gz',
        'https://ftp.ncbi.nlm.nih.gov/pub/mmdb/cdd/cddid_all.tbl.gz',
        'https://ftp.ncbi.nlm.nih.gov/pub/mmdb/cdd/cdd.info',  # only for versions
        # RPSBPROC
        'https://ftp.ncbi.nih.gov/pub/mmdb/cdd/bitscore_specific.txt',
        'https://ftp.ncbi.nih.gov/pub/mmdb/cdd/cddannot.dat.gz',
        'https://ftp.ncbi.nih.gov/pub/mmdb/cdd/cddannot_generic.dat.gz',
        'https://ftp.ncbi.nih.gov/pub/mmdb/cdd/cddid.tbl.gz',
        'https://ftp.ncbi.nih.gov/pub/mmdb/cdd/cdtrack.txt',
        'https://ftp.ncbi.nih.gov/pub/mmdb/cdd/family_superfamily_links',
        # COG categories
        'ftp.ncbi.nlm.nih.gov/pub/COG/COG2020/data/fun-20.tab',
        'ftp.ncbi.nlm.nih.gov/pub/COG/COG2020/data/cog-20.def.tab',
        # COG2EC
        'https://bitbucket.org/scilifelab-lts/lts-workflows-sm-metagenomics/raw/screening_legacy/'
        'lts_workflows_sm_metagenomics/source/utils/cog2ec.py',
        'http://eggnogdb.embl.de/download/eggnog_4.5/eggnog4.protein_id_conversion.tsv.gz',
        'http://eggnogdb.embl.de/download/eggnog_4.5/data/NOG/NOG.members.tsv.gz',
        # NCBIfam, TIGRFAM, Pfam, PRK (protein clusters)
        'https://ftp.ncbi.nlm.nih.gov/hmm/4.0/hmm_PGAP.tsv',
        # SMART
        'https://smart.embl.de/smart/descriptions.pl',
        # KOG
        'https://ftp.ncbi.nlm.nih.gov/pub/COG/KOG/kog'
    ]:
        if os.path.isfile(f"{directory}/{location.split('/')[-1]}"):
            if not skip_downloaded:
                download = str2bool(input(f"{directory}/{location.split('/')[-1]} exists. Overwrite? [Y/N] "))
            else:
                download = False
        else:
            download = True
        if download:
            if quiet:
                print(f'Downloading {location}')        # quiet replaces all wget output with this simple message
            run_command(f'wget {location} -P {directory}{" -q" if quiet else ""}')

    if os.path.isfile(f'{directory}/cdd.tar.gz'):
        os.rename(f'{directory}/cdd.tar.gz', f'{directory}/smps/cdd.tar.gz')

    for file in [
        'cddid_all.tbl', 'eggnog4.protein_id_conversion.tsv', 'NOG.members.tsv', 'cddannot.dat',
            'cddannot_generic.dat', 'cddid.tbl']:
        run_command(f'gunzip {directory}/{file}.gz', print_command=True)

    # Extract the smps
    if sys.platform == "darwin":
        if which('gtar') is None:
            run_command('brew install gnu-tar', print_command=True)
        tool = 'gtar'
    else:
        tool = 'tar'
    wd = os.getcwd()
    os.chdir(f'{directory}/smps')
    run_pipe_command(f'{tool} -xzf cdd.tar.gz --wildcards "*.smp"', print_command=True)
    os.remove('cdd.tar.gz')
    os.chdir(wd)
    get_tabular_taxonomy(f'{directory}/taxonomy.tsv')
    with open(f'{directory}/recognizer_dwnl.timestamp', 'w') as f:
        f.write(strftime("%Y-%m-%d %H:%M:%S", gmtime()))

    generate_cog2ec_df(
        f'{directory}/eggnog4.protein_id_conversion.tsv', f'{directory}/NOG.members.tsv', directory)


def str2bool(v):
    if v.lower() == 'auto':
        return 'auto'
    elif v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise ArgumentTypeError('Boolean value expected.')


def run_rpsblast(query, output, reference, threads='0', max_target_seqs=1, evalue=1e-2):
    # This run_command is different because of reference, which can't be split by space
    bash_command = f'rpsblast -query {query} -db "{reference}" -out {output} -outfmt 11 -num_threads {threads} ' \
                   f'-max_target_seqs {max_target_seqs} -evalue {evalue} 1>recognizer.log 2>recognizer.log'
    run_pipe_command(bash_command)


def parse_cddid(cddid):
    cddid = pd.read_csv(cddid, sep='\t', header=None)[[0, 1, 3]]
    cddid.columns = ['CDD ID', 'DB ID', 'DB description']
    # cddid['CDD ID'] = [f'CDD:{ide}' for ide in cddid['CDD ID']]    # for now, seems to no longer be required
    return cddid


def expand_by_list_column(df, column='COG functional category (letter)'):
    lens = [len(item) for item in df[column]]
    dictionary = dict()
    for col in df.columns:
        dictionary[col] = np.repeat(df[col].values, lens)
    dictionary[column] = np.concatenate(df[column].values)
    return pd.DataFrame(dictionary)


def parse_whog(whog):
    df = pd.read_csv(whog, sep='\t', usecols=[0, 1, 2], header=None, encoding='ISO 8859-1')
    df.columns = ['cog', 'COG functional category (letter)', 'COG protein description']
    df['COG functional category (letter)'] = df['COG functional category (letter)'].apply(lambda x: [i for i in x])
    df = expand_by_list_column(df, column='COG functional category (letter)')
    return df


def parse_kog(kog):
    lines = list()
    for line in [line.rstrip('\n') for line in open(kog).readlines() if line.startswith('[')]:
        line = line.split()
        lines.append([line[0][1], line[1], ' '.join(line[2:])])
    df = pd.DataFrame(lines)
    df.columns = ['KOG functional category (letter)', 'kog', 'KOG protein description']
    df['KOG functional category (letter)'] = df['KOG functional category (letter)'].apply(lambda x: [i for i in x])
    df = expand_by_list_column(df, column='KOG functional category (letter)')
    return df


def parse_blast(file):
    if os.stat(file).st_size != 0:
        blast = pd.read_csv(file, sep='\t', header=None)
        blast.columns = [
            'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue',
            'bitscore']
        return blast
    return pd.DataFrame(columns=[
        'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue',
        'bitscore'])


def pn2database(pn):
    run_command(f"makeprofiledb -in {pn} -title {pn.split('.pn')[0]} -out {pn.split('.pn')[0]}")


def split(a, n):
    k, m = divmod(len(a), n)
    return (a[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(n))


def get_upper_taxids(taxid, tax_df):
    """
    :param taxid: str - taxID to get upper taxIDs from
    :param tax_df: pd.DataFrame - of read taxonomy.tsv (from taxonomy.rdf)
    :returns list - of upper taxIDs
    """
    if taxid == '0':
        return []
    taxids = []
    while taxid != '1' and taxid != 'Taxon':
        taxids.append(taxid)
        taxid = tax_df.loc[taxid]['parent_taxid']
    return taxids


def get_lineages(taxids, taxonomy_df):
    lineages = {}
    all_taxids = []
    for taxid in taxids:
        lineage = get_upper_taxids(taxid, taxonomy_df)
        lineages[taxid] = lineage
        all_taxids += lineage
    return lineages, all_taxids


def get_lineages_multiprocessing(taxids, taxonomy_df, threads=14):
    timed_message(f'Listing all parent tax IDs for {len(taxids)} tax IDs (this may take a while, time for coffee?)')
    taxids_groups = split(list(taxids), threads)
    lineages, res_taxids = ({}, [])
    with Manager() as m:
        with m.Pool() as p:
            result = p.starmap(get_lineages, [(taxids_group, taxonomy_df) for taxids_group in taxids_groups])
    for res in result:
        lineages = {**lineages, **res[0]}
        res_taxids += res[1]
    return lineages, res_taxids


def create_tax_db(smp_directory, db_directory, db_prefix, taxids, hmm_pgap):
    """
    Creates HMM DBs for all required tax IDs, and checks for DBS for cellular organisms and 0 (nan)
    :param smp_directory: (str) - Name of folder with the SMP files
    :param output: (str) - Name of folder to output PN files and databases
    :param db_prefix: (str) - Filename prefix for PN files and databases
    :param taxids: (list) - list of tax ids present in the dataset lacking db
    :param hmm_pgap: (pandas.DataFrame) - df with the information from the hmm_GAP.tsv file
    """
    taxids_with_db = list()
    if len(taxids) == 0:
        return []
    for taxid in tqdm(taxids, desc=f'Organizing PN files for [{len(taxids)}] Tax IDs.'):
        smp_list = [f'{smp_directory}/{source}' for source in hmm_pgap[hmm_pgap['taxonomic_range'] == taxid][
            'source_identifier']]
        with open(f'{db_directory}/{db_prefix}_{taxid}.pn', 'w') as f:
            f.write('\n'.join([f'{file}.smp' for file in set(smp_list)]))
    for taxid in taxids:
        pn2database(f'{db_directory}/{db_prefix}_{taxid}.pn')
        taxids_with_db.append(taxid)
    return taxids_with_db


def is_db_good(database, print_warning=True):
    for ext in ['aux', 'freq', 'loo', 'pdb', 'phr', 'pin', 'pos', 'pot', 'psq', 'ptf', 'pto', 'rps']:
        if not os.path.isfile(f'{database}.{ext}'):
            if print_warning:
                print(f'{database}.{ext} not found! Rebuilding database...')
            return False
    #print(f'{database} seems good!')
    return True


def cog2ec(
        cogblast, table=f'{sys.path[0]}/resources_directory/cog2ec.tsv',
        resources_dir=f'{sys.path[0]}/resources_directory'):
    if not os.path.isfile(table):
        run_command('python {0}/cog2ec.py -c {0}/eggnog4.protein_id_conversion.tsv -m {0}/NOG.members.tsv'.format(
            resources_dir), stdout=open(table, 'w'))
    return pd.merge(cogblast, pd.read_csv(table, sep='\t', names=['cog', 'EC number']), on='cog', how='left')


def cog2ko(cogblast, cog2ko_ssv=f'{sys.path[0]}/resources_directory/cog2ko.ssv'):
    if not os.path.isfile(cog2ko_ssv):
        directory = '/'.join(cog2ko_ssv.split('/')[:-1])
        web_locations = {
            'COG.mappings.v11.0.txt': 'https://stringdb-static.org/download/COG.mappings.v11.0.txt.gz',
            'protein.info.v11.0.txt': 'https://stringdb-static.org/download/protein.info.v11.0.txt.gz'}
        for file, link in web_locations.items():
            if not os.path.isfile(f'{directory}/{file}'):
                run_command(f'wget -P {directory} {link} -q')
                run_command(f'gunzip {directory}/{file}.gz')
        run_pipe_command(
            f"grep -E 'K[0-9]{5}$' {directory}/protein.info.v11.0.txt | awk '{{if (length($NF) == 6) print $1, $NF}}'",
            file=f'{directory}/string2ko.tsv')
        run_pipe_command(
            """awk '{{if (length($4) == 7) print $1"\t"$4}}' {0}/COG.mappings.v11.0.txt | sort | 
            join - {0}/string2ko.tsv""".format(directory), file=f'{directory}/cog2ko.ssv')
        df = pd.read_csv(f'{directory}/cog2ko.ssv', sep=' ', names=['StringDB', 'COG', 'KO'])
        df[['COG', 'KO']].groupby('COG')['KO'].agg([('KO', ','.join)]).reset_index().to_csv(
            f'{directory}/cog2ko.tsv', sep='\t', index=False, header=['cog', 'KO'])
    return pd.merge(cogblast, pd.read_csv(cog2ko_ssv, sep='\t'), on='cog', how='left')


def write_table(table, output, out_format='excel', header=True):
    if out_format == 'excel':
        table.to_excel(f'{output}.xlsx', index=False, header=header)
    elif out_format == 'tsv':
        table.to_csv(f'{output}.tsv', index=False, sep='\t', header=header)


def multi_sheet_excel(writer, data, sheet_name='Sheet', max_lines=1000000, index=False):
    if len(data) < max_lines:
        data.to_excel(writer, sheet_name=sheet_name, index=index)
    else:
        j = 1
        for i in range(0, len(data), max_lines):
            data.iloc[i:(i + max_lines)].to_excel(writer, sheet_name=f'{sheet_name} ({j})', index=index)
            j += 1
    return writer


def create_krona_plot(tsv, output=None, print_command=False):
    if output is None:
        output = tsv.replace('.tsv', '.html')
    run_command(f'ktImportText {tsv} -o {output}', print_command=print_command)


def write_cog_categories(data, output_basename):
    # COG categories quantification
    data = data.groupby(
        ['COG general functional category', 'COG functional category', 'Protein description', 'DB ID']
    ).size().reset_index().rename(columns={0: 'count'})
    data[['count'] + data.columns.tolist()[:-1]].to_csv(
        f'{output_basename}_quantification.tsv', sep='\t', index=False, header=None)
    create_krona_plot(f'{output_basename}_quantification.tsv', f'{output_basename}_quantification.html')


def count_on_file(expression, file, compressed=False):
    return int(check_output(f"{'zgrep' if compressed else 'grep'} -c '{expression}' {file}", shell=True))


def parse_fasta_on_memory(file):
    lines = [line.rstrip('\n') for line in open(file)]
    i = 0
    result = dict()
    while i < len(lines):
        if lines[i].startswith('>'):
            name = lines[i][1:].split()[0]
            result[name] = ''
            i += 1
            while i < len(lines) and not lines[i].startswith('>'):
                result[name] += lines[i]
                i += 1
    return pd.DataFrame.from_dict(result, orient='index', columns=['sequence'])


def write_fasta(data, output, protein_id_col):
    data[protein_id_col] = data[protein_id_col].apply(lambda x: f'>{x}')
    data.to_csv(output, sep='\n', header=False, index=False)


def split_fasta_by_taxid(file, tax_file, protein_id_col, tax_col, output):
    fastas = parse_fasta_on_memory(file)
    fastas.reset_index(inplace=True)
    tax_file = tax_file.reset_index().groupby(protein_id_col)[tax_col].first()
    tax_file = tax_file.reset_index()
    fastas = pd.merge(fastas, tax_file[[protein_id_col, tax_col]], left_on='index', right_on=protein_id_col, how='left')
    cols = fastas.columns.tolist()
    for col in [protein_id_col, 'index']:
        cols.remove(col)
    fastas = fastas.groupby(protein_id_col)[cols].first()
    for taxid in tqdm(set(tax_file[tax_col].tolist()), desc=f'Splitting sequences by taxon'):
        write_fasta(
            fastas[fastas[tax_col] == taxid].reset_index()[[protein_id_col, 'sequence']],
            f'{output}/tmp/{taxid}.fasta', protein_id_col)


def check_tax_databases(smp_directory, db_directory, db_prefix, taxids, hmm_pgap):
    taxids_lacking_db = []
    taxids_with_db = []
    for taxid in list(taxids):
        if not is_db_good(f'{db_directory}/{db_prefix}_{taxid}'):
            taxids_lacking_db.append(taxid)
        else:
            taxids_with_db.append(taxid)
    create_tax_db(smp_directory, db_directory, db_prefix, taxids_lacking_db, hmm_pgap)
    return taxids_with_db + taxids_lacking_db


def get_members_df(resources_directory):
    if os.path.isfile(f'{resources_directory}/members_df.tsv'):
        return pd.read_csv(f'{resources_directory}/members_df.tsv', sep='\t', index_col=0)
    members = pd.read_csv(f'{resources_directory}/NOG.members.tsv', sep='\t', header=None)
    members = members[members[1].str.startswith('COG')]
    members[5] = members[5].apply(lambda x: [name.split('.')[0] for name in x.split(',')])
    members_dict = {}
    for i in tqdm(range(len(members)), desc='Organizing COGs corresponding to each tax ID'):
        for taxid in members.iloc[i, 5]:
            if taxid in members_dict.keys():
                members_dict[taxid] += f',{members.iloc[i, 1]}'
            else:
                members_dict[taxid] = members.iloc[i, 1]
    members_df = pd.DataFrame.from_dict(members_dict, orient='index')
    members_df.columns = ['cogs']
    members_df.to_csv(f'{resources_directory}/members_df.tsv', sep='\t')
    return members_df


def check_cog_tax_database(smp_directory, db_directory):
    smps = glob(f'{smp_directory}/COG*.smp')
    for smp in tqdm(smps, desc=f'Checking split COG database for [{len(smps)}] COGs.'):
        name = smp.split('/')[-1].split('.')[0]
        with open(f'{db_directory}/{name}.pn', 'w') as f:
            f.write(smp)
        if not is_db_good(f'{db_directory}/{name}', print_warning=False):
            pn2database(f'{db_directory}/{name}.pn')


def cog_taxonomic_workflow(
        output, resources_directory, threads, tax_file, tax_col, members_df, max_target_seqs=1, evalue=1e-5):
    check_cog_tax_database(f'{resources_directory}/smps', f'{resources_directory}/dbs')  # for proteins with no taxonomy
    members_df.index = members_df.index.astype(str)
    members_taxids = members_df.index.tolist()
    db_report = pd.DataFrame(columns=['qseqid', 'sseqid', 'SUPERFAMILIES', 'SITES', 'MOTIFS'])
    for taxid in set(tax_file[tax_col].tolist()):
        # Run RPS-BLAST
        if taxid not in members_taxids:
            with Pool(processes=threads) as p:
                p.starmap(run_rpsblast, [(
                    f'{output}/tmp/tmp_{taxid}_{i}.fasta', f'{output}/asn/COG_{taxid}_{i}_aligned.asn',
                    f'{resources_directory}/dbs/COG', '1', max_target_seqs, evalue) for i in range(threads)
                    if os.path.isfile(f'{output}/tmp/tmp_{taxid}_{i}.fasta')])
        else:
            with Pool(processes=threads) as p:
                p.starmap(run_rpsblast, [(
                    f'{output}/tmp/tmp_{taxid}_{i}.fasta', f'{output}/asn/COG_{taxid}_{i}_aligned.asn',
                    ' '.join([f'{resources_directory}/{cog}' for cog in members_df.loc[taxid]['cogs']]), '1',
                    max_target_seqs, evalue) for i in range(threads)
                    if os.path.isfile(f'{output}/tmp/tmp_{taxid}_{i}.fasta')])
        # Convert ASN-11 to TAB-6
        with Pool(processes=threads) as p:
            p.starmap(run_blast_formatter, [(
                f'{output}/asn/COG_{taxid}_{i}_aligned.asn',
                f'{output}/blast/COG_{taxid}_{i}_aligned.blast') for i in range(threads)
                if os.path.isfile(f'{output}/asn/COG_{taxid}_{i}_aligned.asn')])
        # Convert ASN to RPSBPROC
        with Pool(processes=threads) as p:
            p.starmap(run_rpsbproc, [(
                f'{output}/asn/COG_{taxid}_{i}_aligned.asn', resources_directory, evalue) for i in range(threads)
                if os.path.isfile(f'{output}/asn/COG_{taxid}_{i}_aligned.asn')])
        for i in range(threads):
            if os.path.isfile(f'{output}/rpsbproc/COG_{taxid}_{i}_aligned.rpsbproc'):
                rpsbproc_report = get_rpsbproc_info(f'{output}/rpsbproc/COG_{taxid}_{i}_aligned.rpsbproc')
                if len(rpsbproc_report) > 0:
                    db_report = db_report.append(rpsbproc_report)
    db_report.to_csv(f'{output}/COG_report.tsv', sep='\t')


def check_regular_database(smp_directory, db_directory, db_prefix):
    if not is_db_good(f'{db_directory}/{db_prefix}'):
        print(f'Some part of {db_prefix} was not valid! Will rebuild!')
        smp_list = glob(f'{smp_directory}/{db_prefix}*.smp')
        with open(f'{db_directory}/{db_prefix}.pn', 'w') as f:
            f.write('\n'.join(smp_list))
        pn2database(f'{db_directory}/{db_prefix}.pn')
    else:
        print(f'A valid {db_prefix} split database was found!')


def load_relational_tables(resources_directory, tax_file=None):
    timed_message('Loading relational tables')
    cddid = parse_cddid(f'{resources_directory}/cddid_all.tbl')
    cddid['CDD ID'] = cddid['CDD ID'].apply(lambda x: f'CDD:{x}')
    hmm_pgap = pd.read_csv(f'{resources_directory}/hmm_PGAP.tsv', sep='\t', usecols=[1, 10, 12, 14, 15])
    hmm_pgap['source_identifier'] = [ide.split('.')[0] for ide in hmm_pgap['source_identifier']]
    hmm_pgap['source_identifier'] = hmm_pgap['source_identifier'].str.replace('PF', 'pfam')
    smps = [filename.split('/')[-1].rstrip('.smp') for filename in glob(f'{resources_directory}/smps/*.smp')]
    hmm_pgap = hmm_pgap[hmm_pgap['source_identifier'].isin(smps)]
    hmm_pgap['taxonomic_range'] = hmm_pgap['taxonomic_range'].fillna(0.0).apply(
        lambda x: str(int(x)) if type(x) == float else x)
    fun = pd.read_csv(f'{sys.path[0]}/fun.tsv', sep='\t')
    if tax_file is None:
        return cddid, hmm_pgap, fun, None, None
    taxonomy_df = pd.read_csv(f'{resources_directory}/taxonomy.tsv', sep='\t', index_col='taxid',
                              dtype={'taxid': str, 'name': str, 'rank': str, 'parent_taxid': str})
    taxonomy_df['parent_taxid'] = taxonomy_df['parent_taxid'].fillna('0').apply(lambda x: x.split('.')[0])
    members_df = get_members_df(resources_directory)
    members_df['cogs'] = members_df['cogs'].apply(lambda x: set(x.split(',')))
    return cddid, hmm_pgap, fun, taxonomy_df, members_df


def replace_spaces_with_commas(file, tmp_dir):
    timed_message('Replacing spaces for commas')
    run_pipe_command(f"sed -e 's/ /_/g' {file} > {tmp_dir}/tmp.fasta")
    return f'{tmp_dir}/tmp.fasta'


def split_fasta_by_threads(file, output_basename, threads):
    fasta = parse_fasta_on_memory(file)
    keys = list(split(fasta.index, threads))
    for i in range(threads):
        with open(f'{output_basename}_{i}.fasta', 'w') as f:
            for key in keys[i]:
                f.write(f'>{key}\n{fasta.loc[key, "sequence"]}\n')


def taxids_of_interest(tax_file, protein_id_col, tax_col, tax_df):
    tax_file = pd.read_csv(tax_file, sep='\t', index_col=protein_id_col, low_memory=False)
    tax_file[tax_col] = tax_file[tax_col].fillna(0.0).astype(int).astype(str)
    lineages, all_taxids = get_lineages_multiprocessing(set(tax_file[tax_col].tolist()), tax_df)
    return tax_file, lineages, all_taxids


def get_hmm_pgap_taxids(all_taxids, db_prefix, hmm_pgap):
    hmm_pgap = hmm_pgap[hmm_pgap['source_identifier'].str.startswith(db_prefix)]
    hmm_ids = set(hmm_pgap['taxonomic_range'])
    all_taxids_in_hmm_pgap = [tid for tid in all_taxids if tid in hmm_ids]  # each of these parents should have a database if it is possible to have it
    return all_taxids_in_hmm_pgap


def add_sequences(file, report):
    fasta = parse_fasta_on_memory(file)
    return pd.merge(report, fasta, left_on='qseqid', right_index=True, how='left')


def run_rpsbproc(asn_report, resources_directory, evalue):
    run_pipe_command(
        f'rpsbproc -i {asn_report} -o {asn_report.replace("asn", "rpsbproc")} -d {resources_directory} -e {evalue} '
        f'-m rep -f -t both 2>verbose.log')


def parse_rpsbproc_section(handler, line, section_name, i):
    data = list()
    if line.startswith(section_name):
        line = next(handler)
        while not line.startswith(f'END{section_name}'):
            data.append(line.rstrip('\n').split('\t')[i])
            line = next(handler)
        line = next(handler)
    return list(set(data)), line


def parse_rpsbproc(file):
    file = open(file)
    result = []
    try:
        line = [next(file) for i in range(3)][-1]
    except StopIteration:
        return result
    while line.startswith('#'):  # skip first section
        line = next(file)
    line = next(file, None)
    if line is None:
        return result
    while not line.startswith('ENDDATA'):
        line = next(file)
        while not line.startswith('ENDSESSION'):
            query = line.rstrip('\n').split('\t')[4]
            domains, superfamilies, sites, motifs = [], [], [], []
            line = next(file)
            while not line.startswith('ENDQUERY'):
                domains, line = parse_rpsbproc_section(file, line, 'DOMAINS', 3)
                superfamilies, line = parse_rpsbproc_section(file, line, 'SUPERFAMILIES', 3)
                sites, line = parse_rpsbproc_section(file, line, 'SITES', 7)
                motifs, line = parse_rpsbproc_section(file, line, 'MOTIFS', 5)
            result.append([query, domains, superfamilies, sites, motifs])
            line = next(file)
        line = next(file)
    result = pd.DataFrame(result, columns=['qseqid', 'sseqid', 'Superfamilies', 'Sites', 'Motifs'])
    return result


def run_blast_formatter(archive, output, outfmt='6'):
    run_pipe_command(f'blast_formatter -archive {archive} -outfmt {outfmt} -out {output} 2>verbose.log')


def get_rpsbproc_info(rpsbproc_report):
    if not os.path.isfile(rpsbproc_report):
        return pd.DataFrame(columns=['qseqid', 'sseqid', 'SUPERFAMILIES', 'SITES', 'MOTIFS'])
    rpsbproc_report = parse_rpsbproc(rpsbproc_report)
    if len(rpsbproc_report) > 0:
        for col in rpsbproc_report.columns.tolist()[2:]:  # exclude 'qseqid' and 'sseqid'
            rpsbproc_report[col] = rpsbproc_report[col].apply(','.join)
        rpsbproc_report = expand_by_list_column(rpsbproc_report, column='sseqid')
        rpsbproc_report.index = rpsbproc_report.index.astype(str)
        rpsbproc_report.sseqid = rpsbproc_report.sseqid.apply(lambda x: f'CDD:{x}')
        return rpsbproc_report
    else:
        return pd.DataFrame(columns=['qseqid', 'sseqid', 'SUPERFAMILIES', 'SITES', 'MOTIFS'])


def get_cdd_ec(description):
    m = re.compile("EC:([1-6\-].[0-9\-]+.[0-9\-]+.[0-9\-]+)\)").search(description)
    if m is None:
        return np.nan
    return m.group(1)


def add_db_info(report, db, resources_directory, output, hmm_pgap, fun):
    if db in ['CDD', 'Pfam', 'NCBIfam', 'Protein_Clusters', 'TIGRFAM']:
        report = pd.merge(report, hmm_pgap, left_on='DB ID', right_on='source_identifier', how='left')
        report.columns = report.columns.tolist()[:-4] + [
            'Protein description', 'EC number', 'taxonomic_range', 'taxonomic_range_name']
        if db == 'CDD':
            report['EC number'] = report['DB description'].apply(get_cdd_ec)
    elif db == 'Smart':
        smart_table = pd.read_csv(
            f'{resources_directory}/descriptions.pl', sep='\t', skiprows=2, header=None, usecols=[1, 2])
        smart_table.columns = ['Smart ID', 'Smart description']
        smart_table['Smart ID'] = smart_table['Smart ID'].str.replace('SM', 'smart')
        report = pd.merge(report, smart_table, left_on='DB ID', right_on='Smart ID', how='left')
        report.columns = report.columns.tolist()[:-1] + ['Protein description']
    elif db == 'KOG':
        kog_table = parse_kog(f'{resources_directory}/kog')
        kog_table = pd.merge(kog_table, fun, left_on='KOG functional category (letter)',
                             right_on='COG functional category (letter)', how='left')
        report = pd.merge(report, kog_table, left_on='DB ID', right_on='kog', how='left')
        report.columns = report.columns.tolist()[:-4] + ['Protein description'] + report.columns.tolist()[-3:]
        write_cog_categories(report, f'{output}/KOG')
    elif db == 'COG':
        cog_table = parse_whog(f'{resources_directory}/cog-20.def.tab')
        cog_table = pd.merge(cog_table, fun, on='COG functional category (letter)', how='left')
        report = pd.merge(report, cog_table, left_on='DB ID', right_on='cog', how='left')
        # cog2ec
        report = cog2ec(report, table=f'{resources_directory}/cog2ec.tsv', resources_dir=resources_directory)
        # cog2ko
        report = cog2ko(report, cog2ko_ssv=f'{resources_directory}/cog2ko.tsv')
        report.columns = report.columns.tolist()[:-5] + ['Protein description'] + report.columns.tolist()[-4:]
        report.to_csv(f'{output}/COG_report.tsv', sep='\t', index=False)
        write_cog_categories(report, f'{output}/COG')
    else:
        exit('Invalid database for retrieving further information!')
    return report


def custom_database_workflow(file, output, threads, max_target_seqs, evalue, database):
    dbs = database.split(',')
    for db in dbs:
        if not is_db_good(db):
            exit('Some inputted database was not valid!')
    # run annotation with rps-blast and inputted database
    timed_message('Running annotation with RPS-BLAST and inputted database as reference.')
    run_rpsblast(
        file, f'{output}/aligned.blast', ' '.join(dbs), threads=threads, max_target_seqs=max_target_seqs, evalue=evalue)


def clean_intermediates(output, base):
    for report in ['asn', 'blast', 'rpsbproc']:
        run_pipe_command(f'rm {output}/{report}/{base}_*_aligned.{report}')


def taxonomic_workflow(
        output, resources_directory, threads, lineages, all_taxids, databases_prefixes, base, hmm_pgap,
        max_target_seqs=1, evalue=1e-5):
    all_taxids += ['131567', '0']  # cellular organisms and no taxonomy
    hmm_pgap_taxids = get_hmm_pgap_taxids(all_taxids, databases_prefixes[base], hmm_pgap)
    taxids_with_db = check_tax_databases(
        f'{resources_directory}/smps', f'{resources_directory}/dbs', databases_prefixes[base], hmm_pgap_taxids,
        hmm_pgap)
    check_regular_database(f'{resources_directory}/smps', f'{resources_directory}/dbs', databases_prefixes[base])  # for proteins with no taxonomy
    dbs = {taxid: [
        f'{resources_directory}/dbs/{databases_prefixes[base]}_{parent_taxid}' for parent_taxid in
        lineages[taxid] + ['0'] if parent_taxid in taxids_with_db] for taxid in lineages.keys()}
    dbs = {**dbs, **{'0': [f'{resources_directory}/dbs/{databases_prefixes[base]}']}}       # no taxonomy is annotated with all
    db_report = pd.DataFrame(columns=['qseqid', 'sseqid', 'SUPERFAMILIES', 'SITES', 'MOTIFS'])
    for taxid in list(lineages.keys()) + ['0']:
        if os.path.isfile(f'{output}/tmp/{taxid}.fasta'):
            # Run RPS-BLAST
            with Pool(processes=threads) as p:
                p.starmap(run_rpsblast, [(
                    f'{output}/tmp/tmp_{taxid}_{i}.fasta', f'{output}/asn/{base}_{taxid}_{i}_aligned.asn',
                    ' '.join(dbs[taxid]), '1', max_target_seqs, evalue) for i in range(threads)
                    if os.path.isfile(f'{output}/tmp/tmp_{taxid}_{i}.fasta')])
            # Convert ASN-11 to TAB-6
            with Pool(processes=threads) as p:
                p.starmap(run_blast_formatter, [(
                    f'{output}/asn/{base}_{taxid}_{i}_aligned.asn',
                    f'{output}/blast/{base}_{taxid}_{i}_aligned.blast') for i in range(threads)
                    if os.path.isfile(f'{output}/asn/{base}_{taxid}_{i}_aligned.asn')])
            # Convert ASN to RPSBPROC
            with Pool(processes=threads) as p:
                p.starmap(run_rpsbproc, [(
                    f'{output}/asn/{base}_{taxid}_{i}_aligned.asn', resources_directory, evalue) for i in range(threads)
                    if os.path.isfile(f'{output}/asn/{base}_{taxid}_{i}_aligned.asn')])
            for i in range(threads):
                if os.path.isfile(f'{output}/rpsbproc/{base}_{taxid}_{i}_aligned.rpsbproc'):
                    rpsbproc_report = get_rpsbproc_info(f'{output}/rpsbproc/{base}_{taxid}_{i}_aligned.rpsbproc')
                    if len(rpsbproc_report) > 0:
                        db_report = db_report.append(rpsbproc_report)
    db_report.to_csv(f'{output}/{base}_report.tsv', sep='\t')


def multiprocess_workflow(
        output, resources_directory, threads, databases_prefixes, base, max_target_seqs=5, evalue=0.01):
    check_regular_database(f'{resources_directory}/smps', f'{resources_directory}/dbs', databases_prefixes[base])
    # Run RPS-BLAST
    with Pool(processes=threads) as p:
        p.starmap(run_rpsblast, [(
            f'{output}/tmp/tmp_{i}.fasta', f'{output}/asn/{base}_{i}_aligned.asn',
            f'{resources_directory}/dbs/{databases_prefixes[base]}', '1',
            max_target_seqs, evalue) for i in range(threads)])
    # Convert ASN-11 to TAB-6
    with Pool(processes=threads) as p:
        p.starmap(run_blast_formatter, [(
            f'{output}/asn/{base}_{i}_aligned.asn',
            f'{output}/blast/{base}_{i}_aligned.blast') for i in range(threads)
            if os.path.isfile(f'{output}/asn/{base}_{i}_aligned.asn')])
    run_pipe_command(f'cat {output}/blast/{base}_*_aligned.blast', file=f'{output}/blast/{base}_aligned.blast')
    # Convert ASN to RPSBPROC
    with Pool(processes=threads) as p:
        p.starmap(run_rpsbproc, [(
            f'{output}/asn/{base}_{i}_aligned.asn', resources_directory, evalue) for i in range(threads)
            if os.path.isfile(f'{output}/asn/{base}_{i}_aligned.asn')])
    db_report = pd.DataFrame(columns=['qseqid', 'sseqid', 'SUPERFAMILIES', 'SITES', 'MOTIFS'])
    for i in range(threads):
        if os.path.isfile(f'{output}/rpsbproc/{base}_{i}_aligned.rpsbproc'):
            rpsbproc_report = get_rpsbproc_info(f'{output}/rpsbproc/{base}_{i}_aligned.rpsbproc')
            if len(rpsbproc_report) > 0:
                db_report = db_report.append(rpsbproc_report)
    db_report.to_csv(f'{output}/{base}_report.tsv', sep='\t')


def organize_results(file, output, resources_directory, databases, hmm_pgap, cddid, fun, no_output_sequences=False):
    timed_message("Organizing annotation results")
    i = 1
    xlsx_report = pd.ExcelWriter(f'{output}/reCOGnizer_results.xlsx', engine='xlsxwriter')
    all_reports = pd.DataFrame()
    for db in databases:
        run_pipe_command(f'cat {output}/blast/{db}_*_aligned.blast', file=f'{output}/blast/{db}_aligned.blast')
        print(f'[{i}/{len(databases)}] Handling {db} annotation')
        report = pd.read_csv(f'{output}/{db}_report.tsv', sep='\t', index_col=0)
        report = pd.merge(
            report, parse_blast(f'{output}/blast/{db}_aligned.blast'), on=['qseqid', 'sseqid'], how='left')
        report = pd.merge(report, cddid, left_on='sseqid', right_on='CDD ID', how='left')
        del report['CDD ID']
        if not no_output_sequences:
            report = add_sequences(file, report)        # adding protein sequences if requested
        report = add_db_info(report, db, resources_directory, output, hmm_pgap, fun)
        report = report.drop_duplicates()
        # report = report[report['pident'] > pident]  # filter matches by pident - seems no longer implementable after rpsbproc integration
        all_reports = pd.concat([all_reports, report])
        multi_sheet_excel(xlsx_report, report, sheet_name=db)
        i += 1
        #clean_intermediates(output, db)
    all_reports.sort_values(by=['qseqid', 'DB ID']).to_csv(f'{output}/reCOGnizer_results.tsv', sep='\t', index=False)
    xlsx_report.save()


def read_ecmap(fh):
    enzymes = []
    proteins = []
    for line in fh:
        items = line.split("\t")
        m = re.compile("EC:[1-6\-].[0-9\-]+.[0-9\-]+.[0-9\-]+").search(items[2])
        try:
            ec = m.group().split(":")[1]
        except AttributeError:
            continue
        member = f"{items[0]}.{items[1]}"
        proteins.append(member)
        enzymes.append(ec)
    return enzymes, proteins


def ecmap(ec_file):
    with open(ec_file) as handler:
        enzymes, proteins = read_ecmap(handler)
    return enzymes, proteins


def read_cogmap(cogmap_handler):
    cogs = []
    proteins = []
    for line in cogmap_handler:
        items = line.split("\t")
        prots = items[-1].split(",")
        cog = [items[1]] * len(prots)
        cogs += cog
        proteins += prots
    return cogs, proteins


def cogmap(file):
    with open(file) as handler:
        cogs, proteins = read_cogmap(handler)
    return cogs, proteins


def determine_cog2ec(map_df, frac=0.5):
    # Group by cog and enzyme to get number of each EC assignment per cog
    map_df_counts = map_df.groupby(["enzyme", "cog"]).count().reset_index()
    map_df_counts.index = map_df_counts.cog
    map_df_counts.drop("cog", axis=1, inplace=True)
    map_df_counts.sort_index(inplace=True)
    # Count total number of proteins per cog
    cog_counts = map_df_counts.groupby(level=0).sum()
    # Divide enzyme assignment number by total protein number to get fraction of each assignment
    ecfrac = map_df_counts.protein.div(cog_counts.protein).reset_index()
    # Get index of where fraction is above threshold
    index = ecfrac.loc[ecfrac.protein >= frac].index
    # Return mappings where fraction is above threshold
    return map_df_counts.iloc[index]


def generate_cog2ec_df(conversion, members, resources_directory):
    enzymes, proteins = ecmap(conversion)
    ecmap_df = pd.DataFrame(data={"enzyme": enzymes, "protein": proteins})
    cogs, proteins = cogmap(members)
    cogmap_df = pd.DataFrame(data={"cog": cogs, "protein": proteins})
    map_df = pd.merge(ecmap_df, cogmap_df, left_on="protein", right_on="protein")
    cog2ec_df = determine_cog2ec(map_df)
    cog2ec_df.loc[:, "enzyme"].to_csv(f'{resources_directory}/cog2ec.tsv', sep="\t")


def main():
    args = get_arguments()

    if not os.path.isfile(f'{args.resources_directory}/recognizer_dwnl.timestamp') and not args.download_resources:
        args.download_resources = str2bool(input(
            'Resources seem to not have been downloaded for reCOGnizer yet. Do you want to download them? [Y/N] '))

    if args.download_resources:
        download_resources(
            args.resources_directory, quiet=args.quiet, skip_downloaded=args.skip_downloaded)

    if not hasattr(args, "file"):
        exit('No input file specified. Exiting.')

    cddid, hmm_pgap, fun, taxonomy_df, members_df = load_relational_tables(
        args.resources_directory, tax_file=args.tax_file)

    if not args.keep_spaces:
        args.file = replace_spaces_with_commas(args.file, f'{args.output}/tmp')     # if alters input file, internally alters args.file

    if args.database:  # if user database was inputted
        custom_database_workflow(
            args.file, args.output, args.threads, args.max_target_seqs, args.evalue, database=args.database)
    else:
        if args.tax_file is not None:
            tax_file, lineages, all_taxids = taxids_of_interest(
                args.tax_file, args.protein_id_col, args.tax_col, taxonomy_df)
            split_fasta_by_taxid(args.file, tax_file, args.protein_id_col, args.tax_col, args.output)
            # split FASTA for multiprocessing
            for taxid in tqdm(lineages.keys(), desc=timed_message('Splitting FASTA')):
                if os.path.isfile(f'{args.output}/tmp/{taxid}.fasta'):
                    split_fasta_by_threads(
                        f'{args.output}/tmp/{taxid}.fasta', f'{args.output}/tmp/tmp_{taxid}', args.threads)
        split_fasta_by_threads(args.file, f'{args.output}/tmp/tmp', args.threads)     # will likely always need to do this splitting
        databases_prefixes = {
            'CDD': 'cd', 'Pfam': 'pfam', 'NCBIfam': 'NF', 'Protein_Clusters': 'PRK', 'Smart': 'smart',
            'TIGRFAM': 'TIGR', 'COG': 'COG', 'KOG': 'KOG'}
        for base in args.databases:
            db_hmm_pgap = hmm_pgap[hmm_pgap['source_identifier'].str.startswith(databases_prefixes[base])]
            timed_message(f'Running annotation with RPS-BLAST and {base} database as reference.')
            if args.tax_file is not None and base in ['Pfam', 'NCBIfam', 'Protein_Clusters', 'TIGRFAM']:
                taxonomic_workflow(
                    args.output, args.resources_directory, args.threads, lineages, all_taxids, databases_prefixes,
                    base, db_hmm_pgap, max_target_seqs=args.max_target_seqs, evalue=args.evalue)
            elif base in ['COG'] and args.tax_file is not None and args.species_taxids:
                cog_taxonomic_workflow(
                    args.output, args.resources_directory, args.threads, tax_file, args.tax_col, members_df,
                    max_target_seqs=args.max_target_seqs, evalue=args.evalue)
            else:
                multiprocess_workflow(
                    args.output, args.resources_directory, args.threads, databases_prefixes, base,
                    max_target_seqs=args.max_target_seqs, evalue=args.evalue)

        organize_results(
            args.file, args.output, args.resources_directory, args.databases, hmm_pgap, cddid, fun,
            no_output_sequences=args.no_output_sequences)

        for directory in [f'{args.output}/{folder}' for folder in ['fasta', 'asn', 'blast', 'rpsbproc', 'tmp']]:
            files = glob(f'{directory}/*')
            for file in files:
                os.remove(file)


if __name__ == '__main__':
    start_time = time()
    main()
    timed_message(f'reCOGnizer analysis finished in {human_time(time() - start_time)}')
