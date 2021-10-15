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
from multiprocessing import Pool, cpu_count
from time import time, gmtime, strftime
from Bio import Entrez, SeqIO
from requests import get as requests_get
import xml.etree.ElementTree as ET

__version__ = '1.5.2'

Entrez.email = "A.N.Other@example.com"


def get_arguments():
    parser = ArgumentParser(
        description="reCOGnizer - a tool for domain based annotation with the COG database",
        epilog="Input file must be specified.")
    parser.add_argument("-t", "--threads", type=int, default=cpu_count() - 2,
                        help="Number of threads for reCOGnizer to use [max available - 2]")
    parser.add_argument("--evalue", type=float, default=1e-2, help="Maximum e-value to report annotations for [1e-2]")
    parser.add_argument("--pident", type=float, default=0, help="Minimum pident to report annotations for [0]")
    parser.add_argument(
        "-o", "--output", type=str, help="Output directory [reCOGnizer_results]", default='reCOGnizer_results')
    parser.add_argument(
        "-dr", "--download-resources", default=False, action="store_true",
        help='If resources for reCOGnizer are not available at "resources_directory" [false]')
    parser.add_argument(
        "-rd", "--resources-directory", type=str, default=os.path.expanduser('~/recognizer_resources'),
        help="Output directory for storing databases and other resources [~/recognizer_resources]")
    parser.add_argument(
        "-dbs", "--databases", type=str, nargs='+',
        choices=["CDD", "Pfam", "NCBIfam", "Protein_Clusters", "Smart", "TIGRFAM", "COG", "KOG"],
        default=["CDD", "Pfam", "NCBIfam", "Protein_Clusters", "Smart", "TIGRFAM", "COG", "KOG"],
        help="Databases to include in functional annotation [all available]")
    parser.add_argument(
        "-db", "--database", type=str,
        help="Basename of database for annotation. If multiple databases, use comma separated list (db1,db2,db3)")
    parser.add_argument(
        "--custom-database", action="store_true", default=False, help="If database was NOT produced by reCOGnizer")
    parser.add_argument(
        "-mts", "--max-target-seqs", default=1, help="Number of maximum identifications for each protein [1]")
    parser.add_argument(
        "--remove-spaces", action="store_true", default=False,
        help="BLAST ignores sequences IDs after the first space. "
             "This option changes all spaces to underscores to keep the full IDs.")
    parser.add_argument(
        "--no-output-sequences", action="store_true", default=False,
        help="Protein sequences from the FASTA input will be stored in their own column.")
    parser.add_argument(
        "--no-blast-info", action="store_true", default=False,
        help="Information from the alignment will be stored in their own columns.")
    parser.add_argument('-v', '--version', action='version', version=f'reCOGnizer {__version__}')

    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument(
        "-f", "--file", type=str, required=True, help="Fasta file with protein sequences for annotation")

    taxArguments = parser.add_argument_group('taxonomy arguments')
    taxArguments.add_argument(
        "--tax-file", help="File with taxonomic identification of proteins inputted (TSV). "
                           "Must have one line per query, query name on first column, taxid on second.")
    taxArguments.add_argument(
        "--protein-id-col", default='qseqid',
        help="Name of column with protein headers as in supplied FASTA file [qseqid]")
    taxArguments.add_argument(
        "--tax-col", default='Taxonomic identifier (GENUS)',
        help="Name of column with tax IDs of proteins [Taxonomic identifier (GENUS)]")

    args = parser.parse_args()

    args.output = args.output.rstrip('/')
    args.resources_directory = args.resources_directory.rstrip('/')

    for directory in [args.output, args.resources_directory]:
        if not os.path.isdir(directory):
            Path(directory).mkdir(parents=True, exist_ok=True)
            print(f'Created {directory}')

    return args


def timed_message(message):
    print(f'{strftime("%Y-%m-%d %H:%M:%S", gmtime())}: {message}')


def run_command(bash_command, print_command=True, stdout=None, stderr=None):
    if print_command:
        print(bash_command)
    run(bash_command.split(), stdout=stdout, stderr=stderr)


def run_pipe_command(bash_command, file='', mode='w', print_message=True):
    if print_message:
        print(bash_command)
    if file == '':
        Popen(bash_command, stdin=PIPE, shell=True).communicate()
    elif file == 'PIPE':
        return Popen(bash_command, stdin=PIPE, shell=True, stdout=PIPE).communicate()[0].decode('utf8')
    else:
        with open(file, mode) as output_file:
            Popen(bash_command, stdin=PIPE, shell=True, stdout=output_file).communicate()


def parse_fasta(filename):
    return SeqIO.parse(filename, "fasta")


def get_tabular_taxonomy(output):
    res = requests_get('https://ftp.uniprot.org/pub/databases/uniprot/current_release/rdf/taxonomy.rdf.xz')
    with open('taxonomy.rdf.xz', 'wb') as f:
        f.write(res.content)
    run_command(f'unxz taxonomy.rdf.xz')
    print('Reading RDF taxonomy')
    root = ET.parse('taxonomy.rdf').getroot()
    elems = root.findall('{http://www.w3.org/1999/02/22-rdf-syntax-ns#}Description')
    with open(output, 'w') as f:
        written = f.write('\t'.join(['taxid','name','rank','parent_taxid']) + '\n')
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


def download_resources(directory):
    for location in [
        # Download CDD
        'ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/cdd.tar.gz',
        'https://ftp.ncbi.nlm.nih.gov/pub/mmdb/cdd/cddid_all.tbl.gz',
        'https://ftp.ncbi.nlm.nih.gov/pub/mmdb/cdd/cdd.info',   # only for versions
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
            if str2bool(input(f"{directory}/{location.split('/')[-1]} exists. Overwrite? [Y/N] ")):
                run_command(f'wget {location} -P {directory}')

    for file in ['cddid_all.tbl', 'eggnog4.protein_id_conversion.tsv', 'NOG.members.tsv']:
        run_command(f'gunzip {directory}/{file}.gz')

    # Extract the smps
    if sys.platform == "darwin":
        if which('gtar') is None:
            run_command('brew install gnu-tar')
        tool = 'gtar'
    else:
        tool = 'tar'
    wd = os.getcwd()
    os.chdir(directory)
    run_pipe_command(f'{tool} -xzf cdd.tar.gz --wildcards "*.smp"')
    os.chdir(wd)

    get_tabular_taxonomy(f'{directory}/taxonomy.tsv')


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
    bashCommand = ['rpsblast', '-query', query, '-db', reference, '-out', output, '-outfmt', '11',
                   '-num_threads', str(threads), '-max_target_seqs', str(max_target_seqs), '-evalue', str(evalue)]
    print(' '.join(bashCommand))
    run(bashCommand, stdout=DEVNULL)


def parse_cddid(cddid):
    cddid = pd.read_csv(cddid, sep='\t', header=None)[[0, 1, 3]]
    cddid.columns = ['CDD ID', 'DB ID', 'DB description']
    #cddid['CDD ID'] = [f'CDD:{ide}' for ide in cddid['CDD ID']]
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


def create_split_db(smp_directory, output, db_prefix, threads=12):
    """
    Input:
        smp_directory: foldername where the SMP files are. These files are
        obtained from ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/cdd.tar.gz
        output: basename for PN and databases
        threads: STR, number of threads that the workflow will use
    Output:
        threads - 1 databases will be outputed, each with a consecutive part of
        the list of SMP files available. These databases are formated for RPS-BLAST
        search
    """
    print(f'Generating databases for [{threads}] threads.')
    smp_list = glob(f'{smp_directory}/{db_prefix}*.smp')
    parts = list(split(smp_list, threads))
    for i in range(len(parts)):
        with open(f'{output}/{db_prefix}_{threads}_{i}.pn', 'w') as f:
            f.write('\n'.join(parts[i]))
    pool = Pool()
    pool.map(pn2database, [f'{output}/{db_prefix}_{threads}_{i}.pn' for i in range(len(parts))])


def get_upper_taxids(taxid, tax_df):
    """
    :param taxid: str - taxID to get upper taxIDs from
    :param tax_df: pd.DataFrame - of read taxonomy.tsv (from taxonomy.rdf)
    :returns list - of upper taxIDs
    """
    if taxid == '0':
        return list
    taxids = list()
    while taxid != '1' and taxid != 'Taxon':
        taxids.append(taxid)
        taxid = tax_df.loc[taxid]['parent_taxid']
    return taxids


def get_all_taxids(taxids, taxonomy_df):
    all_taxids = list()
    for taxid in tqdm(taxids, desc=f'Listing all parent tax IDs for {len(taxids)} tax IDs'):
        all_taxids += get_upper_taxids(taxid, taxonomy_df)
    return set(all_taxids + ['0'])


def create_tax_db(smp_directory, output, db_prefix, taxids, hmm_pgap):
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
        with open(f'{output}/{db_prefix}_{taxid}.pn', 'w') as f:
            f.write('\n'.join([f'{file}.smp' for file in set(smp_list)]))
    '''
    pool = Pool()
    pool.map(pn2database, [f'{output}/{db_prefix}_{tax_id}.pn' for tax_id in taxids])
    '''
    for taxid in taxids:
        pn = f'{output}/{db_prefix}_{taxid}.pn'
        run_command(f"makeprofiledb -in {pn} -title {pn.split('.pn')[0]} -out {pn.split('.pn')[0]}")
        taxids_with_db.append(taxid)
    return taxids_with_db


def is_db_good(database):
    for ext in ['aux', 'freq', 'loo', 'pdb', 'phr', 'pin', 'pos', 'pot', 'psq', 'ptf', 'pto', 'rps']:
        if not os.path.isfile(f'{database}.{ext}'):
            print(f'{database}.{ext} not found!')
            return False
    print(f'{database} seems good!')
    return True


def cog2ec(cogblast, table=f'{sys.path[0]}/resources_directory/cog2ec.tsv',
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
        for file in web_locations.keys():
            if not os.path.isfile(f'{directory}/{file}'):
                run_command(f'wget -P {directory} {web_locations[file]}')
                run_command(f'gunzip {directory}/{file}.gz')
        run_pipe_command(
            """grep -E 'K[0-9]{5}$' """ + directory + """/protein.info.v11.0.txt | 
            awk '{{if (length($NF) == 6) print $1, $NF}}'""",
            file=f'{directory}/string2ko.tsv')
        run_pipe_command(
            """awk '{{if (length($4) == 7) print $1"\t"$4}}' {0}/COG.mappings.v11.0.txt | sort | 
            join - {0}/string2ko.tsv""".format(directory), file=f'{directory}/cog2ko.ssv')
        df = pd.read_csv(f'{directory}/cog2ko.ssv', sep=' ', names=['StringDB', 'COG', 'KO'])
        df[['COG', 'KO']].groupby('COG')['KO'].agg([('KO', ','.join)]).reset_index().to_csv(
            'f{directory}/cog2ko.tsv', sep='\t', index=False, header=['cog', 'KO'])
    return pd.merge(cogblast, pd.read_csv(cog2ko_ssv, sep='\t'), on='cog', how='left')


def write_table(table, output, out_format='excel', header=True):
    if out_format == 'excel':
        table.to_excel(f'{output}.xlsx', index=False, header=header)
    elif out_format == 'tsv':
        table.to_csv(f'{output}.tsv', index=False, sep='\t', header=header)


def multi_sheet_excel(writer, data, sheet_name='Sheet', lines=1000000, index=False):
    if len(data) < lines:
        data.to_excel(writer, sheet_name=sheet_name, index=index)
    else:
        for i in range(0, len(data), lines):
            j = min(i + lines, len(data))
            data.iloc[i:(i + lines)].to_excel(writer, sheet_name=f'{sheet_name} ({j})', index=index)
    return writer


def create_krona_plot(tsv, output=None):
    if output is None:
        output = tsv.replace('.tsv', '.html')
    run_command(f'ktImportText {tsv} -o {output}')


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
    print(f'Parsing {file}')
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
    tax_file.reset_index(inplace=True)
    fastas = parse_fasta_on_memory(file)
    fastas = pd.merge(fastas, tax_file[[protein_id_col, tax_col]], left_index=True, right_on=protein_id_col, how='left')
    cols = fastas.columns.tolist()
    cols.remove(protein_id_col)
    fastas = fastas.groupby(protein_id_col)[cols].first()
    for taxid in tqdm(set(tax_file[tax_col].tolist()), desc=f'Splitting sequences by taxon'):
        write_fasta(
            fastas[fastas[tax_col] == taxid].reset_index()[[protein_id_col, 'sequence']],
            f'{output}/tmp/{taxid}.fasta', protein_id_col)


def check_tax_databases(smp_directory, output, db_prefix, taxids, hmm_pgap):
    taxids_lacking_db = list()
    taxids_with_db = list()
    for taxid in list(taxids):
        if not is_db_good(f'{smp_directory}/{db_prefix}_{taxid}'):
            taxids_lacking_db.append(taxid)
        else:
            taxids_with_db.append(taxid)
    create_tax_db(smp_directory, output, db_prefix, taxids_lacking_db, hmm_pgap)
    return taxids_with_db + taxids_lacking_db


def check_threads_database(smp_directory, db_prefix, threads):
    dbs = [f'{smp_directory}/{db_prefix}_{threads}_{i}' for i in range(threads)]
    i = 0
    remake_db = False
    while not remake_db and i < len(dbs):
        if not is_db_good(dbs[i]):
            print(f'Some part of {db_prefix} was not valid!')
            remake_db = True
        i += 1
    if remake_db:
        create_split_db(smp_directory, smp_directory, db_prefix, threads)
    else:
        print(f'A valid {db_prefix} split database for [{threads}] threads was found!')


def load_relational_tables(resources_directory):
    timed_message('Loading relational tables')
    cddid = parse_cddid(f'{resources_directory}/cddid_all.tbl')
    cddid['CDD ID'] = cddid['CDD ID'].apply(lambda x: f'CDD:{x}')
    hmm_pgap = pd.read_csv(f'{resources_directory}/hmm_PGAP.tsv', sep='\t', usecols=[1, 10, 12, 14, 15])
    hmm_pgap['source_identifier'] = [ide.split('.')[0] for ide in hmm_pgap['source_identifier']]
    hmm_pgap['source_identifier'] = hmm_pgap['source_identifier'].str.replace('PF', 'pfam')
    smps = [filename.split('/')[-1].rstrip('.smp') for filename in glob(f'{resources_directory}/*.smp')]
    hmm_pgap = hmm_pgap[hmm_pgap['source_identifier'].isin(smps)]
    hmm_pgap['taxonomic_range'] = hmm_pgap['taxonomic_range'].fillna(0.0).apply(
        lambda x: str(int(x)) if type(x) == float else x)
    fun = pd.read_csv(f'{sys.path[0]}/fun.tsv', sep='\t')
    taxonomy_df = pd.read_csv(f'{resources_directory}/taxonomy.tsv', sep='\t', index_col='taxid',
                              dtype={'taxid': str, 'name': str, 'rank': str, 'parent_taxid': str})
    taxonomy_df['parent_taxid'] = taxonomy_df['parent_taxid'].fillna('0')
    return cddid, hmm_pgap, fun, taxonomy_df


def replace_spaces_with_commas(file):
    timed_message('Replacing spaces for commas')
    run_pipe_command(f"sed -i -e 's/ /_/g' {file}")


def taxids_of_interest(tax_file, protein_id_col, tax_col):
    tax_file = pd.read_csv(tax_file, sep='\t', index_col=protein_id_col, low_memory=False)
    tax_file[tax_col] = tax_file[tax_col].fillna(0.0).astype(int).astype(str)
    return tax_file, set(tax_file[tax_col].tolist())


def get_hmm_pgap_taxids(main_taxids, db_prefix, hmm_pgap, taxonomy_df):
    all_taxids = get_all_taxids(main_taxids, taxonomy_df)                   # of the supplied taxids, their parents are calculated
    all_taxids = list(all_taxids) + ['131567', '0']                         # cellular organisms and no taxonomy
    hmm_pgap = hmm_pgap[hmm_pgap['source_identifier'].str.startswith(db_prefix)]
    hmm_ids = set(hmm_pgap['taxonomic_range'])
    all_taxids_in_hmm_pgap = [tid for tid in all_taxids if tid in hmm_ids]  # each of these parents should have a database if it is possible to have it
    return all_taxids_in_hmm_pgap


def add_sequences(file, report):
    fasta = parse_fasta(file)
    report['Sequence'] = pd.Series(dtype=str)
    report.set_index('qseqid', inplace=True)
    for entry in fasta:
        report.at[entry.id, 'Sequence'] = str(entry.seq)
    report.reset_index(inplace=True)
    return report


def run_rpsbproc(asn_report, resources_directory, evalue):
    run_pipe_command(
        f'rpsbproc -i {asn_report} -o {asn_report.replace(".asn", ".better")} -d {resources_directory} -e {evalue} '
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
    next(file)
    next(file)
    line = next(file)
    result = pd.DataFrame(columns=['DOMAINS', 'SUPERFAMILIES', 'SITES', 'MOTIFS'])
    while line.startswith('#'):     # skip first section
        line = next(file)
    line = next(file)
    while not line.startswith('ENDDATA'):
        line = next(file)
        while not line.startswith('ENDSESSION'):
            query = line.rstrip('\n').split('\t')[4]
            domains, superfamilies, sites, motifs = list(), list(), list(), list()
            line = next(file)
            while not line.startswith('ENDQUERY'):
                domains, line = parse_rpsbproc_section(file, line, 'DOMAINS', 3)
                superfamilies, line = parse_rpsbproc_section(file, line, 'SUPERFAMILIES', 3)
                sites, line = parse_rpsbproc_section(file, line, 'SITES', 7)
                motifs, line = parse_rpsbproc_section(file, line, 'MOTIFS', 5)
            result = result.append(pd.DataFrame(
                [[domains, superfamilies, sites, motifs]], columns=result.columns, index=[query]))
            line = next(file)
        line = next(file)
    result.reset_index(inplace=True)
    result.columns = ['qseqid', 'sseqid'] + result.columns.tolist()[2:]
    return result


def run_blast_formatter(archive, output, outfmt='6'):
    run_pipe_command(f'blast_formatter -archive {archive} -outfmt {outfmt} -out {output} 2>verbose.log')


def get_post_processing(asn_report, resources_directory, evalue):
    run_rpsbproc(asn_report, resources_directory, evalue)
    rpsbproc_report = parse_rpsbproc(asn_report.replace(".asn", ".better"))
    if len(rpsbproc_report) > 0:
        for col in rpsbproc_report.columns.tolist()[2:]:    # exclude 'qseqid' and 'sseqid'
            rpsbproc_report[col] = rpsbproc_report[col].apply(','.join)
        rpsbproc_report = expand_by_list_column(rpsbproc_report, column='sseqid')
        rpsbproc_report.index = rpsbproc_report.index.astype(str)
        rpsbproc_report.sseqid = rpsbproc_report.sseqid.apply(lambda x: f'CDD:{x}')
        return rpsbproc_report
    else:
        return pd.DataFrame(columns=['qseqid', 'sseqid', 'SUPERFAMILIES', 'SITES', 'MOTIFS'])


def add_db_info(report, db, resources_directory, output, hmm_pgap, fun):
    if db in ['CDD', 'Pfam', 'NCBIfam', 'Protein_Clusters', 'TIGRFAM']:
        report = pd.merge(report, hmm_pgap, left_on='DB ID', right_on='source_identifier', how='left')
        report.columns = report.columns.tolist()[:-4] + [
            'Protein description', 'EC number', 'taxonomic_range', 'taxonomic_range_name']
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


def taxonomic_db_workflow(
        output, resources_directory, threads, main_taxids, databases_prefixes, base, hmm_pgap, taxonomy_df,
        max_target_seqs=1, evalue=1e-5):
    hmm_pgap_taxids = get_hmm_pgap_taxids(main_taxids, databases_prefixes[base], hmm_pgap, taxonomy_df)
    taxids_with_db = check_tax_databases(
        resources_directory, resources_directory, databases_prefixes[base], hmm_pgap_taxids, hmm_pgap)
    dbs = {taxid: [
        f'{resources_directory}/{databases_prefixes[base]}_{parent_taxid}' for parent_taxid in
        get_upper_taxids(taxid, taxonomy_df) + ['0'] if parent_taxid in taxids_with_db] for taxid in main_taxids}

    with Pool(processes=threads) as p:
        p.starmap(run_rpsblast, [(
            f'{output}/{taxid}.fasta', f'{output}/{base}_{taxid}_aligned.asn', ' '.join(dbs[taxid]), '1',
            max_target_seqs, evalue) for taxid in main_taxids])

    with Pool(processes=threads) as p:
        p.starmap(run_blast_formatter, [(
            f'{output}/{base}_{taxid}_aligned.asn', f'{output}/{base}_{taxid}_aligned.blast') for taxid in main_taxids])

    db_report = pd.DataFrame(columns=['qseqid', 'sseqid', 'SUPERFAMILIES', 'SITES', 'MOTIFS'])
    for taxid in main_taxids:
        report_6 = get_post_processing(f'{output}/{base}_{taxid}_aligned.asn', resources_directory, evalue)
        blast = parse_blast(f'{output}/{base}_{taxid}_aligned.blast')
        db_report = db_report.append(pd.merge(report_6, blast, on=['qseqid', 'sseqid'], how='left'))
    db_report.to_csv(f'{output}/{base}_report.tsv', sep='\t')


def multithreaded_db_workflow(
        file, output, resources_directory, threads, databases_prefixes, base, max_target_seqs=5, evalue=0.01):
    check_threads_database(resources_directory, databases_prefixes[base], threads)
    run_rpsblast(
        file, f'{output}/{base}_aligned.asn',
        ' '.join([f'{resources_directory}/{databases_prefixes[base]}_{threads}_{i}' for i in range(threads)]),
        threads=threads, max_target_seqs=max_target_seqs, evalue=evalue)
    run_blast_formatter(f'{output}/{base}_aligned.asn', f'{output}/{base}_aligned.blast')
    report_6 = get_post_processing(f'{output}/{base}_aligned.asn', resources_directory, evalue)
    pd.merge(report_6, parse_blast(f'{output}/{base}_aligned.blast'), on=['qseqid', 'sseqid'], how='left').to_csv(
        f'{output}/{base}_report.tsv', sep='\t')


def main():
    args = get_arguments()

    if args.download_resources:
        download_resources(args.resources_directory)

    cddid, hmm_pgap, fun, taxonomy_df = load_relational_tables(args.resources_directory)

    if args.remove_spaces:
        replace_spaces_with_commas(args.file)

    if args.database:  # if database was inputted
        custom_database_workflow(
            args.file, args.output, args.threads, args.max_target_seqs, args.evalue, database=args.database)
    else:
        if hasattr(args, "tax_file"):
            tax_file, main_taxids = taxids_of_interest(args.tax_file, args.protein_id_col, args.tax_col)
            split_fasta_by_taxid(args.file, tax_file, args.protein_id_col, args.tax_col, args.output)
        databases_prefixes = {
            'CDD': 'cd', 'Pfam': 'pfam', 'NCBIfam': 'NF', 'Protein_Clusters': 'PRK', 'Smart': 'smart',
            'TIGRFAM': 'TIGR', 'COG': 'COG', 'KOG': 'KOG'}

        for base in args.databases:
            db_hmm_pgap = hmm_pgap[hmm_pgap['source_identifier'].str.startswith(databases_prefixes[base])]
            timed_message(f'Running annotation with RPS-BLAST and {base} database as reference.')
            if hasattr(args, "tax_file") and base in ['Pfam', 'NCBIfam', 'Protein_Clusters', 'TIGRFAM']:
                taxonomic_db_workflow(
                    args.output, args.resources_directory, args.threads, main_taxids, databases_prefixes, base,
                    db_hmm_pgap, taxonomy_df, max_target_seqs=args.max_target_seqs, evalue=args.evalue)
            else:
                multithreaded_db_workflow(
                    args.file, args.output, args.resources_directory, args.threads, databases_prefixes, base,
                    max_target_seqs=args.max_target_seqs, evalue=args.evalue)

        timed_message("Organizing annotation results")
        i = 1
        xlsx_report = pd.ExcelWriter(f'{args.output}/reCOGnizer_results.xlsx', engine='xlsxwriter')
        all_reports = pd.DataFrame()
        for db in args.databases:
            print(f'[{i}/{len(args.databases)}] Handling {db} annotation')
            report = pd.read_csv(f'{args.output}/{db}_report.tsv', sep='\t', index_col=0)
            report = pd.merge(report, cddid, left_on='sseqid', right_on='CDD ID', how='left')
            del report['CDD ID']
            if not args.no_output_sequences:
                report = add_sequences(args.file, report)       # adding protein sequences if requested
            report = add_db_info(report, db, args.resources_directory, args.output, hmm_pgap, fun)
            report = report[report['pident'] > args.pident]     # filter matches by pident
            all_reports = pd.concat([all_reports, report])
            multi_sheet_excel(xlsx_report, report, sheet_name=db)
            i += 1
        all_reports.sort_values(by=['qseqid', 'DB ID']).to_csv(
            f'{args.output}/reCOGnizer_results.tsv', sep='\t', index=False)
        xlsx_report.save()


if __name__ == '__main__':
    start_time = time()
    main()
    print(f'reCOGnizer analysis finished in {strftime("%Hh%Mm%Ss", gmtime(time() - start_time))}')
