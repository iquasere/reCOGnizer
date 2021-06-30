#!/usr/bin/env python
"""
reCOGnizer - a tool for functional annotation with COGs

By João Sequeira

Nov 2019
"""

import argparse
import glob
import multiprocessing
import os
import pathlib
import shutil
import subprocess
import sys
import time
import numpy as np
import pandas as pd
from multiprocessing import Pool
from time import gmtime, strftime
from progressbar import ProgressBar

__version__ = '1.4.9'


def get_arguments():
    parser = argparse.ArgumentParser(
        description="reCOGnizer - a tool for domain based annotation with the COG database",
        epilog="Input file must be specified.")
    parser.add_argument("-t", "--threads", type=str,
                        default=str(multiprocessing.cpu_count() - 2),
                        help="Number of threads for reCOGnizer to use [max available - 2]")
    parser.add_argument("--evalue", type=float, default=1e-2, help="Maximum e-value to report annotations for [1e-2]")
    parser.add_argument("--pident", type=float, default=0, help="Minimum pident to report annotations for [0]")
    parser.add_argument("-o", "--output", type=str, help="Output directory [reCOGnizer_results]",
                        default='reCOGnizer_results')
    parser.add_argument("-dr", "--download-resources", default=False, action="store_true",
                        help='If resources for reCOGnizer are not available at "resources_directory" [false]')
    parser.add_argument("-rd", "--resources-directory", type=str,
                        help="Output directory for storing databases and other resources [~/resources_directory]",
                        default=os.path.expanduser('~/resources_directory'))
    parser.add_argument("-dbs", "--databases", type=str, nargs='+',
                        choices=["CDD", "Pfam", "NCBIfam", "Protein_Clusters", "Smart", "TIGRFAM", "COG", "KOG"],
                        default=["CDD", "Pfam", "NCBIfam", "Protein_Clusters", "Smart", "TIGRFAM", "COG", "KOG"],
                        help="Databases to include in functional annotation")
    parser.add_argument("-db", "--database", type=str, help="""Basename of database for annotation. 
                        If multiple databases, use comma separated list (db1,db2,db3)""")
    parser.add_argument("--custom-database", action="store_true", default=False,
                        help="If database was NOT produced by reCOGnizer")
    parser.add_argument("-seqs", "--max-target-seqs", type=str,
                        help="Number of maximum identifications for each protein [1]", default="1")
    parser.add_argument("--remove-spaces", action="store_true", default=False,
                        help="""BLAST ignores sequences IDs after the first space.
                        This option changes all spaces to underscores to keep the full IDs.""")
    parser.add_argument("--no-output-sequences", action="store_true", default=False,
                        help="Protein sequences from the FASTA input will be stored in their own column.")
    parser.add_argument("--no-blast-info", action="store_true", default=False,
                        help="Information from the alignment will be stored in their own columns.")
    parser.add_argument('-v', '--version', action='version', version=f'reCOGnizer {__version__}')

    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument("-f", "--file", type=str, required=True,
                               help="Fasta file with protein sequences for annotation")

    args = parser.parse_args()

    args.output = args.output.rstrip('/')
    args.resources_directory = args.resources_directory.rstrip('/')

    for directory in [args.output, args.resources_directory]:
        if not os.path.isdir(directory):
            pathlib.Path(directory).mkdir(parents=True, exist_ok=True)
            print(f'Created {directory}')

    return args


def timed_message(message):
    print(f'{strftime("%Y-%m-%d %H:%M:%S", gmtime())}:{message}')


def run_command(bashCommand, print_command=True, stdout=None):
    if print_command:
        print(bashCommand)
    subprocess.run(bashCommand.split(), stdout=stdout)


def run_pipe_command(bashCommand, file='', mode='w', print_message=True):
    if print_message:
        print(bashCommand)
    if file == '':
        subprocess.Popen(bashCommand, stdin=subprocess.PIPE, shell=True).communicate()
    elif file == 'PIPE':
        return subprocess.Popen(bashCommand, stdin=subprocess.PIPE, shell=True,
                                stdout=subprocess.PIPE).communicate()[0].decode('utf8')
    else:
        with open(file, mode) as output_file:
            subprocess.Popen(bashCommand, stdin=subprocess.PIPE, shell=True, stdout=output_file).communicate()


def parse_fasta(file):
    lines = [line.rstrip('\n') for line in open(file)]
    i = 0
    sequences = dict()
    while i < len(lines):
        if lines[i].startswith('>'):
            name = lines[i][1:]
            sequences[name] = ''
            i += 1
            while i < len(lines) and not lines[i].startswith('>'):
                sequences[name] += lines[i]
                i += 1
    return sequences


def download_resources(directory):
    for location in [
        # Download CDD
        'ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/cdd.tar.gz',
        'https://ftp.ncbi.nlm.nih.gov/pub/mmdb/cdd/cddid_all.tbl.gz',
        # COG categories
        'ftp.ncbi.nlm.nih.gov/pub/COG/COG2020/data/fun-20.tab',
        'ftp.ncbi.nlm.nih.gov/pub/COG/COG2020/data/cog-20.def.tab',
        # COG2EC
        'https://bitbucket.org/scilifelab-lts/lts-workflows-sm-metagenomics/raw/screening_legacy/lts_workflows_sm_metagenomics/source/utils/cog2ec.py',
        'http://eggnogdb.embl.de/download/eggnog_4.5/eggnog4.protein_id_conversion.tsv.gz',
        'http://eggnogdb.embl.de/download/eggnog_4.5/data/NOG/NOG.members.tsv.gz',
        # NCBIfam, TIGRFAM, Pfam, PRK (protein clusters)
        'https://ftp.ncbi.nlm.nih.gov/hmm/4.0/hmm_PGAP.tsv',
        # SMART
        'https://smart.embl.de/smart/descriptions.pl',
        # KOG
        'https://ftp.ncbi.nlm.nih.gov/pub/COG/KOG/kog'
    ]:
        run_command(f'wget {location} -P {directory}')

    for file in ['cddid_all.tbl', 'eggnog4.protein_id_conversion.tsv', 'NOG.members.tsv']:
        run_command(f'gunzip {directory}/{file}.gz')

    # Extract the smps
    if sys.platform == "darwin":
        if shutil.which('gtar') is None:
            run_command('brew install gnu-tar')
        tool = 'gtar'
    else:
        tool = 'tar'
    wd = os.getcwd()
    os.chdir(directory)
    run_pipe_command(f'{tool} -xzf cdd.tar.gz --wildcards "*.smp"')
    os.chdir(wd)


def run_rpsblast(query, output, reference, threads='0', max_target_seqs='1', evalue=10e-2):
    # This run_command is different because of reference, which can't be split by space
    bashCommand = ['rpsblast', '-query', query, '-db', reference, '-out', output, '-outfmt', '6',
                   '-num_threads', threads, '-max_target_seqs', max_target_seqs, '-evalue', str(evalue)]
    print(' '.join(bashCommand))
    subprocess.run(bashCommand)


def parse_cddid(cddid):
    cddid = pd.read_csv(cddid, sep='\t', header=None)[[0, 1, 3]]
    cddid.columns = ['CDD ID', 'DB ID', 'DB description']
    cddid['CDD ID'] = [f'CDD:{ide}' for ide in cddid['CDD ID']]
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
    blast = pd.read_csv(file, sep='\t', header=None)
    blast.columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
                     'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
    return blast


def cdd2id(blast, cddid=f'{sys.path[0]}/resources_directory/cddid_all.tbl'):
    """
    Input:
        blast: str - filename of blast outputted by RPS-BLAST with CDD identifications
    Output:
        pandas.DataFrame with CDD IDs converted to COGs and respective categories
    """
    blast = parse_blast(blast)
    cddid = parse_cddid(cddid)
    return pd.merge(blast, cddid, left_on='sseqid', right_on='CDD ID', how='left')


def pn2database(pn):
    run_command(f"makeprofiledb -in {pn} -title {pn.split('.pn')[0]} -out {pn.split('.pn')[0]}")


def split(a, n):
    k, m = divmod(len(a), n)
    return (a[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(n))


def create_split_db(smp_directory, output, db_prefix, threads='12'):
    """
    Input:
        smp_directory: foldername where the SMP files are. These files are
        obtained from ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/cdd.tar.gz
        output: basename for PN and databases
        threads: STR, number of threads that the workflow will use
        step: number of SMP files per database
    Output:
        threads - 1 databases will be outputed, each with a consecutive part of
        the list of SMP files available. These databases are formated for RPS-BLAST
        search
    """
    print(f'Generating databases for [{threads}] threads.')
    smp_list = glob.glob('{}/{}*.smp'.format(smp_directory, db_prefix))
    parts = list(split(smp_list, int(threads)))
    for i in range(len(parts)):
        with open(f'{output}/{db_prefix}_{threads}_{i}.pn', 'w') as f:
            f.write('\n'.join(parts[i]))

    pool = Pool()
    pool.map(pn2database, [f'{output}/{db_prefix}_{threads}_{i}.pn' for i in range(len(parts))])


def validate_database(database):
    for ext in ['aux', 'freq', 'loo', 'phr', 'pin', 'pn', 'psd', 'psi', 'psq', 'rps']:
        if not os.path.isfile(f'{database}.{ext}'):
            print(f'{database}.{ext} not found!')
            return False
    return True


def cog2ec(cogblast, table=f'{sys.path[0]}/resources_directory/cog2ec.tsv',
           resources_dir=f'{sys.path[0]}/resources_directory'):
    if not os.path.isfile(table):
        run_command('python {0}/cog2ec.py -c {0}/eggnog4.protein_id_conversion.tsv -m {0}/NOG.members.tsv'.format(
            resources_dir), stdout=open(table, 'w'))
    return pd.merge(cogblast, pd.read_csv(table, sep='\t', names=['cog', 'EC number']), on='cog', how='left')


def cog2ko(cogblast, cog2ko=f'{sys.path[0]}/resources_directory/cog2ko.ssv'):
    if not os.path.isfile(cog2ko):
        directory = '/'.join(cog2ko.split('/')[:-1])
        web_locations = {
            'COG.mappings.v11.0.txt': 'https://stringdb-static.org/download/COG.mappings.v11.0.txt.gz',
            'protein.info.v11.0.txt': 'https://stringdb-static.org/download/protein.info.v11.0.txt.gz'}
        for file in web_locations.keys():
            if not os.path.isfile(f'{directory}/{file}'):
                run_command(f'wget -P {directory} {web_locations[file]}')
                run_command(f'gunzip {directory}/{file}.gz')
        run_pipe_command(
            """grep -E 'K[0-9]{5}$' """ + directory + """/protein.info.v11.0.txt | awk '{{if (length($NF) == 6) print $1, $NF}}'""",
            file=f'{directory}/string2ko.tsv')
        run_pipe_command(
            """awk '{{if (length($4) == 7) print $1"\t"$4}}' {0}/COG.mappings.v11.0.txt | sort | join - {0}/string2ko.tsv""".format(
                directory),
            file=f'{directory}/cog2ko.ssv')
        df = pd.read_csv(f'{directory}/cog2ko.ssv', sep=' ', names=['StringDB', 'COG', 'KO'])
        df[['COG', 'KO']].groupby('COG')['KO'].agg([('KO', ','.join)]).reset_index().to_csv(
            'f{directory}/cog2ko.tsv', sep='\t', index=False, header=['cog', 'KO'])
    return pd.merge(cogblast, pd.read_csv(cog2ko, sep='\t'), on='cog', how='left')


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


def write_cog_categories(data, output_basename, db='COG'):
    # COG categories quantification
    data = data.groupby(
        ['COG general functional category', 'COG functional category', 'Protein description', 'DB ID']
    ).size().reset_index().rename(columns={0: 'count'})
    data.to_excel(f'{output_basename}_quantification.xlsx')
    data[['count'] + data.columns.tolist()[:-1]].to_csv(
        f'{output_basename}_quantification.tsv', sep='\t', index=False, header=None)
    create_krona_plot(f'{output_basename}_quantification.tsv', f'{output_basename}_quantification.html')


def main():
    # get arguments
    args = get_arguments()

    if args.download_resources:
        download_resources(args.resources_directory)

    if args.database:  # if database was inputted
        inputted_db = True
        database_groups = args.database.split(',')
        for database in database_groups:
            if not validate_database(database):
                exit('Some inputted database was not valid!')
    else:
        inputted_db = False

        database_groups = list()
        databases_prefixes = {'CDD': 'cd',
                              'Pfam': 'pfam',
                              'NCBIfam': 'NF',
                              'Protein_Clusters': 'PRK',
                              'Smart': 'smart',
                              'TIGRFAM': 'TIGR',
                              'COG': 'COG',
                              'KOG': 'KOG'}

        for base in args.databases:
            database_group = [f'{args.resources_directory}/{databases_prefixes[base]}_{args.threads}_{i}'
                              for i in range(int(args.threads))]
            remake_dbs = False

            i = 0
            while not remake_dbs and i < len(database_group):
                if not validate_database(database_group[i]):
                    print(f'Some part of {base} was not valid!')
                    remake_dbs = True
                i += 1

            if remake_dbs:
                # create database if it doesn't exit
                create_split_db(args.resources_directory, args.resources_directory,
                                databases_prefixes[base], args.threads)
            else:
                print(f'A valid {base} split database for [{args.threads}] threads was found!')

            database_groups.append((base, database_group))

    # Replacing spaces for commas
    timed_message('Replacing spaces for commas')
    if args.remove_spaces:
        run_pipe_command(f"sed -i -e 's/ /_/g' {args.file}")

    if inputted_db:
        # run annotation with rps-blast and database
        timed_message('Running annotation with RPS-BLAST and inputted database as reference.')
        run_rpsblast(args.file, f'{args.output}/aligned.blast', ' '.join(database_groups),
                     threads=args.threads, max_target_seqs=args.max_target_seqs, evalue=args.evalue)
    else:
        for db_group in database_groups:
            # run annotation with rps-blast and database
            timed_message(f'Running annotation with RPS-BLAST and {db_group[0]} database as reference.')
            run_rpsblast(args.file, f'{args.output}/{db_group[0]}_aligned.blast', ' '.join(db_group[1]),
                         threads=args.threads, max_target_seqs=args.max_target_seqs, evalue=args.evalue)

    if inputted_db:
        exit()

    # Load the relational tables
    cddid = parse_cddid(f'{args.resources_directory}/cddid_all.tbl')
    ncbi_table = pd.read_csv(f'{args.resources_directory}/hmm_PGAP.tsv', sep='\t', usecols=[1, 10, 12, 15])
    ncbi_table['source_identifier'] = [ide.split('.')[0] for ide in ncbi_table['source_identifier']]
    ncbi_table['source_identifier'] = ncbi_table['source_identifier'].str.replace('PF', 'pfam')
    fun = pd.read_csv(f'{sys.path[0]}/fun.tsv', sep='\t')

    timed_message("Organizing annotation results")

    i = 1
    j = len(args.databases)
    for db in args.databases:
        print(f'[{i}/{j}] Handling {db} identifications')
        report = parse_blast(f'{args.output}/{db}_aligned.blast')
        report = pd.merge(report, cddid, left_on='sseqid', right_on='CDD ID', how='left')
        del report['CDD ID']

        # adding protein sequences if requested
        if not args.no_output_sequences:
            fasta = parse_fasta(args.file)
            fasta = pd.DataFrame.from_dict(fasta, orient='index')
            fasta.columns = ['Sequence']
            report = pd.merge(report, fasta, left_on='qseqid', right_index=True, how='left')

        if db in ['CDD', 'Pfam', 'NCBIfam', 'Protein_Clusters', 'TIGRFAM']:
            report = pd.merge(report, ncbi_table, left_on='DB ID', right_on='source_identifier', how='left')
            report.columns = report.columns.tolist()[:-3] + ['Protein description', 'EC number', 'taxonomic_range_name']

        elif db == 'Smart':
            smart_table = pd.read_csv(f'{args.resources_directory}/descriptions.pl',
                                      sep='\t', skiprows=2, header=None, usecols=[1, 2])
            smart_table.columns = ['Smart ID', 'Smart description']
            smart_table['Smart ID'] = smart_table['Smart ID'].str.replace('SM', 'smart')
            report = pd.merge(report, smart_table, left_on='DB ID', right_on='Smart ID', how='left')
            report.columns = report.columns.tolist()[:-1] + ['Protein description']

        elif db == 'KOG':
            kog_table = parse_kog(f'{args.resources_directory}/kog')
            kog_table = pd.merge(kog_table, fun, left_on='KOG functional category (letter)',
                                 right_on='COG functional category (letter)', how='left')
            report = pd.merge(report, kog_table, left_on='DB ID', right_on='kog', how='left')
            report.columns = report.columns.tolist()[:-4] + ['Protein description'] + report.columns.tolist()[-3:]

            write_cog_categories(report, f'{args.output}/KOG', db='KOG')

        else:
            cog_table = parse_whog(f'{args.resources_directory}/cog-20.def.tab')
            cog_table = pd.merge(cog_table, fun, on='COG functional category (letter)', how='left')
            report = pd.merge(report, cog_table, left_on='DB ID', right_on='cog', how='left')
            # cog2ec
            report = cog2ec(report, table=f'{args.resources_directory}/cog2ec.tsv',
                            resources_dir=args.resources_directory)
            # cog2ko
            report = cog2ko(report, cog2ko=f'{args.resources_directory}/cog2ko.tsv')
            report.columns = report.columns.tolist()[:-5] + ['Protein description'] + report.columns.tolist()[-4:]
            report.to_csv(f'{args.output}/COG_report.tsv', sep='\t', index=False)

            write_cog_categories(report, f'{args.output}/COG', db='COG')

        report.to_csv(f'{args.output}/{db}_report.tsv', sep='\t', index=False)
        i += 1

    timed_message(f'Filtering matches for: pident < {args.pident}')
    timed_message(f'Organizing all results at {args.output}/reCOGnizer_results.xlsx')
    writer = pd.ExcelWriter(f'{args.output}/reCOGnizer_results.xlsx', engine='xlsxwriter')

    pbar = ProgressBar()
    for base in pbar(args.databases):
        report = pd.read_csv(f'{args.output}/{base}_report.tsv', sep='\t')
        report = report[report['pident'] > args.pident]
        multi_sheet_excel(writer, report, sheet_name=base)
    writer.save()

    timed_message(f'Organizing all results at {args.output}/reCOGnizer_results.tsv')
    all = pd.DataFrame()
    for base in args.databases:
        report = pd.read_csv(f'{args.output}/{base}_report.tsv', sep='\t')
        report = report[report['pident'] > args.pident]
        all = pd.concat([all, report])
    all.sort_values(by=['qseqid', 'DB ID']).to_csv(f'{args.output}/reCOGnizer_results.tsv', sep='\t', index=False)


if __name__ == '__main__':
    start_time = time.time()
    main()
    print(f'reCOGnizer analysis finished in {strftime("%Hh%Mm%Ss", gmtime(time.time() - start_time))}')
