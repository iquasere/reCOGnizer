#!/usr/bin/env python
"""
reCOGnizer - a tool for functional annotation with COGs

By João Sequeira

Nov 2019
"""

import pandas as pd
from time import gmtime, strftime
import argparse, sys, os, multiprocessing, glob, subprocess, pathlib

__version__ = '1.2.2'

def get_arguments():    
    parser = argparse.ArgumentParser(description="reCOGnizer - a tool for domain based annotation with the COG database",
        epilog="Input file must be specified.")
    parser.add_argument("-t", "--threads", type = str, 
                        default = str(multiprocessing.cpu_count() - 2),
                        help = """Number of threads for reCOGnizer to use. 
                        Default is number of CPUs available minus 2.""")
    parser.add_argument("-o", "--output", type = str, help = "Output directory",
                        default = 'reCOGnizer_results'),
    parser.add_argument("-rd", "--resources-directory", type = str, 
                        help = "Output directory for storing COG databases and other resources",
                        default = sys.path[0] + '/Databases')
    parser.add_argument("-db", "--database", type = str,
                        help = """Basename of COG database for annotation. 
                        If multiple databases, use comma separated list (db1,db2,db3)""")
    parser.add_argument("--custom-database", action = "store_true", default = False,
                        help = "If database was NOT produced by reCOGnizer")
    parser.add_argument("-seqs", "--max-target-seqs", type = str,
                        help="""Number of maximum identifications for each protein.
                        Default is 1.""", default = "1")
    parser.add_argument("--tsv", action = "store_true", default = False,
                        help="Tables will be produced in TSV format (and not EXCEL).")
    parser.add_argument('-v', '--version', action='version', version='reCOGnizer ' + __version__)
    
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument("-f", "--file", type = str, required = True,
                        help="Fasta file with protein sequences for annotation")

    args = parser.parse_args()
    
    args.output = args.output.rstrip('/')
    args.resources_directory = args.resources_directory.rstrip('/')

    for directory in [args.output, args.resources_directory]:
        if not os.path.isdir(directory):
            pathlib.Path(directory).mkdir(parents=True, exist_ok=True)
            print('Created ' + directory)

    return args

'''
Input:
    message: a message to be printed
Output:
    will print the message with the time in human readable format
'''
def timed_message(message):
    print(strftime("%Y-%m-%d %H:%M:%S", gmtime()) + ': ' + message)

'''
Input:
    bashCommand: str - command to run
    print_command: boolean - if it should print the command
Output:
    command will be run... hopefully
'''
def run_command(bashCommand, print_command = True, stdout = None):
    if print_command:
        print(bashCommand)
    subprocess.run(bashCommand.split(), stdout = stdout)
    
'''
Input:
    database_directory: str - directory to store files for database construction
        and reCOGnizer workflow
Output:
    files to construct database will be downloaded to database_directory
'''
def download_resources(database_directory):
    if not os.path.isfile('{}/COG0001.smp'.format(database_directory)):
        print('{}/COG0001.smp not found! Retrieving from cdd.tar.gz...'.format(database_directory))
        if not os.path.isfile('{}/cdd.tar.gz'.format(database_directory)):
            print('{}/cdd.tar.gz not found! Downloading...'.format(database_directory))
            run_command('wget ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/cdd.tar.gz -P {}'.format(database_directory))
        wd = os.getcwd()
        os.chdir(database_directory)
        print('tar -xzf {}/cdd.tar.gz --wildcards "COG*.smp"'.format(database_directory))
        subprocess.Popen('tar -xzf cdd.tar.gz --wildcards "COG*.smp"', shell = True).communicate() # I couldn't, for the life of me, put the -C or --directory flags to work. No idea what happened, this just works
        os.chdir(wd)
    if not os.path.isfile('{}/cddid.tbl'.format(database_directory)):
        print('{}/cddid.tbl not found! Downloading...'.format(database_directory))
        run_command('wget ftp://ftp.ncbi.nlm.nih.gov/pub/mmdb/cdd/cddid.tbl.gz -P {}'.format(database_directory))
        run_command('gunzip {}/cddid.tbl.gz'.format(database_directory))
    if not os.path.isfile('{}/fun.txt'.format(database_directory)):
        print('{}/fun.txt not found! Downloading...'.format(database_directory))
        run_command('wget ftp://ftp.ncbi.nlm.nih.gov/pub/COG/COG/fun.txt -P {}'.format(database_directory))
    if not os.path.isfile('{}/whog'.format(database_directory)):
        print('{}/whog not found! Downloading...'.format(database_directory))
        run_command('wget ftp://ftp.ncbi.nlm.nih.gov/pub/COG/COG/whog -P {}'.format(database_directory))
    
'''
Input: 
    fasta: str - name of a fasta file of proteins to be annotated
    output: str - filename of rps_blast to be created
    cog: str - COG blast DB basename from ftp://ftp.ncbi.nlm.nih.gov/pub/mmdb/cdd/little_endian/Cog_LE.tar.gz
    threads: str - number of threads to use
    max_target_seqs: str - number of identifications to obtain for each protein
Output: 
    annotated file with CDD IDs
'''
def run_rpsblast(fasta, output, cog, threads = '0', max_target_seqs = '1'):
    bashCommand = ['rpsblast', '-query', fasta,  '-db', cog, '-out', output, 
                   '-outfmt', '6', '-num_threads', threads, '-max_target_seqs', 
                   max_target_seqs]
    print(' '.join(bashCommand))
    subprocess.run(bashCommand)
    
'''
Input:
    cddid: str - filename of cddid.tbl file
Output:
    returns pandas.DataFrame relating CDD ID to COG function
'''
def parse_cddid(cddid):
    cddid = pd.read_csv(cddid, sep = '\t', header = None)[[0,1]]
    cddid.columns = ['CDD ID', 'cog']
    cddid['CDD ID'] = ['CDD:{}'.format(str(ide)) for ide in cddid['CDD ID']]
    return cddid[cddid['cog'].str.startswith('COG')]

'''
Input: 
    fun: str - the fun.txt filename
Output: 
    returns result: dict - in the form {COG category (letter): (COG category (name), COG supercategory)}
'''      
def parse_fun(fun):
    lines = open(fun).readlines()
    result = dict()
    supercategory = lines[0]
    i = 1
    while i < len(lines):
        line = lines[i].rstrip('\n')
        if '[' in line:
            letter = line.split('[')[-1].split(']')[0]
            name = line.split('] ')[-1]
            result[letter] = [supercategory.rstrip('\n'), name]
        else:
            supercategory = line
        i += 1
    result = pd.DataFrame.from_dict(result).transpose().reset_index()
    result.columns = ['COG functional category (letter)', 
        'COG general functional category', 'COG functional category']
    return result

'''
Input:
    whog: str - filename of whog file
Output:
    returns df: pandas.DataFrame - COG, COG functional category (letter), COG protein description
'''
def parse_whog(whog):
    lines = list()
    for line in [line.rstrip('\n') for line in open(whog).readlines() if line.startswith('[')]:
        line = line.split()
        lines.append([line[0][1], line[1], ' '.join(line[2:])])
    df = pd.DataFrame(lines)
    df.columns = ['COG functional category (letter)', 'cog', 'COG protein description']
    return df

'''
Input:
    blast: str - filename of blast outputted by RPS-BLAST with CDD identifications
Output:
    pandas.DataFrame with CDD IDs converted to COGs and respective categories
'''
def cdd2cog(blast, cddid = sys.path[0] + '/Databases/cddid.tbl', 
            fun = sys.path[0] + '/Databases/fun.txt', 
            whog = sys.path[0] + '/Databases/whog'):
    result = pd.read_csv(blast, sep = '\t', header = None)[[0,1]]
    result.columns = ['qseqid', 'CDD ID']
    cddid = parse_cddid(cddid)
    fun = parse_fun(fun)
    whog = parse_whog(whog)
    whog = pd.merge(whog, fun, on = 'COG functional category (letter)',  how = 'left')
    result = pd.merge(result, cddid, on = 'CDD ID', how = 'left')
    result = pd.merge(result, whog, on = 'cog', how = 'left')
    return result[['qseqid', 'COG general functional category', 
                   'COG functional category', 'COG protein description', 'cog']]

'''
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
'''
def create_split_cog_db(smp_directory, output, threads = '12'):
    '''
    Input:
        a: list - list to be split
        n: int - number of parts into
    Output:
        list - the list split in n parts
    '''
    def split(a, n):
        k, m = divmod(len(a), n)
        return (a[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(n))
    
    database_reporter = '/'.join(output.split('/')[:-1]) + '/databases.txt'
    dbs = (open(database_reporter).read().split('\n') if
           os.path.isfile(database_reporter) else list())
    if threads in dbs:
        print('Already built COG database for [' + threads + '] threads.')
    else:
        print('Generating COG databases for [' + threads + '] threads.')
        smp_list = glob.glob(smp_directory + '/COG*.smp')
        parts = list(split(smp_list, int(threads)))
        for i in range(len(parts)):
            with open('{}_{}_{}.pn'.format(output, threads, str(i)), 'w') as f:
                f.write('\n'.join(parts[i]))

        for file in ['{}_{}_{}.pn'.format(output, threads, str(i)) for i in range(len(parts))]:
            run_command('makeprofiledb -in {0} -title {1} -out {1}'.format(     # -title and -out options are defaulted as input file name to -in argument; -dbtype default is 'rps'
                    file, file.split('.pn')[0]))
            
        with open(database_reporter,'w') as dr:
            dr.write('\n'.join(dbs + [threads]))
        
'''
Input: 
    database: str - database basename
Output:
    boolean - True if it seems a valid database, false otherwise
'''
def validate_database(database):
    for ext in ['aux','freq','loo','phr','pin','pn','psd','psi','psq','rps']:
        if not os.path.isfile('{}.{}'.format(database, ext)):
            return False
    return True

'''
Input:
    files_directory: str
Output:
    files will be downloaded if missing
'''
def download_eggnog_files(directory = sys.path[0] + '/Databases'):
    web_locations = {'cog2ec.py': 'https://bitbucket.org/scilifelab-lts/lts-workflows-sm-metagenomics/raw/screening_legacy/lts_workflows_sm_metagenomics/source/utils/cog2ec.py',
                     'eggnog4.protein_id_conversion.tsv': 'http://eggnogdb.embl.de/download/eggnog_4.5/eggnog4.protein_id_conversion.tsv.gz',
                     'NOG.members.tsv': 'http://eggnogdb.embl.de/download/eggnog_4.5/data/NOG/NOG.members.tsv.gz'}
    for file in ['cog2ec.py', 'eggnog4.protein_id_conversion.tsv', 'NOG.members.tsv']:
        if not os.path.isfile('{}/{}'.format(directory, file)):
            print('{}/{} not found! Downloading...'.format(directory, file))
            run_command('wget -P {} {}'.format(directory, web_locations[file]))
            if web_locations[file][-3:] == '.gz':
                run_command('gunzip {}/{}.gz'.format(directory, file))
        else:
            print('{}/{} found!'.format(directory, file))

'''
Input:
    cogblast: pandas.DataFrame - qseqid, cog categories
Output:
    pandas.DataFrame - same dataframe with added "EC number" column
'''
def cog2ec(cogblast, table = sys.path[0] + '/Databases/cog2ec.tsv', 
           resources_dir = sys.path[0] + '/Databases'):
    if not os.path.isfile(table):
        download_eggnog_files(directory = resources_dir)
        run_command('python {0}/cog2ec.py -c {0}/eggnog4.protein_id_conversion.tsv -m {0}/NOG.members.tsv'.format(
                resources_dir), stdout = open(table, 'w'))
    cog2ec = pd.read_csv(table, sep = '\t', names = ['cog', 'EC number'])
    return pd.merge(cogblast, cog2ec, on = 'cog', how = 'left')

'''
Input:
    table: pandas.DataFrame - table to write
    output: str - filename of output
    out_format: str - 'tsv' or 'excel' format to write as
Output:
    table will be written in the format specified
'''
def write_table(table, output, out_format = 'excel', header = True):
    if out_format == 'excel':
        table.to_excel(output + '.xlsx', index = False, header = header)
    elif out_format == 'tsv':
        table.to_csv(output + '.tsv', index = False, sep = '\t', header = header)

'''
Input:
    tsv: filename of TSV file to be inputed. Must have the format 
    value\tcategorie1\tcategorie2\t..., with no header
    output: filename of HTML krona plot to output
Output:
    A krona plot will be created at output if it has been specified, else
    at tsv.replace('.tsv','.html')
'''
def create_krona_plot(tsv, output = None):
    if output is None:
        output = tsv.replace('.tsv','.html')
    conda_exec = subprocess.check_output('which conda'.split()).decode('utf8')
    run_command('{}/bin/ktImportText {} -o {}'.format(
            conda_exec.split('/bin')[0], tsv, output))

def main():
    
    # get arguments
    args = get_arguments()
    
    if args.database:
        if not args.custom_database:                                            # if database was built by reCOGnizer
            args.threads = int(args.database.split('_')[-1])
            databases = ['{}_{}'.format(args.database, str(i)) for i in range(args.threads)]
        else:
            databases = args.database.split(',')
        for database in databases:
            if not validate_database(args.database):
                print('Database not valid!')
                exit()
    else:
        # check if necessary files exist to build database
        download_resources(args.resources_directory) 
        
        # create database if it doesn't exit
        timed_message('Checking if database exists for {} threads.'.format(args.threads))
        create_split_cog_db(args.resources_directory, 
                            args.resources_directory + '/COG', args.threads)
    
        # set database(s)
        databases = [pn.split('.pn')[0] for pn in glob.glob('{}/COG_{}_*.pn'.format(
                args.resources_directory, args.threads))]
            
    # run annotation with psi-blast and COG database
    timed_message('Running annotation with PSI-BLAST and COG database as reference.')
    run_rpsblast(args.file, args.output + '/cdd_aligned.blast', ' '.join(databases),
                 threads = args.threads, max_target_seqs = args.max_target_seqs)
    
    # convert CDD IDs to COGs
    timed_message('Converting CDD IDs to respective COG IDs.')
    cogblast = cdd2cog(args.output + '/cdd_aligned.blast',
                       cddid = args.resources_directory + '/cddid.tbl', 
                       fun = args.resources_directory + '/fun.txt', 
                       whog = args.resources_directory + '/whog')
    
    # cog2ec
    timed_message('Converting COG IDs to EC numbers.')
    cogblast = cog2ec(cogblast, table = args.resources_directory + '/cog2ec.tsv',
                      resources_dir = args.resources_directory)
    
    # write protein COG assignment
    out_format = 'tsv' if args.tsv else 'excel'
    write_table(cogblast,
                args.output + '/protein2cog', 
                out_format = out_format)
    
    timed_message('Protein ID to COG and EC number is available at {}.'.format(
            args.output + '/protein2cog'))
    
    # quantify COG categories
    timed_message('Quantifying COG categories.')
    cog_quantification = cogblast.groupby(['COG general functional category',
            'COG functional category', 'COG protein description', 'cog']).size(
            ).reset_index().rename(columns={0:'count'})
    write_table(cog_quantification, 
                args.output + '/cog_quantification', 
                out_format = out_format)
    timed_message('COG categories quantification is available at {}.'.format(
            args.output + '/cog_quantification'))
    
    # represent that quantification in krona plot
    timed_message('Creating Krona plot representation.')
    write_table(cog_quantification[['count'] + cog_quantification.columns.tolist()[:-1]], 
                         args.output + '/krona',
                         header = False,
                         out_format = 'tsv')
    create_krona_plot(args.output + '/krona.tsv', args.output + '/cog_quantification.html')
    
if __name__ == '__main__':
    main()