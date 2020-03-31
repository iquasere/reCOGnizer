# -*- coding: utf-8 -*-
"""
reCOGnizer - a tool for functional annotation with COGs

By João Sequeira

Nov 2019
"""

import pandas as pd
from time import gmtime, strftime
import argparse, shutil, sys, os, multiprocessing, glob, subprocess, pathlib

def get_arguments():    
    parser = argparse.ArgumentParser(description="reCOGnizer - a tool for domain based annotation with the COG database",
        epilog="Input file must be specified.")
    parser.add_argument("-f", "--file", type = str, required = True,
                        help="Fasta file with protein sequences for annotation")
    parser.add_argument("-t", "--threads", type = str, 
                        default = str(multiprocessing.cpu_count() - 2),
                        help = """Number of threads for reCOGnizer to use. 
                        Default is number of CPUs available minus 2.""")
    parser.add_argument("-o", "--output", type = str, help = "Output directory",
                        default = 'reCOGnizer_results'),
    parser.add_argument("-odb", "--output-databases", type = str, 
                        help = "Output directory for storing COG databases",
                        default = sys.path[0] + '/Databases')
    parser.add_argument("-db", "--database", type = str,
                        help = """Basename of COG database for annotation. 
                        If multiple databases, use comma separated list (db1,db2,db3)""")
    parser.add_argument("--database-by-recognizer", type = str,
                        help = "If inputed database was produced by reCOGnizer")
    parser.add_argument("-seqs", "--max-target-seqs", type = str,
                        help="""Number of maximum identifications for each protein.
                        Default is 1.""", default = "1")
    parser.add_argument("--tsv", action = "store_true", default = False,
                        help="Tables will be produced in TSV format (and not EXCEL).")

    args = parser.parse_args()
    
    args.output = args.output.rstrip('/')
    args.output_databases = args.output_databases.rstrip('/')

    for directory in [args.output, args.output_databases]:
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
def run_command(bashCommand, print_command = True):
    if print_command:
        print(bashCommand)
    subprocess.run(bashCommand.split())
    
'''
Input:
Output:
'''
def download_resources(database_directory):
    run_command('bash {}/download_resources.bash {}'.format(sys.path[0], database_directory))
    
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
    blast: str - name of blast output with CDD IDs
    output: str - output foldername where to store the resuls folder
    cddid: str - name of cddid summary file from ftp://ftp.ncbi.nlm.nih.gov/pub/mmdb/cdd/cddid.tbl.gz
    fun: str - name of fun file available at ftp://ftp.ncbi.nlm.nih.gov/pub/COG/COG/fun.txt
    whog: str - name of whog file available at ftp://ftp.ncbi.nlm.nih.gov/pub/COG/COG/whog
Output: 
    CDD IDs will be converted to respective COGs, with results stored in 'output'
    folder
'''
def annotate_cogs(blast, output, cddid, fun, whog):
    if os.path.isdir(os.getcwd() + '/results'):                                 # the cdd2cog tool does not overwrite, and fails if results directory already exists
        print('Eliminating ' + os.getcwd() + '/results')
        shutil.rmtree(os.getcwd() + '/results', ignore_errors=True)             # is not necessary when running the tool once, but better safe then sorry!
    run_command('perl reCOGnizer/cdd2cog.pl -r {} -c {} -f {} -w {}'.format(
            blast, cddid, fun, whog))
    if os.path.isdir(output + '/results'):
        shutil.rmtree(output + '/results')
    shutil.copytree('results', output + '/results')

'''
Input: 
    name of cddblast to parse
Output: 
    pandas.DataFrame object
'''      
def parse_cogblast(cogblast):
    cogblast = pd.read_csv(cogblast, header=None, skiprows = 1, sep = '\t', low_memory=False)
    cogblast = cogblast[list(range(0,14))+[18]]                             #several extra columns are produced because of bad formatting
    cogblast.columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 
               'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'cog',
               'functional categories', 'COG protein description']
    return cogblast
    
'''
Input: 
    cogblast: the output from cdd2go, a blast file with CDD and COG annotations
    fun: the fun.txt file available at ftp://ftp.ncbi.nih.gov/pub/COG/COG/fun.txt
Output: 
    returns pandas.DataFrame with the functional categories intrisic levels 
    reorganized into corresponding columns
'''      
def organize_cdd_blast(cogblast, fun = sys.path[0] + '/Databases/fun.txt'):
    cogblast = parse_cogblast(cogblast)
    cogblast = cogblast[cogblast['functional categories'].notnull()]
    cog_relation = parse_fun(fun)
    data = [cog_relation[functional_category] for functional_category in cogblast['functional categories']]
    result = pd.DataFrame(data)
    result.columns = ['COG general functional category','COG functional category']
    result = pd.concat([result[['COG general functional category','COG functional category']], 
                        cogblast[['COG protein description','cog','qseqid']]], axis = 1)
    return result

'''
Input: the fun.txt file available at ftp://ftp.ncbi.nih.gov/pub/COG/COG/fun.txt
Output: a dictionary in the form {COG category (letter): (COG category (name), COG supercategory)}
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
    return result

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
def create_split_cog_db(smp_directory, output, threads = '6'):
    '''
    Input:
        a: list - list to be splited
        n: int - number of parts into
    Output:
        list - a splited in n parts
        output + /databases.txt will have added the threads for new database
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
intput: 
    database: str - database basename
output:
    boolean - True if it seems a valid database, false otherwise
'''
def validate_database(database):
    for ext in ['aux','freq','loo','phr','pin','pn','psd','psi','psq','rps']:
        if not os.path.isfile('{}.{}'.format(database, ext)):
            return False
    return True

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
    run_command('perl Krona/KronaTools/scripts/ImportText.pl {} -o {}'.format(tsv, output))
        
def main():
    
    # get arguments
    args = get_arguments()
    
    if args.database:
        if args.database_by_recognizer:                                         # if database was built by reCOGnizer
            args.threads = int(args.database.split('_')[-1])
            databases = [args.database + '_' + i for i in range(args.threads)]
        else:
            databases = args.database.split(',')
        for database in databases:
            if not validate_database(args.database):
                print('Database not valid!')
                exit()
    else:
        # check if necessary files exist to build database
        for file in ['cddid.tbl', 'fun.txt', 'whog']:
            if not os.path.isfile('{}/Databases/{}'.format(sys.path[0], file)):
                print('{}/Databases/{} not found!'.format(sys.path[0], file))
                download_resources('{}/Databases'.format(sys.path[0]))
        
        # create database if it doesn't exit
        timed_message('Checking if database exists for {} threads.'.format(args.threads))
        create_split_cog_db('{}/Databases'.format(sys.path[0]), args.output_databases + '/COG', args.threads)
    
        # set database(s)
        databases = [pn.split('.pn')[0] for pn in glob.glob('{}/COG_{}_*.pn'.format(
                args.output_databases, args.threads))]
            
    # run annotation with psi-blast and COG database
    timed_message('Running annotation with PSI-BLAST and COG database as reference.')
    run_rpsblast(args.file, args.output + '/cdd_aligned.blast', ' '.join(databases),
                 threads = args.threads, max_target_seqs = args.max_target_seqs)
    
    # convert CDD IDs to COGs
    timed_message('Converting CDD IDs to respective COG IDs.')
    annotate_cogs(args.output + '/cdd_aligned.blast', args.output, 
        sys.path[0] + '/Databases/cddid.tbl', sys.path[0] + '/Databases/fun.txt', 
        sys.path[0] + '/Databases/whog')
    
    # organize the results from cdd2cog and write protein COG assignment
    out_format = 'tsv' if args.tsv else 'excel'
    timed_message('Retrieving COG categories from COGs.')
    cogblast = organize_cdd_blast(args.output + '/results/rps-blast_cog.txt')
    write_table(cogblast[['qseqid'] + cogblast.columns.tolist()[:-1]],
                         args.output + '/protein2cog', 
                         out_format = out_format)
    
    # quantify COG categories
    timed_message('Quantifying COG categories.')
    del cogblast['qseqid']
    cogblast = cogblast.groupby(cogblast.columns.tolist()).size().reset_index().rename(columns={0:'count'})
    write_table(cogblast, 
                args.output + '/cog_quantification', 
                out_format = out_format)
    timed_message('COG categories quantification is available at {}.'.format(
            args.output + '/cog_quantification'))
    
    # represent that quantification in krona plot
    timed_message('Creating Krona plot representation.')
    write_table(cogblast[['count'] + cogblast.columns.tolist()[:-1]], 
                         args.output + '/cog_quantification.tsv',
                         header = False,
                         out_format = 'tsv')
    create_krona_plot(args.output + '/cog_quantification.tsv')
            
if __name__ == '__main__':
    main()