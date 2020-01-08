"""
reCOGnizer - a tool for functional annotation with COGs

By João Sequeira

Nov 2019
"""

import pandas as pd
from time import gmtime, strftime
import argparse, shutil, os, multiprocessing, glob, subprocess, pathlib

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
                        default = 'Databases')
    parser.add_argument("-db", "--database", type = str,                        # TODO - Add option to use already existing database by inputing its preffix as argument
                        help = "Basename of COG database for annotation")
    parser.add_argument("-seqs", "--max-target-seqs", type = str,
                        help="""Number of maximum identifications for each protein.
                        Default is 1.""", default = "1")

    args = parser.parse_args()
    
    args.output = args.output.rstrip('/'); args.output_databases = args.output_databases.rstrip('/'); 

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
    fasta: str - name of a fasta file of proteins to be annotated
    output: str - filename of rps_blast to be created
    cog: str - COG blast DB basename from ftp://ftp.ncbi.nlm.nih.gov/pub/mmdb/cdd/little_endian/Cog_LE.tar.gz
    threads: str - number of threads to use
    max_target_seqs: str - number of identifications to obtain for each protein
Output: 
    annotated file with CDD IDs
'''
def run_rpsblast(fasta, output, cog, threads = '0', max_target_seqs = '1'):
    bashCommand = 'rpsblast -query {} -db "{}" -out {} -outfmt 6 -num_threads {} -max_target_seqs {}'.format(
            fasta, cog, output, threads, max_target_seqs)
    open('Databases/command.bash', 'w').write(bashCommand + '\n') # subprocess was not handling well running this command, so an intermediate file is written with the command to run
    print(bashCommand)
    run_command('bash Databases/command.bash', print_command = False)
    os.remove('Databases/command.bash')
    
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
    os.rename('results', output + '/results')
    
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
def organize_cdd_blast(cogblast, fun = 'Databases/fun.txt'):
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
def create_split_cog_db(smp_directory, output, threads = '6', step = None):
    database_reporter = '/'.join(output.split('/')[:-1]) + '/databases.txt'
    dbs = (open(database_reporter).read().split('\n') if
    os.path.isfile(database_reporter) else list())
    if threads in dbs:
        print('Already built COG database for [' + threads + '] threads.')
    else:
        print('Generating COG databases for [' + threads + '] threads.')
        smp_list = glob.glob(smp_directory + '/COG*.smp')
        if step is None:
            step = round(len(smp_list) / float(threads))
        i = 0
        pn_files = list()
        output += '_' + threads + '_'
        while i + step < len(smp_list):
            pn_files.append(output + str(int(i / step)) + '.pn')
            open(output + str(int(i / step)) + '.pn', 'w').write('\n'.join(smp_list[i:i + step]))
            i += step
        open(output + str(int(i / step)) + '.pn', 'w').write('\n'.join(smp_list[i:len(smp_list)]))
        pn_files.append(output + str(int(i / step)) + '.pn')
        for file in pn_files:
            run_command('makeprofiledb -in {0} -title {1} -out {1}'.format(     # -title and -out options are defaulted as input file name to -in argument; -dbtype default is 'rps'
                    file, file.split('.pn')[0]))
        open(database_reporter,'w').write('\n'.join(dbs + [threads]))
        
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
    run_command('ktImportText {} -o {}'.format(tsv, output))
        
def main():
    
    # get arguments
    args = get_arguments()
    
    # create database if it doesn't exit
    timed_message('Checking if database exists for {} threads.'.format(args.threads))
    create_split_cog_db('Databases', args.output_databases + '/COG', args.threads)
    
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
        'Databases/cddid.tbl', 'Databases/fun.txt', 'Databases/whog')
    
    # organize the results from cdd2cog and write protein COG assignment
    timed_message('Retrieving COG categories from COGs.')
    cogblast = organize_cdd_blast(args.output + '/results/rps-blast_cog.txt')
    cogblast.to_excel(args.output + '/protein2cog.xlsx')
    
    # quantify COG categories
    timed_message('Quantifying COG categories.')
    del cogblast['qseqid']
    cogblast = cogblast.groupby(cogblast.columns.tolist()).size().reset_index().rename(columns={0:'count'})
    cogblast.to_csv(args.output + '/cog_quantification.tsv', sep = '\t', index = False)
    timed_message('COG categories quantification is available at {}.'.format(
            args.output + '/cog_quantification.tsv'))
    
    # represent that quantification in krona plot
    timed_message('Creating Krona plot representation.')
    create_krona_plot(args.output + '/cog_quantification.tsv')
            
if __name__ == '__main__':
    main()
