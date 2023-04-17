# reCOGnizer

A tool for domain-based annotation with databases from the [Conserved Domains Database](https://www.ncbi.nlm.nih.gov/Structure/cdd/cdd.shtml).

* [Features](https://github.com/iquasere/reCOGnizer#features)
* [Installing reCOGnizer](https://github.com/iquasere/reCOGnizer#installing-recognizer)
* [Annotation with reCOGnizer](https://github.com/iquasere/reCOGnizer#annotation-with-recognizer)
* [Output](https://github.com/iquasere/reCOGnizer#output)
* [Other parameters](https://github.com/iquasere/reCOGnizer#other-parameters)
* [Referencing reCOGnizer](https://github.com/iquasere/reCOGnizer#referencing-recognizer)


## Features

reCOGnizer performs domain-based annotation with RPS-BLAST and databases from CDD as reference.
* Reference databases currently implemented: CDD, NCBIfam, Pfam, TIGRFAM, Protein Clusters, SMART, COG and KOG.
* reCOGnizer builds split versions of these databases with which RPS-BLAST can run in multithread, significantly increasing the speed of annotation.
* After domain assignment to proteins, reCOGnizer converts CDD IDs to the IDs of the respective DBs, and obtains domain descriptions available at CDD.
* Further information is retrieved depending on the database in question:
    * NCBIfam, Pfam, TIGRFAM and Protein Clusters annotations are complemented with taxonomic classifications and EC numbers
    * SMART annotations are complemented with SMART descriptions
    * COG and KOG annotations are complemented with COG categories and EC numbers and KEGG Orthologs (for COG)

A detailed representation of reCOGnizer's workflow is presented in Fig. 1.

![ScreenShot](recognizer_workflow.jpg)
Fig. 1. Workflow of reGOGnizer, which includes the pre-analysis step of constructing databases, the domain-based annotation of inputted protein sequences, the interconversion of CDD IDs to other databases IDs, the ID mapping through several databases for obtaining information about the obtained IDs, and the output of information into TSV, Excel and HTML reports.

## Installing reCOGnizer

To install reCOGnizer, simply run: ```conda install -c conda-forge -c bioconda recognizer```

### From source

reCOGnizer can also be installed from source by running:
```
git clone https://github.com/iquasere/reCOGnizer.git
sudo reCOGnizer/install.bash
```

To test that reCOGnizer was correctly installed, run: ```recognizer.py -v```

## Annotation with reCOGnizer

The simplest way to run reCOGnizer is to just specify the fasta filename and an output directory - though even the output directory is not mandatory. It is recommended that a "resources" directory is specified to store the databases that reCOGnizer requires.
```
recognizer.py -f input_file.fasta -o output_directory
```

## Output

reCOGnizer takes a FASTA file as input and produces two main outputs into the output directory:
* ```reCOGnizer_results.tsv```, a table with the annotations for each protein
* ```cog_quantification``` and respective Krona representation (Fig. 2), which describes the functional landscape of the proteins in the input file

[![Image Alt Text](krona_plot.png)](https://iquasere.github.io/reCOGnizer)

Fig. 2. Krona plot with the quantification of COGs identified in the simulated dataset used to test [MOSCA](https://github.com/iquasere/MOSCA) and reCOGnizer.

## Using previously gathered taxonomic information

reCOGnizer can make use of taxonomic information by filtering Markov Models for the specific taxa of interest. 
This can be done by providing a file with the taxonomic information of the proteins.


## Other parameters

```
  -h, --help            show this help message and exit
  -f FILE, --file FILE  Fasta file with protein sequences for annotation
  -t THREADS, --threads THREADS
                        Number of threads for reCOGnizer to use [max available - 1]
  --evalue EVALUE       Maximum e-value to report annotations for [1e-2]
  --pident PIDENT       [DEPRECATED] Minimum pident to report annotations for [0]
  -o OUTPUT, --output OUTPUT
                        Output directory [reCOGnizer_results]
  -dr, --download-resources
                        If resources for reCOGnizer are not available at "resources_directory" [false]
  -rd RESOURCES_DIRECTORY, --resources-directory RESOURCES_DIRECTORY
                        Output directory for storing databases and other resources [~/recognizer_resources]
  -dbs DATABASES, --databases DATABASES
                        Databases to include in functional annotation (comma-separated) [all available]
  -db DATABASE, --database DATABASE
                        Basename of database for annotation. If multiple databases, use comma separated list (db1,db2,db3)
  --custom-database     If database was NOT produced by reCOGnizer
  -mts MAX_TARGET_SEQS, --max-target-seqs MAX_TARGET_SEQS
                        Number of maximum identifications for each protein [1]
  --keep-spaces         BLAST ignores sequences IDs after the first space. This option changes all spaces to underscores to keep the full IDs.
  --no-output-sequences
                        Protein sequences from the FASTA input will be stored in their own column.
  --no-blast-info       Information from the alignment will be stored in their own columns.
  --quiet               Don't output download information, used mainly for CI.
  -sd, --skip-downloaded
                        Skip download of resources detected as already downloaded.
  --keep-intermediates  Keep intermediate annotation files generated in reCOGnizer's workflow, i.e., ASN, RPSBPROC and BLAST reports and split FASTA inputs.
  -v, --version         show program's version number and exit

Taxonomy Arguments:
  --tax-file TAX_FILE   File with taxonomic identification of proteins inputted (TSV). Must have one line per query, query name on first column, taxid on second.
  --protein-id-col PROTEIN_ID_COL
                        Name of column with protein headers as in supplied FASTA file [qseqid]
  --tax-col TAX_COL     Name of column with tax IDs of proteins [Taxonomic identifier (SPECIES)]
  --species-taxids      If tax col contains Tax IDs of species (required for running COG taxonomic)
```

## Referencing reCOGnizer

If you use reCOGnizer, please cite its [publication](https://www.sciencedirect.com/science/article/pii/S2001037022001179).
